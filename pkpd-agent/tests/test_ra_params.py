"""Stage-1 Layer 5 parameter-estimation: key loading, scoring, baseline, loop tools.

All synthetic (no MATLAB, no LLM): the pure-Python layer that loads the ESM2 parameter
key, scores predicted values order-of-magnitude against it, computes the naive baseline,
and drives the estimate/finalize loop.
"""

import math
import unittest

from pkpd_agent.engines import ra_params as P
from pkpd_agent.tools.registry import ToolRegistry
from pkpd_agent.tools.ra_params_loop_tools import register_ra_params_loop_tools


class TestKey(unittest.TestCase):
    def test_loads_130(self):
        t = P.load_truth()
        self.assertGreaterEqual(len(t), 120)
        self.assertTrue(all(p.value != 0 for p in t))

    def test_prompt_view_hides_value(self):
        t = P.load_truth()
        view = P.prompt_view(t)
        self.assertEqual(set(view[0]), {"name", "units", "cell_context"})
        self.assertNotIn("value", view[0])


class TestScoring(unittest.TestCase):
    def setUp(self):
        self.truth = [P.Param("k1", "1/day", "FLS", 0.1),
                      P.Param("k2", "1/day", "FLS", 1.0),
                      P.Param("f1", "dimensionless", "FLS", 1.5),
                      P.Param("f2", "dimensionless", "FLS", 2.0)]

    def test_perfect(self):
        pred = {p.name: p.value for p in self.truth}
        s = P.score_params(pred, self.truth)
        self.assertEqual(s["overall"]["median_log10_err"], 0.0)

    def test_order_of_magnitude_buckets(self):
        # k1 off by 10x, k2 exact, f1/f2 exact
        pred = {"k1": 1.0, "k2": 1.0, "f1": 1.5, "f2": 2.0}
        s = P.score_params(pred, self.truth)
        self.assertEqual(s["overall"]["n"], 4)
        self.assertAlmostEqual(s["overall"]["within_10x"], 1.0)   # 10x counts as <=1 dex

    def test_split_dimensionless_dimensional(self):
        pred = {p.name: p.value for p in self.truth}
        s = P.score_params(pred, self.truth)
        self.assertEqual(s["dimensionless"]["n"], 2)
        self.assertEqual(s["dimensional"]["n"], 2)

    def test_missing_predictions_not_scored(self):
        s = P.score_params({"k1": 0.1}, self.truth)
        self.assertEqual(s["n_scored"], 1)

    def test_baseline_is_computed(self):
        s = P.score_params({p.name: p.value for p in self.truth}, self.truth)
        self.assertIn("n", s["naive_unit_geomean_baseline"])


class TestBaseline(unittest.TestCase):
    def test_unit_geomean(self):
        truth = [P.Param("a", "1/day", "x", 0.1), P.Param("b", "1/day", "x", 10.0)]
        base = P.unit_geomean_baseline(truth)
        # geomean(0.1, 10) = 1.0 for both
        self.assertAlmostEqual(base["a"], 1.0, places=6)
        self.assertAlmostEqual(base["b"], 1.0, places=6)

    def test_physiological_split(self):
        # rates/concentrations are physiological; per-molecule secretion is model-scaling
        t = [P.Param("kcl", "1/day", "x", 1.0), P.Param("kd", "sec-1", "x", 1e-4),
             P.Param("KD", "M", "x", 1e-10),
             P.Param("sec", "nanogram/(molecule*day)", "x", 1e-9),
             P.Param("f", "dimensionless", "x", 1.5)]
        s = P.score_params({p.name: p.value for p in t}, t)
        self.assertEqual(s["physiological"]["n"], 3)     # kcl, kd, KD
        self.assertEqual(s["model_scaling"]["n"], 1)      # sec
        self.assertEqual(s["dimensionless"]["n"], 1)

    def test_real_key_dimensional_is_hard(self):
        # sanity: on the real key the naive baseline should be far from perfect on the
        # dimensional params (that is why the benchmark is meaningful)
        t = P.load_truth()
        s = P.score_params(P.unit_geomean_baseline(t), t)
        self.assertLess(s["dimensional"]["within_3x"], 0.7)


class TestCleanPredictions(unittest.TestCase):
    def test_list_and_dict(self):
        self.assertEqual(P.clean_predictions([{"name": "k", "value": 0.5}]), {"k": 0.5})
        self.assertEqual(P.clean_predictions({"k": 0.5}), {"k": 0.5})

    def test_bad_rows_skipped(self):
        self.assertEqual(P.clean_predictions([{"name": "k"}, {"value": 3}]), {})


class _FakeSession:
    def __init__(self):
        self._d = {}

    def get(self, k, default=None):
        return self._d.get(k, default)

    def put(self, k, v):
        self._d[k] = v


class TestLoopTools(unittest.TestCase):
    def _reg(self, truth):
        reg = ToolRegistry()
        register_ra_params_loop_tools(reg, None, {"truth": truth})
        return reg

    def test_registers(self):
        reg = self._reg([])
        for t in ("param_inspect", "param_estimate", "param_finalize"):
            self.assertIn(t, reg)

    def test_inspect_hides_values(self):
        reg = self._reg([P.Param("k1", "1/day", "FLS", 0.1)])
        res = reg.dispatch("param_inspect", {}, _FakeSession())
        self.assertTrue(res.ok)
        self.assertIn("parameters", res.data)
        self.assertNotIn("value", res.data["parameters"][0])

    def test_estimate_accumulates(self):
        reg = self._reg([P.Param("k1", "1/day", "FLS", 0.1),
                         P.Param("k2", "1/day", "FLS", 1.0)])
        sess = _FakeSession()
        reg.dispatch("param_estimate", {"predictions": [{"name": "k1", "value": 0.2}]}, sess)
        res = reg.dispatch("param_estimate",
                           {"predictions": [{"name": "k2", "value": 1.0}]}, sess)
        self.assertEqual(res.data["n_estimated"], 2)
        self.assertEqual(res.data["n_remaining"], 0)

    def test_finalize_requires_estimates(self):
        reg = self._reg([P.Param("k1", "1/day", "FLS", 0.1)])
        self.assertFalse(reg.dispatch("param_finalize", {}, _FakeSession()).ok)

    def test_finalize_scores(self):
        truth = [P.Param("k1", "1/day", "FLS", 0.1), P.Param("k2", "1/day", "FLS", 1.0)]
        reg = self._reg(truth)
        sess = _FakeSession()
        reg.dispatch("param_estimate",
                     {"predictions": [{"name": "k1", "value": 0.1},
                                      {"name": "k2", "value": 1.0}]}, sess)
        res = reg.dispatch("param_finalize", {}, sess)
        self.assertTrue(res.ok)
        self.assertEqual(sess.get("param_final")["overall"]["median_log10_err"], 0.0)


if __name__ == "__main__":
    unittest.main()
