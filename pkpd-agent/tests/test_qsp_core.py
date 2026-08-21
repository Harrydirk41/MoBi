"""Model-agnostic core: Edge, edge scoring, sign scoring, param scoring (no vocab)."""

import unittest

from pkpd_agent.engines import qsp_core as C


class TestEdgeScoring(unittest.TestCase):
    def setUp(self):
        self.truth = [C.Edge("A", 1, "B"), C.Edge("A", 1, "C"), C.Edge("D", -1, "B")]

    def test_perfect(self):
        s = C.score_network(self.truth, self.truth)
        self.assertEqual((s["precision"], s["recall"], s["f1"]), (1.0, 1.0, 1.0))

    def test_sign_aware_vs_topology(self):
        prop = [C.Edge("D", 1, "B")]                 # right pair, wrong sign
        self.assertEqual(C.score_network(prop, self.truth, sign_aware=True)["hit"], 0)
        self.assertEqual(C.score_network(prop, self.truth, sign_aware=False)["hit"], 1)


class TestSignScoring(unittest.TestCase):
    def setUp(self):
        self.truth = [C.Edge("A", 1, "B"), C.Edge("A", 1, "C"), C.Edge("D", -1, "B")]

    def test_majority_baseline(self):
        s = C.score_signs({}, self.truth)
        self.assertAlmostEqual(s["majority_baseline"], 2 / 3, places=2)

    def test_all_positive_ties_majority(self):
        s = C.score_signs({e.pair(): 1 for e in self.truth}, self.truth)
        self.assertEqual(s["correct"], 2)
        self.assertFalse(s["beats_majority"])


class TestParamScoring(unittest.TestCase):
    def test_split_and_baseline(self):
        t = [C.Param("k", "1/day", "x", 1.0), C.Param("KD", "M", "x", 1e-10),
             C.Param("s", "nanogram/(molecule*day)", "x", 1e-9),
             C.Param("f", "dimensionless", "x", 1.5)]
        s = C.score_params({p.name: p.value for p in t}, t)
        self.assertEqual(s["physiological"]["n"], 2)      # 1/day, M
        self.assertEqual(s["model_scaling"]["n"], 1)
        self.assertEqual(s["dimensionless"]["n"], 1)
        self.assertEqual(s["overall"]["median_log10_err"], 0.0)

    def test_geomean_baseline(self):
        t = [C.Param("a", "1/day", "x", 0.1), C.Param("b", "1/day", "x", 10.0)]
        base = C.unit_geomean_baseline(t)
        self.assertAlmostEqual(base["a"], 1.0)

    def test_clean_predictions(self):
        self.assertEqual(C.clean_predictions([{"name": "k", "value": 2}]), {"k": 2.0})
        self.assertEqual(C.clean_predictions({"k": 2}), {"k": 2.0})
        self.assertEqual(C.clean_predictions([{"name": "k"}]), {})


if __name__ == "__main__":
    unittest.main()
