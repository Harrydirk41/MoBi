"""osp_sweep_methods: the deterministic partition x permeability grid sweep that re-fits physchem
under each method and adopts the best by GMFE. This is the general fix for the failure where the
agent tries a method at frozen (wrong) physchem, sees no change, and wrongly abandons it (the
Vancomycin case: it even tried Schmitt but with logP=0 and gave up). The engine run is mocked, so no
PK-Sim is needed; the test asserts the sweep covers the whole grid, picks the best method by DATA,
and hardcodes no answer. Also covers the GFR-fraction plausibility guard."""

import unittest

from pkpd_agent.tools.registry import ToolRegistry
from pkpd_agent.tools.osp_loop_tools import register_osp_loop_tools
from pkpd_agent.config import AgentConfig
from pkpd_agent.state import ModelingSession
from pkpd_agent.engines import osp_optimize as OO
from pkpd_agent.engines import osp_score
from pkpd_agent.engines.snapshot_edit import PARTITION_METHODS, PERMEABILITY_METHODS


class TestSweepMethods(unittest.TestCase):
    def setUp(self):
        self.seen = []
        # a mock optimizer where ONE method combo (Schmitt / Charge dependent Schmitt) fits best -
        # mirroring the Vancomycin ground truth. The sweep must FIND it from the GMFE alone.
        def fake_run(cli, snap, observed, estimate, fix=None, structure=None, max_evals=30, **kw):
            cm = (structure or {}).get("calculation_methods", {})
            part, perm = cm.get("partition"), cm.get("permeability")
            self.seen.append((part, perm))
            good = (part == "Schmitt" and perm == "Charge dependent Schmitt")
            gmfe = 1.40 if good else 2.35
            return {"ok": True, "optimized": {k: sum(v) / 2 for k, v in estimate.items()},
                    "fit": {"gmfe": gmfe}, "by_route": {}, "worst_datasets": [],
                    "params_at_bound": [], "sensitivity": {}, "n_evals": max_evals,
                    "fit_simulations": ["s1"]}
        self._orig = OO.run_optimization
        OO.run_optimization = fake_run
        self.reg = ToolRegistry()
        register_osp_loop_tools(self.reg, AgentConfig(mock=False),
                                {"cli": object(), "snapshot_path": "x", "observed": [], "input": {}})
        self.sweep = self.reg.get("osp_sweep_methods").handler

    def tearDown(self):
        OO.run_optimization = self._orig

    def test_sweeps_full_grid_and_picks_best_by_data(self):
        sess = ModelingSession(goal="g")
        r = self.sweep({"estimate": {"Lipophilicity": [-2, 3],
                                     "Fraction unbound (plasma, reference value)": [0.3, 0.8]}}, sess)
        self.assertTrue(r.ok)
        # covered EVERY partition x permeability combo (the enumerable 15-grid), not a sample
        self.assertEqual(set(self.seen), {(p, q) for p in PARTITION_METHODS
                                          for q in PERMEABILITY_METHODS})
        self.assertEqual(len(self.seen), len(PARTITION_METHODS) * len(PERMEABILITY_METHODS))
        # picked the best-fitting method purely from GMFE - no answer was passed in
        best = r.data["best"]
        self.assertEqual((best["partition"], best["permeability"]),
                         ("Schmitt", "Charge dependent Schmitt"))
        self.assertEqual(best["gmfe"], 1.40)
        # adopted into the session so the report re-runs the winning structure
        self.assertEqual(sess.get("osp_best_gmfe"), 1.40)
        self.assertEqual(sess.get("osp_best_edits")["calculation_methods"]["partition"], "Schmitt")

    def test_refits_physchem_under_each_method(self):
        # every combo re-runs the optimizer with the physchem in 'estimate' (not frozen)
        self.sweep({"estimate": {"Lipophilicity": [-2, 3]}}, ModelingSession(goal="g"))
        self.assertEqual(len(self.seen), 15)   # 15 independent fits, one per method combo

    def test_requires_estimate(self):
        r = self.sweep({}, ModelingSession(goal="g"))
        self.assertFalse(r.ok)
        self.assertIn("estimate", r.message)


class TestGFRGuard(unittest.TestCase):
    def test_gfr_out_of_range_flagged(self):
        flags = osp_score.plausibility([{"parameter": "GFR fraction", "value": 1.4}])
        self.assertTrue(any("GFR fraction" in f["message"] for f in flags))

    def test_gfr_pushed_low_flagged(self):
        flags = osp_score.plausibility([{"parameter": "GFR fraction", "value": 0.3}])
        self.assertTrue(any("reabsorption" in f["message"] for f in flags))

    def test_gfr_physiological_ok(self):
        self.assertEqual(osp_score.plausibility([{"parameter": "GFR fraction", "value": 1.0}]), [])


if __name__ == "__main__":
    unittest.main()
