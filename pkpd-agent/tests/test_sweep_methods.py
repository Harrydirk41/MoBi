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

    def test_ties_with_same_physchem_do_NOT_early_stop(self):
        # the mock ties at 2.35 for most methods but returns the SAME optimized value for all - a
        # later method (Schmitt) is better, so the sweep MUST NOT stop early and must find it
        seen = self.seen
        seen.clear()
        r = self.sweep({"estimate": {"Lipophilicity": [-2, 3]}}, ModelingSession(goal="g"))
        self.assertEqual(len(seen), 15)                       # full grid, no early stop
        self.assertEqual(r.data["best"]["gmfe"], 1.40)        # found the later better method


class TestSweepEarlyStopCompensation(unittest.TestCase):
    """Safe early stop: methods tie in GMFE AND a free param shifts across them (compensation) ->
    the degeneracy is structural, stop early. This is the Vancomycin case."""
    def setUp(self):
        self.seen = []

        def fake_run(cli, snap, observed, estimate, fix=None, structure=None, max_evals=30, **kw):
            part = (structure or {}).get("calculation_methods", {}).get("partition")
            self.seen.append(part)
            # every method fits equally (GMFE 1.85) but with a DIFFERENT logP - the free param
            # absorbs the method difference (unidentifiable)
            logp = {"PK-Sim Standard": 0.5, "Rodgers and Rowland": 1.2, "Schmitt": 2.1}.get(part, 3.0)
            return {"ok": True, "optimized": {"Lipophilicity": logp}, "fit": {"gmfe": 1.85},
                    "by_route": {}, "worst_datasets": [], "params_at_bound": [],
                    "sensitivity": {}, "n_evals": max_evals, "fit_simulations": ["s1"]}
        self._orig = OO.run_optimization
        OO.run_optimization = fake_run
        self.reg = ToolRegistry()
        register_osp_loop_tools(self.reg, AgentConfig(mock=False),
                                {"cli": object(), "snapshot_path": "x", "observed": [], "input": {}})
        self.sweep = self.reg.get("osp_sweep_methods").handler

    def tearDown(self):
        OO.run_optimization = self._orig

    def test_stops_early_on_compensation(self):
        self.sweep({"estimate": {"Lipophilicity": [-5.0, 3.0]}}, ModelingSession(goal="g"))
        # stopped after covering >=2 partition methods (not the full 15) because a free param
        # compensated the tie - far fewer than 15 combos
        n_partitions = len(set(self.seen))
        self.assertLessEqual(n_partitions, 3)
        self.assertLess(len(self.seen), 15)


class TestSweepWidening(unittest.TestCase):
    """A FREE (non-given) parameter that rails to a too-tight self-imposed bound is widened to its
    physical range and re-swept - the Vancomycin case: effective Lipophilicity capped at 1.5 (a
    measured logP the input never gave) blocks the fit; the true value 2.23 needs the wider range."""
    def setUp(self):
        self.calls = []

        def fake_run(cli, snap, observed, estimate, fix=None, structure=None, max_evals=30, **kw):
            self.calls.append(dict(estimate))
            hi = estimate["Lipophilicity"][1]
            # good fit only reachable if the upper bound allows ~2.2 (the withheld fitted value)
            if hi >= 2.0:
                return {"ok": True, "optimized": {"Lipophilicity": 2.2}, "fit": {"gmfe": 1.40},
                        "by_route": {}, "worst_datasets": [], "params_at_bound": [],
                        "sensitivity": {}, "n_evals": max_evals, "fit_simulations": ["s1"]}
            # capped too low -> rails to the upper bound, mediocre fit (all methods identical)
            return {"ok": True, "optimized": {"Lipophilicity": hi}, "fit": {"gmfe": 1.93},
                    "by_route": {}, "worst_datasets": [],
                    "params_at_bound": [{"parameter": "Lipophilicity", "value": hi,
                                         "bound": "upper"}],
                    "sensitivity": {}, "n_evals": max_evals, "fit_simulations": ["s1"]}
        self._orig = OO.run_optimization
        OO.run_optimization = fake_run
        self.reg = ToolRegistry()
        # input gives NO lipophilicity -> it is a free fit-target the sweep may widen
        register_osp_loop_tools(self.reg, AgentConfig(mock=False),
                                {"cli": object(), "snapshot_path": "x", "observed": [],
                                 "input": {"given_data": {"literature_physicochemical": [
                                     {"parameter": "Molecular weight", "value": 1449.3}]}}})
        self.sweep = self.reg.get("osp_sweep_methods").handler

    def tearDown(self):
        OO.run_optimization = self._orig

    def test_widens_free_railed_param_and_recovers_fit(self):
        sess = ModelingSession(goal="g")
        r = self.sweep({"estimate": {"Lipophilicity": [-4.0, 1.5]}}, sess)   # too-tight upper bound
        self.assertTrue(r.ok)
        # it re-swept with a widened upper bound (physical ceiling 7.0) after the rail
        self.assertTrue(any(e["Lipophilicity"][1] >= 7.0 for e in self.calls),
                        "should have widened the railed free param to its physical range")
        self.assertEqual(r.data["best"]["gmfe"], 1.40)                       # recovered the good fit
        self.assertEqual(sess.get("osp_best_gmfe"), 1.40)

    def test_given_measurement_is_not_widened(self):
        # if the input GAVE a lipophilicity, the tight bound is respected (not widened past it)
        self.reg = ToolRegistry()
        register_osp_loop_tools(self.reg, AgentConfig(mock=False),
                                {"cli": object(), "snapshot_path": "x", "observed": [],
                                 "input": {"given_data": {"literature_physicochemical": [
                                     {"parameter": "logP", "value": 1.4}]}}})
        sweep = self.reg.get("osp_sweep_methods").handler
        self.calls.clear()
        sweep({"estimate": {"Lipophilicity": [-4.0, 1.5]}}, ModelingSession(goal="g"))
        self.assertFalse(any(e["Lipophilicity"][1] >= 7.0 for e in self.calls),
                         "a GIVEN measurement must not be widened")


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
