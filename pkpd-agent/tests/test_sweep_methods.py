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
        # 15 coarse-rank fits cover every combo, + 1 fine refine of the winner
        self.assertEqual(len(self.seen), len(PARTITION_METHODS) * len(PERMEABILITY_METHODS) + 1)
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
        # 15 coarse-rank fits (one per combo) + 1 fine refine of the winner
        self.assertEqual(len(self.seen), 16)

    def test_requires_estimate(self):
        r = self.sweep({}, ModelingSession(goal="g"))
        self.assertFalse(r.ok)
        self.assertIn("estimate", r.message)

    def test_forwards_add_processes_to_every_combo(self):
        """REGRESSION: the sweep once read only structure.processes and dropped
        structure.add_processes, so every method ran on a no-clearance model and the
        clearance estimate froze. It must carry the added enzyme process into each combo."""
        seen_struct = []
        orig = OO.run_optimization

        def capture(cli, snap, observed, estimate, fix=None, structure=None, max_evals=30, **kw):
            seen_struct.append(structure)
            return orig(cli, snap, observed, estimate, fix=fix, structure=structure,
                        max_evals=max_evals, **kw)
        OO.run_optimization = capture
        try:
            self.sweep({"estimate": {"Intrinsic clearance": [0.01, 50]},
                        "structure": {"add_processes": [
                            {"type": "metabolization_first_order", "molecule": "CYP1A2",
                             "parameters": {"Intrinsic clearance": 1.0}}]}},
                       ModelingSession(goal="g"))
        finally:
            OO.run_optimization = orig
        # EVERY optimizer call in the sweep must carry the added process (not just methods)
        self.assertTrue(seen_struct)
        for st in seen_struct:
            ap = (st or {}).get("add_processes") or []
            self.assertTrue(any(p.get("molecule") == "CYP1A2" for p in ap),
                            "sweep dropped add_processes - the method grid ran with no clearance")
            self.assertIn("calculation_methods", st)   # and it still varies the methods

    def test_identical_resweep_is_memoized(self):
        # the agent re-sweeping the SAME grid must NOT re-run PK-Sim (it wasted hours) - the second
        # identical call returns the cached result and steers to osp_optimize
        sess = ModelingSession(goal="g")
        self.sweep({"estimate": {"Lipophilicity": [-2, 3]}}, sess)
        n_after_first = len(self.seen)
        r2 = self.sweep({"estimate": {"Lipophilicity": [-2, 3]}}, sess)   # identical
        self.assertEqual(len(self.seen), n_after_first)          # no additional PK-Sim runs
        self.assertTrue(r2.data.get("cached"))
        self.assertIn("do NOT re-sweep", r2.message)

    def test_ties_with_same_physchem_do_NOT_early_stop(self):
        # the mock ties at 2.35 for most methods but returns the SAME optimized value for all - a
        # later method (Schmitt) is better, so the sweep MUST NOT stop early and must find it
        seen = self.seen
        seen.clear()
        r = self.sweep({"estimate": {"Lipophilicity": [-2, 3]}}, ModelingSession(goal="g"))
        self.assertEqual(len(seen), 16)                       # full 15-grid coarse + 1 refine
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


class TestSweepCoordinateDescentCheapPath(unittest.TestCase):
    """When the partition sweep at the anchor permeability ALREADY fits well (perfusion-limited,
    small molecule), the permeability axis is inert and must be SKIPPED - the whole point of
    coordinate descent. Only the 5 partition methods (at the anchor permeability) are fitted, plus
    the winner's full-fidelity refine - never the other permeabilities."""
    def setUp(self):
        from pkpd_agent.engines.snapshot_edit import PARTITION_METHODS, PERMEABILITY_METHODS
        self.PART, self.PERM = PARTITION_METHODS, PERMEABILITY_METHODS
        self.anchor = PERMEABILITY_METHODS[0]
        self.seen = []

        def fake_run(cli, snap, observed, estimate, fix=None, structure=None, max_evals=30, **kw):
            cm = (structure or {}).get("calculation_methods", {})
            part, perm = cm.get("partition"), cm.get("permeability")
            self.seen.append((part, perm))
            # Rodgers and Rowland already fits great (1.20); others fit fine too. All GOOD (<=1.5),
            # so permeability is inert and should never be tried beyond the anchor.
            gmfe = 1.20 if part == "Rodgers and Rowland" else 1.45
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

    def test_permeability_axis_skipped_when_partition_fit_is_good(self):
        r = self.sweep({"estimate": {"Lipophilicity": [-2, 3]}}, ModelingSession(goal="g"))
        self.assertTrue(r.ok)
        used_perms = {pe for _, pe in self.seen}
        self.assertEqual(used_perms, {self.anchor})            # ONLY the anchor permeability
        # 5 partition coarse fits + 1 winner refine = 6, never the 15+ of a full grid
        self.assertLessEqual(len(self.seen), len(self.PART) + 1)
        self.assertEqual(r.data["best"]["partition"], "Rodgers and Rowland")
        self.assertEqual(r.data["best"]["permeability"], self.anchor)

    def test_full_grid_flag_forces_exhaustive(self):
        self.sweep({"estimate": {"Lipophilicity": [-2, 3]}, "full_grid": True},
                   ModelingSession(goal="g"))
        combos = {(p, q) for p in self.PART for q in self.PERM}
        self.assertTrue(combos <= set(self.seen))              # every combo covered


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


class TestFixGivenPhyschem(unittest.TestCase):
    """A physchem value the input GAVE constrains the fit. A TRUE measurement (fu, pKa)
    is FIXED; but LIPOPHILICITY is an effective distribution parameter, so a given logP
    is turned into a bounded WINDOW around the measured value, not pinned at it (pinning
    it at the measured logP biases Vd and made the alfentanil web-lookup run unfittable)."""
    def test_given_lipophilicity_is_windowed_not_pinned(self):
        from pkpd_agent.tools.osp_loop_tools import _fix_given_physchem
        lit = [{"parameter": "Lipophilicity", "value": 1.4}]
        est, fix, notes = _fix_given_physchem(
            {"Lipophilicity": [0.5, 4.0], "Intrinsic clearance": [0.01, 5.0]}, {}, lit)
        self.assertIn("Lipophilicity", est)             # still estimated (effective param)
        self.assertNotIn("Lipophilicity", fix)          # NOT pinned at the measured value
        lo, hi = est["Lipophilicity"]
        self.assertLessEqual(lo, 1.4)                   # the window brackets the measured 1.4
        self.assertGreaterEqual(hi, 1.4)
        self.assertGreaterEqual(lo, 0.4 - 1e-9)         # ~measured +/- 1 log unit
        self.assertLessEqual(hi, 2.4 + 1e-9)            # (2.4 = 1.4 + 1, tighter than the requested 4.0)
        self.assertIn("Intrinsic clearance", est)       # the genuinely-uncertain one still fitted
        self.assertTrue(notes)

    def test_true_measurement_fraction_unbound_is_still_pinned(self):
        # a REAL measurement (not an effective parameter) is still fixed at its value.
        # (_given_physchem_value currently only surfaces lipophilicity, so this guards the
        # general branch: whatever it surfaces that is NOT lipophilicity gets pinned.)
        from pkpd_agent.tools import osp_loop_tools as T
        orig = T._given_physchem_value
        T._given_physchem_value = lambda base, lit: (0.1 if "unbound" in base.lower() else None)
        try:
            est, fix, notes = T._fix_given_physchem(
                {"Fraction unbound": [0.05, 0.2], "Intrinsic clearance": [0.01, 5.0]}, {}, [])
            self.assertNotIn("Fraction unbound", est)   # a true measurement is fixed
            self.assertEqual(fix["Fraction unbound"], 0.1)
            self.assertIn("Intrinsic clearance", est)
        finally:
            T._given_physchem_value = orig

    def test_not_given_lipophilicity_stays_free(self):
        from pkpd_agent.tools.osp_loop_tools import _fix_given_physchem
        # no lipophilicity in the input -> it remains a free (effective) fit-target
        est, fix, notes = _fix_given_physchem(
            {"Lipophilicity": [-2, 7]}, {}, [{"parameter": "Molecular weight", "value": 300}])
        self.assertIn("Lipophilicity", est)
        self.assertEqual(fix, {})
        self.assertFalse(notes)

    def test_invitro_clearance_is_not_fixed(self):
        from pkpd_agent.tools.osp_loop_tools import _fix_given_physchem
        # an in-vitro kinetic input is NOT a direct physchem measurement - it stays refinable
        est, fix, _ = _fix_given_physchem(
            {"Intrinsic clearance": [0.01, 5.0]}, {},
            [{"parameter": "In vitro CL/recombinant enzyme", "value": 2.0}])
        self.assertIn("Intrinsic clearance", est)
        self.assertEqual(fix, {})


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


class TestSweepScreenVerdict(unittest.TestCase):
    def setUp(self):
        self._orig = OO.run_optimization

    def tearDown(self):
        OO.run_optimization = self._orig

    def _sweep(self, gmfe):
        def fake(cli, snap, observed, estimate, fix=None, structure=None, max_evals=30, **kw):
            return {"ok": True, "optimized": {k: sum(v) / 2 for k, v in estimate.items()},
                    "fit": {"gmfe": gmfe}, "by_route": {}, "worst_datasets": [],
                    "params_at_bound": [], "sensitivity": {}, "n_evals": max_evals,
                    "fit_simulations": ["s1"]}
        OO.run_optimization = fake
        reg = ToolRegistry()
        register_osp_loop_tools(reg, AgentConfig(mock=False),
                                {"cli": object(), "snapshot_path": "x", "observed": [], "input": {}})
        return reg.get("osp_sweep_methods").handler(
            {"estimate": {"Lipophilicity": [-2, 3], "Intrinsic clearance": [0.1, 5]}},
            ModelingSession(goal="g"))

    def test_adequate_when_fit_is_good(self):
        r = self._sweep(1.30)
        self.assertTrue(r.data["screen_adequate"])
        self.assertIsNone(r.data["advice"])

    def test_not_adequate_when_all_poor(self):
        r = self._sweep(3.10)                    # every method fits badly
        self.assertFalse(r.data["screen_adequate"])
        self.assertIn("not the lever here", r.data["advice"])


class TestSweepParallel(unittest.TestCase):
    def setUp(self):
        self._orig = OO.run_optimization

    def tearDown(self):
        OO.run_optimization = self._orig

    def test_parallel_runs_all_combos_and_picks_best(self):
        import threading
        seen, tids, lock = [], set(), threading.Lock()

        def fake(cli, snap, obs, estimate, fix=None, structure=None, max_evals=30, **kw):
            import time
            time.sleep(0.01)
            cm = (structure or {}).get("calculation_methods", {})
            with lock:
                seen.append((cm.get("partition"), cm.get("permeability")))
                tids.add(threading.get_ident())
            good = cm.get("partition") == "Schmitt"
            return {"ok": True, "optimized": {k: sum(v) / 2 for k, v in estimate.items()},
                    "fit": {"gmfe": 1.4 if good else 2.3}, "by_route": {}, "params_at_bound": [],
                    "sensitivity": {}, "n_evals": max_evals, "fit_simulations": ["s1"]}
        OO.run_optimization = fake
        cfg = AgentConfig(mock=False)
        cfg.parallel_jobs = 4
        reg = ToolRegistry()
        register_osp_loop_tools(reg, cfg, {"cli": object(), "snapshot_path": "x",
                                           "observed": [], "input": {}})
        r = reg.get("osp_sweep_methods").handler(
            {"estimate": {"Lipophilicity": [-2, 3], "Intrinsic clearance": [0.1, 5]}},
            ModelingSession(goal="g"))
        self.assertTrue(r.ok)
        self.assertEqual(r.data["best"]["partition"], "Schmitt")
        self.assertGreater(len(tids), 1)              # actually ran on multiple threads


class TestMethodGuidance(unittest.TestCase):
    """Data-shape guidance (not answers): staged IV->PO and saturable identifiability."""
    def _g(self, observed):
        from pkpd_agent.tools.osp_loop_tools import _method_guidance
        return _method_guidance(observed)

    def test_staging_when_iv_and_po_present(self):
        obs = [{"route": "IV", "dose": "1 mg"}, {"route": "PO", "dose": "1 mg"}]
        g = self._g(obs)
        self.assertIn("route_staging", g)
        self.assertIn("IV", g["route_staging"])
        self.assertIn("absorption", g["route_staging"].lower())

    def test_no_staging_when_iv_only(self):
        g = self._g([{"route": "IV", "dose": "1 mg"}, {"route": "IV", "dose": "2 mg"}])
        self.assertNotIn("route_staging", g)          # single route -> nothing to stage
        self.assertEqual(g["dose_levels"], 2)
        self.assertNotIn("saturable_identifiability", g)   # >=2 doses -> identifiable

    def test_saturable_warned_with_single_dose(self):
        g = self._g([{"route": "IV", "dose": "1 mg"}, {"route": "PO", "dose": "1 mg"}])
        self.assertEqual(g["dose_levels"], 1)
        self.assertIn("saturable_identifiability", g)
