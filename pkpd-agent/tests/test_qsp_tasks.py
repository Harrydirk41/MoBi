"""Model-agnostic downstream-task engine, config loading, and loop-tool registration.

All synthetic (no MATLAB/SimBiology): the pure-Python layer that turns a Vpop run
into arm response rates via the config's column map, builds dose/override/sample
strings, scores predictions, and registers the general loop tools off a QSPTaskConfig
loaded from projects/vantage_ra/tasks.json.
"""

import unittest

from pkpd_agent.engines import qsp_tasks as T
from pkpd_agent.engines import qsp_config as C
from pkpd_agent.tools.registry import ToolRegistry


CFG = C.get("vantage_ra")
COLS = CFG.run_columns


def _run(**cols):
    return {"columns": cols}


class TestConfigLoading(unittest.TestCase):
    def test_loads_and_validates(self):
        self.assertEqual(CFG.validate(), [])
        self.assertEqual(len(CFG.readout_states), len(C.READOUT_ROLES))

    def test_alias_and_catalogs(self):
        self.assertEqual(C.get("ra").name, CFG.name)
        self.assertIn("TCZ", CFG.drugs)
        self.assertIn("F_IL6", CFG.design_targets)
        self.assertIn("KD_TCZ", CFG.fit_params)

    def test_unknown_project_raises(self):
        with self.assertRaises(KeyError):
            C.get("no_such_model")

    def test_clinical_weeks_are_int_keyed(self):
        self.assertEqual(CFG.trial_target("TCZ", 24, "raw"),
                         {"ACR20": 45.0, "ACR50": 29.0, "ACR70": 13.9, "rem": 38.0})


class TestClinicalReference(unittest.TestCase):
    def test_raw_and_pcorr(self):
        raw = CFG.trial_target("TCZ", 24, "raw")
        self.assertEqual(raw, {"ACR20": 45.0, "ACR50": 29.0, "ACR70": 13.9, "rem": 38.0})
        pcorr = CFG.trial_target("TCZ", 24, "pcorr")
        self.assertEqual(pcorr, {"ACR20": 20.0, "ACR50": 19.0, "ACR70": 12.0, "rem": 36.0})

    def test_pcorr_floored(self):
        mtx = CFG.trial_target("MTX", 12, "pcorr")
        self.assertTrue(all(v >= 0 for v in mtx.values()))

    def test_unknown_returns_none(self):
        self.assertIsNone(CFG.trial_target("TCZ", 99, "raw"))
        self.assertIsNone(CFG.trial_target("XYZ", 24, "raw"))


class TestSummarizeRun(unittest.TestCase):
    def test_first_and_second_line_rates(self):
        res = _run(
            patient=[1, 2, 3, 4],
            DAS28_base=[5.0, 5.2, 4.8, 6.0], DAS28_read=[3.0, 4.9, 4.7, 3.2],
            ACR20=[1, 0, 0, 1], ACR50=[1, 0, 0, 0], ACR70=[0, 0, 0, 0],
            Rem=[1, 0, 0, 1], MTX_NonResp=[0, 1, 1, 1],
            TCZ_ACR20=[0, 1, 0, 1], TCZ_ACR50=[0, 1, 0, 0],
            TCZ_ACR70=[0, 0, 0, 0], TCZ_Rem=[0, 1, 0, 0])
        s = CFG.summarize_run(res)
        self.assertEqual(s["first_line"]["n"], 4)
        self.assertEqual(s["first_line"]["ACR20"], 50.0)
        self.assertEqual(s["second_line"]["n_subgroup"], 3)
        self.assertEqual(s["second_line"]["ACR20"], 66.7)
        self.assertEqual(s["second_line"]["ACR70"], 0.0)
        self.assertEqual(s["severity"]["baseline_mean"], 5.25)

    def test_empty_second_line(self):
        res = _run(patient=[1, 2], ACR20=[1, 1], ACR50=[1, 0], ACR70=[0, 0],
                   Rem=[1, 0], MTX_NonResp=[0, 0], TCZ_ACR20=[0, 0], TCZ_ACR50=[0, 0],
                   TCZ_ACR70=[0, 0], TCZ_Rem=[0, 0])
        s = CFG.summarize_run(res)
        self.assertEqual(s["second_line"]["n_subgroup"], 0)
        self.assertIsNone(s["second_line"]["ACR20"])

    def test_nan_flags_ignored(self):
        nan = float("nan")
        res = _run(patient=[1, 2, 3], ACR20=[1, 0, nan], ACR50=[0, 0, 0],
                   ACR70=[0, 0, 0], Rem=[0, 0, 0], MTX_NonResp=[1, 1, nan],
                   TCZ_ACR20=[1, nan, 0], TCZ_ACR50=[0, 0, 0],
                   TCZ_ACR70=[0, 0, 0], TCZ_Rem=[0, 0, 0])
        s = CFG.summarize_run(res)
        self.assertEqual(s["second_line"]["n_subgroup"], 2)
        self.assertEqual(s["second_line"]["ACR20"], 100.0)


class TestDoseSpecAndScore(unittest.TestCase):
    def test_dose_spec_switch(self):
        self.assertEqual(
            T.build_dose_spec(["MTX_15mg_Q1W_SC_t200"], ["TCZ8mgkg_Q4W_IV_t200"], 285),
            "MTX_15mg_Q1W_SC_t200;TCZ8mgkg_Q4W_IV_t200@285")

    def test_dose_spec_scale(self):
        self.assertEqual(
            T.build_dose_spec(["A"], ["B"], 285, 0.5), "A;B*0.5@285")
        self.assertEqual(T.build_dose_spec(["A"], ["B"], 285, 1.0), "A;B@285")

    def test_min_dose(self):
        self.assertTrue(T.score_min_dose({"ACR20": 40.0}, 0.7, 35.0)["target_met"])
        self.assertFalse(T.score_min_dose({"ACR20": 20.0}, 0.5, 35.0)["target_met"])

    def test_flagship_mae(self):
        sc = T.score_flagship({"ACR20": 45.0, "ACR50": 24.0}, {"ACR20": 45.0, "ACR50": 24.0})
        self.assertEqual(sc["mae_pp"], 0.0)
        sc2 = T.score_flagship({"ACR20": 40.0, "ACR50": None},
                               {"ACR20": 50.0, "ACR50": 24.0})
        self.assertEqual(sc2["n_endpoints"], 1)
        self.assertEqual(sc2["mae_pp"], 10.0)


class TestCalibration(unittest.TestCase):
    def test_override_spec(self):
        self.assertEqual(T.build_override_spec({"KD_TCZ": 2.5e-12}), "KD_TCZ=2.5e-12")
        self.assertEqual(T.build_override_spec({}), "")

    def test_score_fit_logfold(self):
        sc = T.score_fit({"ACR20": 40.0}, {"ACR20": 45.0},
                         {"KD_TCZ": 2.5e-11}, {"KD_TCZ": 2.5e-12})
        self.assertEqual(sc["acr_mae_pp"], 5.0)
        self.assertEqual(sc["parameters"]["KD_TCZ"]["log10_fold_from_ref"], 1.0)

    def test_numeric_fit_log(self):
        import math
        target = 1.65e-12
        res = T.numeric_fit_1d(lambda v: abs(math.log10(v) - math.log10(target)) * 10.0,
                               1e-13, 1e-9, log=True, max_evals=25)
        self.assertLess(abs(math.log10(res["fitted"]) - math.log10(target)), 0.1)


class TestVpop(unittest.TestCase):
    def test_sample_spec(self):
        s = T.build_sample_spec({"F_TNFa": (0.1, 50, "log"), "kg": (2e4, 4e6)})
        self.assertEqual(s, "F_TNFa,0.1,50,log;kg,20000,4e+06,lin")

    def test_score_vpop(self):
        gs = T.score_vpop([5.0, 5.2, 4.8, 6.0, 3.5, 5.5, 4.5, 5.1], CFG.vpop_target)
        self.assertEqual(gs["yield_pct"], 100.0)
        self.assertLess(gs["distribution_distance"], 1.0)

    def test_select_to_moments(self):
        import random
        random.seed(1)
        pool = [random.uniform(2.0, 8.5) for _ in range(150)]
        sel = T.select_to_moments(pool, CFG.vpop_target)
        self.assertTrue(sel["ok"])
        self.assertLess(abs(sel["weighted_mean"] - 5.12), 0.6)

    def test_multi_anchor_selection(self):
        # a cohort with severity + two therapy-response flags; select weights to match
        # THREE anchors at once (a moment + two response rates)
        import random
        random.seed(3)
        cohort = []
        for _ in range(200):
            sev = random.uniform(2.0, 8.0)
            cohort.append({"severity": sev,
                           "MTX_ACR20": 1 if sev < 4.5 else 0,   # milder respond to MTX
                           "TCZ_ACR20": 1 if sev < 6.5 else 0})  # more respond to TCZ
        anchors = [{"key": "severity", "mean": 5.0, "sd": 1.3},
                   {"key": "MTX_ACR20", "target": 35.0},
                   {"key": "TCZ_ACR20", "target": 70.0}]
        r = T.select_multi_anchor(cohort, anchors)
        self.assertTrue(r["ok"])
        by = {a["key"]: a["achieved"] for a in r["anchors"]}
        # each anchor is matched within a few units (a feasible weighting exists)
        self.assertLess(abs(by["severity"] - 5.0), 0.6)
        self.assertLess(abs(by["MTX_ACR20"] - 35.0), 8.0)
        self.assertLess(abs(by["TCZ_ACR20"] - 70.0), 8.0)
        self.assertGreater(r["effective_sample_size"], 5)

    def test_multi_anchor_needs_candidates(self):
        self.assertFalse(T.select_multi_anchor([], [{"key": "x", "target": 1}])["ok"])

    def test_drivers_have_spans(self):
        for name, p in CFG.vpop_drivers.items():
            self.assertEqual(len(p["span"]), 2, name)
            self.assertLess(p["span"][0], p["span"][1], name)


class TestValidationHelpers(unittest.TestCase):
    def test_ir_mask_and_intersection(self):
        mtx = _run(patient=[1, 2, 3, 4, 5], ACR50=[1, 0, 0, 0, 1],
                   DAS28_read=[2.5, 4.0, 5.0, 3.5, 3.0])
        ada = _run(patient=[1, 2, 3, 4, 5], ACR50=[0, 0, 1, 0, 0],
                   DAS28_read=[3.0, 4.5, 2.8, 3.6, 3.1])
        mm = CFG.ir_mask(mtx, response_key="ACR50")
        am = CFG.ir_mask(ada, response_key="ACR50")
        self.assertTrue(mm[2] and mm[3] and mm[4])
        self.assertFalse(mm[1] or mm[5])
        self.assertEqual({p for p in mm if mm[p] and am.get(p)}, {2, 4})

    def test_response_in_subgroup(self):
        tcz = _run(patient=[1, 2, 3, 4, 5], ACR20=[1, 1, 1, 0, 1],
                   ACR50=[0, 1, 0, 0, 1], ACR70=[0, 0, 0, 0, 0])
        resp = CFG.response_in_subgroup(tcz, {2, 4}, roles=["ACR20", "ACR50", "ACR70"])
        self.assertEqual(resp["n"], 2)
        self.assertEqual(resp["ACR20"], 50.0)
        self.assertEqual(resp["ACR70"], 0.0)


class _FakeSession:
    def __init__(self):
        self._d = {}

    def get(self, k, default=None):
        return self._d.get(k, default)

    def put(self, k, v):
        self._d[k] = v


class TestLoopToolRegistration(unittest.TestCase):
    def test_trial_tools(self):
        from pkpd_agent.tools.qsp_trial_loop_tools import register_qsp_trial_loop_tools
        reg = ToolRegistry()
        register_qsp_trial_loop_tools(reg, None, {"cfg": CFG, "sb": None, "vpop": "V"})
        for n in ("trial_inspect", "trial_run", "trial_finalize"):
            self.assertIn(n, reg)
        res = reg.dispatch("trial_inspect", {}, _FakeSession())
        self.assertTrue(res.ok)
        self.assertIn("drug_formulary", res.data)
        self.assertFalse(reg.dispatch("trial_run", {}, _FakeSession()).ok)

    def test_fit_tools(self):
        from pkpd_agent.tools.qsp_fit_loop_tools import register_qsp_fit_loop_tools
        reg = ToolRegistry()
        register_qsp_fit_loop_tools(reg, None, {
            "cfg": CFG, "sb": None, "vpop": "V", "arm": "X",
            "target": {"ACR20": 45.0}, "enable_optimize": True})
        for n in ("fit_inspect", "fit_try", "fit_optimize", "fit_finalize"):
            self.assertIn(n, reg)
        res = reg.dispatch("fit_inspect", {}, _FakeSession())
        self.assertEqual(res.data["parameters_to_fit"][0]["name"], "KD_TCZ")

    def test_fit_optimize_gated(self):
        from pkpd_agent.tools.qsp_fit_loop_tools import register_qsp_fit_loop_tools
        reg = ToolRegistry()
        register_qsp_fit_loop_tools(reg, None, {
            "cfg": CFG, "sb": None, "vpop": "V", "arm": "X",
            "target": {"ACR20": 45.0}, "enable_optimize": False})
        self.assertNotIn("fit_optimize", reg)

    def test_vpop_tools(self):
        from pkpd_agent.tools.qsp_vpop_loop_tools import register_qsp_vpop_loop_tools
        reg = ToolRegistry()
        register_qsp_vpop_loop_tools(reg, None, {"cfg": CFG, "sb": None})
        for n in ("vpop_inspect", "vpop_sample", "vpop_select", "vpop_finalize"):
            self.assertIn(n, reg)
        res = reg.dispatch("vpop_inspect", {}, _FakeSession())
        self.assertTrue(res.data["disease_driver_parameters"])
        self.assertIn("mean", res.data["clinical_target"])

    def test_design_tools(self):
        from pkpd_agent.tools.qsp_design_loop_tools import register_qsp_design_loop_tools
        reg = ToolRegistry()
        register_qsp_design_loop_tools(reg, None, {
            "cfg": CFG, "sb": None, "sbproj": "p", "vpop": "v"})
        for n in ("design_inspect", "design_try", "design_finalize"):
            self.assertIn(n, reg)
        res = reg.dispatch("design_inspect", {}, _FakeSession())
        self.assertEqual(len(res.data["targetable_pathways"]), len(CFG.design_targets))
        self.assertFalse(reg.dispatch("design_try", {}, _FakeSession()).ok)
        self.assertFalse(reg.dispatch(
            "design_try", {"target": "F_XYZ", "efficacy": 0.8}, _FakeSession()).ok)

    def test_enrichment_elite_and_refit(self):
        from pkpd_agent.engines import qsp_tasks
        # 4 candidates; F over a wide log span; only in-band responders are elite
        cols = {"sev_base": [2.0, 4.0, 5.0, 9.0],   # #0 too mild, #3 too severe
                "MTX":      [1,   1,   0,   1],       # #0,1,3 respond; #2 doesn't
                "F":        [0.1, 8.0, 50.0, 90.0]}
        elite = qsp_tasks.elite_mask(cols, ["MTX"], [3.2, 8.0])
        self.assertEqual(elite, [1])                 # only #1 is in-band AND responds
        bounds = {"F": [0.01, 100.0, "log"]}
        nb = qsp_tasks.refit_bounds(cols, bounds, [1, 2])  # elites span F 8..50
        lo, hi, scale = nb["F"]
        self.assertEqual(scale, "log")
        self.assertLess(lo, 8.0)                      # padded below the elite min
        self.assertGreater(hi, 50.0)                  # padded above the elite max
        self.assertGreaterEqual(lo, 0.01)             # clipped to original
        self.assertLessEqual(hi, 100.0)
        # too few elites -> bounds unchanged
        self.assertEqual(qsp_tasks.refit_bounds(cols, bounds, [1])["F"], [0.01, 100.0, "log"])

    def test_concat_columns(self):
        from pkpd_agent.engines import qsp_tasks
        a = {"x": [1, 2], "y": [3, 4]}
        b = {"x": [5], "y": [6]}
        c = qsp_tasks.concat_columns(a, b)
        self.assertEqual(c["x"], [1, 2, 5])
        self.assertEqual(c["y"], [3, 4, 6])
        self.assertEqual(qsp_tasks.concat_columns({}, b)["x"], [5])

    def test_realize_vpop(self):
        from pkpd_agent.engines import qsp_tasks
        # patient 3 carries almost all the weight -> the realized Vpop is enriched in it
        w = [0.0, 0.01, 0.01, 0.98]
        vp = qsp_tasks.realize_vpop(w, size=200, seed=1)
        self.assertEqual(vp["size"], 200)
        self.assertTrue(vp["unique"] <= 4)
        frac3 = sum(1 for i in vp["indices"] if i == 3) / len(vp["indices"])
        self.assertGreater(frac3, 0.8)          # high-weight patient dominates
        self.assertNotIn(0, vp["indices"])      # zero-weight patient never drawn
        # degenerate input -> empty, no crash
        self.assertEqual(qsp_tasks.realize_vpop([0, 0], 10)["size"], 0)

    def test_filter_columns_to_band(self):
        from pkpd_agent.engines import qsp_tasks
        cols = {"sev_base": [2.0, 3.5, 5.0, 8.5, 6.0], "MTX": [1, 0, 1, 1, 0],
                "sample": [0, 1, 2, 3, 4]}
        f = qsp_tasks.filter_columns_to_band(cols, "sev_base", [3.2, 8.0])
        # 2.0 (too low) and 8.5 (too high) dropped; 3.5/5.0/6.0 kept, in lockstep
        self.assertEqual(f["n_kept"], 3)
        self.assertEqual(f["n_total"], 5)
        self.assertEqual(f["columns"]["sev_base"], [3.5, 5.0, 6.0])
        self.assertEqual(f["columns"]["MTX"], [0, 1, 0])
        self.assertEqual(f["columns"]["sample"], [1, 2, 4])

    def test_validate_tools(self):
        from pkpd_agent.tools.qsp_validate_loop_tools import \
            register_qsp_validate_loop_tools
        reg = ToolRegistry()
        register_qsp_validate_loop_tools(reg, None, {"cfg": CFG, "sb": None, "vpop": "v"})
        for n in ("validate_inspect", "validate_run", "validate_finalize"):
            self.assertIn(n, reg)
        res = reg.dispatch("validate_inspect", {}, _FakeSession())
        self.assertIn("available_arms", res.data)
        self.assertIn("comparator", res.data)


if __name__ == "__main__":
    unittest.main()
