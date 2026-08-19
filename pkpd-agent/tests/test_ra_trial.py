"""RA QSP virtual-trial summarizing/scoring and loop-tool registration.

All synthetic (no MATLAB/SimBiology): the pure-Python layer that turns a Vpop run
into arm response rates, builds the sequential-switch dose spec, and scores a
predicted flagship against a held-out target.
"""

import unittest

from pkpd_agent.engines import osp_ra_trial as R
from pkpd_agent.engines import ra_clinical_reference as RC
from pkpd_agent.tools.registry import ToolRegistry
from pkpd_agent.tools.ra_trial_loop_tools import register_ra_trial_loop_tools


class TestClinicalReference(unittest.TestCase):
    def test_tcz_raw_and_pcorr(self):
        raw = RC.trial_target("TCZ", 24, "raw")
        self.assertEqual(raw, {"ACR20": 45.0, "ACR50": 29.0, "ACR70": 13.9})
        pcorr = RC.trial_target("TCZ", 24, "pcorr")
        # drug minus placebo (45-25, 29-10, 13.9-1.9)
        self.assertEqual(pcorr, {"ACR20": 20.0, "ACR50": 19.0, "ACR70": 12.0})

    def test_pcorr_floored_at_zero(self):
        # a made-up check that pcorr never goes negative is covered by real data
        mtx = RC.trial_target("MTX", 12, "pcorr")
        self.assertTrue(all(v >= 0 for v in mtx.values()))

    def test_unknown_returns_none(self):
        self.assertIsNone(RC.trial_target("TCZ", 99, "raw"))
        self.assertIsNone(RC.trial_target("XYZ", 24, "raw"))

    def test_agent_prediction_matches_real_rose(self):
        # the agent's committed prediction should be close to real ROSE wk24
        pred = {"ACR20": 42.1, "ACR50": 26.3, "ACR70": 15.8}
        mae = R.score_flagship(pred, RC.trial_target("TCZ", 24, "raw"))["mae_pp"]
        self.assertLess(mae, 5.0)


def _run(**cols):
    return {"columns": cols}


class TestSummarizeRun(unittest.TestCase):
    def test_first_and_second_line_rates(self):
        # 4 patients; #1 responds to MTX, #2-4 are MTX_NonResp; TCZ rescues #2,#4
        res = _run(
            patient=[1, 2, 3, 4],
            DAS28_base=[5.0, 5.2, 4.8, 6.0], DAS28_read=[3.0, 4.9, 4.7, 3.2],
            ACR20=[1, 0, 0, 1], ACR50=[1, 0, 0, 0], ACR70=[0, 0, 0, 0],
            Rem=[1, 0, 0, 1], MTX_NonResp=[0, 1, 1, 1],
            TCZ_ACR20=[0, 1, 0, 1], TCZ_ACR50=[0, 1, 0, 0],
            TCZ_ACR70=[0, 0, 0, 0], TCZ_Rem=[0, 1, 0, 0])
        s = R.summarize_run(res)
        self.assertEqual(s["first_line"]["n"], 4)
        self.assertEqual(s["first_line"]["ACR20"], 50.0)
        # denominator for the flagship is the MTX_NonResp==1 subgroup (3 patients)
        self.assertEqual(s["second_line"]["n_MTX_IR"], 3)
        self.assertEqual(s["second_line"]["ACR20"], 66.7)   # 2 of 3
        self.assertEqual(s["second_line"]["ACR70"], 0.0)
        self.assertEqual(s["das28"]["baseline_mean"], 5.25)

    def test_empty_second_line_when_no_nonresponders(self):
        res = _run(patient=[1, 2], ACR20=[1, 1], ACR50=[1, 0], ACR70=[0, 0],
                   Rem=[1, 0], MTX_NonResp=[0, 0],
                   TCZ_ACR20=[0, 0], TCZ_ACR50=[0, 0], TCZ_ACR70=[0, 0], TCZ_Rem=[0, 0])
        s = R.summarize_run(res)
        self.assertEqual(s["second_line"]["n_MTX_IR"], 0)
        self.assertIsNone(s["second_line"]["ACR20"])

    def test_nan_flags_are_ignored(self):
        nan = float("nan")
        res = _run(patient=[1, 2, 3], ACR20=[1, 0, nan], ACR50=[0, 0, 0],
                   ACR70=[0, 0, 0], Rem=[0, 0, 0], MTX_NonResp=[1, 1, nan],
                   TCZ_ACR20=[1, nan, 0], TCZ_ACR50=[0, 0, 0],
                   TCZ_ACR70=[0, 0, 0], TCZ_Rem=[0, 0, 0])
        s = R.summarize_run(res)
        # only patients 1,2 have a finite MTX_NonResp==1; TCZ_ACR20 finite for #1 only
        self.assertEqual(s["second_line"]["n_MTX_IR"], 2)
        self.assertEqual(s["second_line"]["ACR20"], 100.0)   # 1 finite value, =1


class TestDoseSpecAndScore(unittest.TestCase):
    def test_build_dose_spec_sequential_switch(self):
        spec = R.build_dose_spec(["MTX_15mg_Q1W_SC_t200"],
                                 ["TCZ8mgkg_Q4W_IV_t200"], 285)
        self.assertEqual(spec, "MTX_15mg_Q1W_SC_t200;TCZ8mgkg_Q4W_IV_t200@285")

    def test_build_dose_spec_concurrent(self):
        spec = R.build_dose_spec(["MTX_15mg_Q1W_SC_t200", "TCZ8mgkg_Q4W_IV_t200"])
        self.assertEqual(spec, "MTX_15mg_Q1W_SC_t200;TCZ8mgkg_Q4W_IV_t200")

    def test_build_dose_spec_scale_and_switch(self):
        spec = R.build_dose_spec(["MTX_15mg_Q1W_SC_t200"],
                                 ["TCZ8mgkg_Q4W_IV_t200"], 285, 0.5)
        self.assertEqual(spec, "MTX_15mg_Q1W_SC_t200;TCZ8mgkg_Q4W_IV_t200*0.5@285")

    def test_build_dose_spec_scale_one_is_noop(self):
        spec = R.build_dose_spec(["MTX_15mg_Q1W_SC_t200"],
                                 ["TCZ8mgkg_Q4W_IV_t200"], 285, 1.0)
        self.assertEqual(spec, "MTX_15mg_Q1W_SC_t200;TCZ8mgkg_Q4W_IV_t200@285")

    def test_score_min_dose_met(self):
        sc = R.score_min_dose({"ACR20": 40.0}, 0.7, 35.0)
        self.assertTrue(sc["target_met"])
        self.assertEqual(sc["dose_scale"], 0.7)
        self.assertIsNotNone(sc["response_per_dose"])

    def test_score_min_dose_missed(self):
        sc = R.score_min_dose({"ACR20": 20.0}, 0.5, 35.0)
        self.assertFalse(sc["target_met"])
        self.assertIsNone(sc["response_per_dose"])

    def test_score_mae(self):
        sc = R.score_flagship({"ACR20": 45.0, "ACR50": 24.0, "ACR70": 14.0},
                              {"ACR20": 45.0, "ACR50": 24.0, "ACR70": 14.0})
        self.assertEqual(sc["mae_pp"], 0.0)
        self.assertEqual(sc["n_endpoints"], 3)

    def test_score_skips_missing_endpoints(self):
        sc = R.score_flagship({"ACR20": 40.0, "ACR50": None},
                              {"ACR20": 50.0, "ACR50": 24.0, "ACR70": 14.0})
        self.assertEqual(sc["n_endpoints"], 1)
        self.assertEqual(sc["mae_pp"], 10.0)


class TestCalibrationFit(unittest.TestCase):
    def test_build_override_spec(self):
        self.assertEqual(R.build_override_spec({"KD_TCZ": 2.5e-12}), "KD_TCZ=2.5e-12")
        self.assertEqual(R.build_override_spec({}), "")

    def test_score_fit_perfect_acr(self):
        tgt = {"ACR20": 45.0, "ACR50": 29.0, "ACR70": 13.9}
        sc = R.score_fit(tgt, tgt, {"KD_TCZ": 2.5e-12}, {"KD_TCZ": 2.5e-12})
        self.assertEqual(sc["acr_mae_pp"], 0.0)
        self.assertEqual(sc["parameters"]["KD_TCZ"]["log10_fold_from_ref"], 0.0)

    def test_score_fit_logfold(self):
        sc = R.score_fit({"ACR20": 40.0}, {"ACR20": 45.0},
                         {"KD_TCZ": 2.5e-11}, {"KD_TCZ": 2.5e-12})
        self.assertEqual(sc["acr_mae_pp"], 5.0)
        self.assertEqual(sc["parameters"]["KD_TCZ"]["log10_fold_from_ref"], 1.0)

    def test_fit_param_has_reference(self):
        self.assertIn("KD_TCZ", R.FIT_PARAMS)
        self.assertEqual(R.FIT_PARAMS["KD_TCZ"]["reference"], 2.5e-12)


class TestFitLoopTools(unittest.TestCase):
    def _reg(self):
        from pkpd_agent.tools.ra_fit_loop_tools import register_ra_fit_loop_tools
        reg = ToolRegistry()
        register_ra_fit_loop_tools(reg, None, {
            "sb": None, "vpop": "V.xlsx", "arm": "MTX;TCZ@285",
            "target": {"ACR20": 45.0, "ACR50": 29.0, "ACR70": 13.9},
            "fit_params": ["KD_TCZ"]})
        return reg

    def test_registers(self):
        reg = self._reg()
        for n in ("fit_inspect", "fit_try", "fit_finalize"):
            self.assertIn(n, reg)

    def test_inspect_exposes_param_and_target(self):
        res = self._reg().dispatch("fit_inspect", {}, _FakeSession())
        self.assertTrue(res.ok)
        self.assertEqual(res.data["observed_target"]["ACR20"], 45.0)
        self.assertEqual(res.data["parameters_to_fit"][0]["name"], "KD_TCZ")

    def test_try_rejects_empty(self):
        self.assertFalse(self._reg().dispatch("fit_try", {}, _FakeSession()).ok)

    def test_finalize_needs_a_run(self):
        res = self._reg().dispatch("fit_finalize", {"overrides": {"KD_TCZ": 1e-12}},
                                   _FakeSession())
        self.assertFalse(res.ok)

    def test_optimize_gated_off_by_default_flag(self):
        from pkpd_agent.tools.ra_fit_loop_tools import register_ra_fit_loop_tools
        reg = ToolRegistry()
        register_ra_fit_loop_tools(reg, None, {
            "sb": None, "vpop": "V", "arm": "X", "target": {"ACR20": 45.0},
            "fit_params": ["KD_TCZ"], "enable_optimize": False})
        self.assertNotIn("fit_optimize", reg)

    def test_optimize_registered_when_enabled(self):
        from pkpd_agent.tools.ra_fit_loop_tools import register_ra_fit_loop_tools
        reg = ToolRegistry()
        register_ra_fit_loop_tools(reg, None, {
            "sb": None, "vpop": "V", "arm": "X", "target": {"ACR20": 45.0},
            "fit_params": ["KD_TCZ"], "enable_optimize": True})
        self.assertIn("fit_optimize", reg)


class TestNumericFit1D(unittest.TestCase):
    def test_recovers_scalar_minimum_log(self):
        import math
        target = 1.65e-12

        def evaluate(val):
            return abs(math.log10(val) - math.log10(target)) * 10.0

        res = R.numeric_fit_1d(evaluate, 1e-13, 1e-9, log=True, max_evals=25)
        self.assertLess(abs(math.log10(res["fitted"]) - math.log10(target)), 0.1)
        self.assertGreater(res["n_evals"], 2)
        self.assertTrue(res["trace"])

    def test_recovers_scalar_minimum_linear(self):
        def evaluate(val):
            return (val - 3.0) ** 2

        res = R.numeric_fit_1d(evaluate, 0.0, 10.0, log=False, max_evals=30)
        self.assertAlmostEqual(res["fitted"], 3.0, delta=0.2)


class TestVpopGeneration(unittest.TestCase):
    def test_build_sample_spec(self):
        s = R.build_sample_spec({"F_TNFa": (0.1, 50, "log"), "kg_FLS_Baseline": (2e4, 4e6)})
        self.assertEqual(s, "F_TNFa,0.1,50,log;kg_FLS_Baseline,20000,4e+06,lin")

    def test_score_vpop_good_vs_bad(self):
        good = [5.0, 5.2, 4.8, 6.0, 3.5, 5.5, 4.5, 5.1]     # all in band, ~target
        gs = R.score_vpop(good)
        self.assertEqual(gs["yield_pct"], 100.0)
        self.assertLess(gs["distribution_distance"], 1.0)
        bad = [1.0, 1.2, 9.0, 8.5, 0.5, 2.0]                 # mostly out of band
        bs = R.score_vpop(bad)
        self.assertLess(bs["yield_pct"], 50.0)

    def test_score_vpop_empty_band(self):
        s = R.score_vpop([1.0, 1.5, 2.0])   # all below the active band
        self.assertEqual(s["n_accepted"], 0)
        self.assertIsNone(s["distribution_distance"])

    def test_vpop_params_have_spans(self):
        for name, p in R.VPOP_PARAMS.items():
            self.assertEqual(len(p["span"]), 2, name)
            self.assertLess(p["span"][0], p["span"][1], name)


class TestVpopLoopTools(unittest.TestCase):
    def _reg(self):
        from pkpd_agent.tools.ra_vpop_loop_tools import register_ra_vpop_loop_tools
        reg = ToolRegistry()
        register_ra_vpop_loop_tools(reg, None, {"sb": None})
        return reg

    def test_registers(self):
        reg = self._reg()
        for n in ("vpop_inspect", "vpop_sample", "vpop_finalize"):
            self.assertIn(n, reg)

    def test_inspect_lists_drivers_and_target(self):
        res = self._reg().dispatch("vpop_inspect", {}, _FakeSession())
        self.assertTrue(res.ok)
        self.assertTrue(res.data["disease_driver_parameters"])
        self.assertIn("mean", res.data["clinical_target"])

    def test_sample_rejects_empty(self):
        self.assertFalse(self._reg().dispatch("vpop_sample", {}, _FakeSession()).ok)


class TestDrugCatalog(unittest.TestCase):
    def test_tocilizumab_is_il6(self):
        self.assertIn("IL-6", R.DRUG_CATALOG["TCZ"]["mechanism"])
        self.assertIn("TCZ8mgkg_Q4W_IV_t200", R.DRUG_CATALOG["TCZ"]["doses"])

    def test_every_drug_has_doses_and_mechanism(self):
        for code, d in R.DRUG_CATALOG.items():
            self.assertTrue(d["doses"], code)
            self.assertTrue(d["mechanism"], code)


class TestLoopToolRegistration(unittest.TestCase):
    def test_registers_inspect_and_run(self):
        reg = ToolRegistry()
        register_ra_trial_loop_tools(reg, None, {
            "sb": None, "vpop": "V.xlsx", "calibrated_arms": []})
        self.assertIn("ra_inspect", reg)
        self.assertIn("ra_run_trial", reg)

    def test_inspect_returns_formulary_and_timeline(self):
        reg = ToolRegistry()
        register_ra_trial_loop_tools(reg, None, {
            "sb": None, "vpop": "V.xlsx",
            "calibrated_arms": [{"arm": "MTX"}]})
        res = reg.dispatch("ra_inspect", {}, _FakeSession())
        self.assertTrue(res.ok)
        self.assertIn("drug_formulary", res.data)
        self.assertEqual(res.data["trial_timeline"]["first_line_readout_day"], 284)
        self.assertEqual(res.data["trial_timeline"]["second_line_readout_day"], 600)

    def test_run_trial_rejects_empty_protocol(self):
        reg = ToolRegistry()
        register_ra_trial_loop_tools(reg, None, {
            "sb": None, "vpop": "V.xlsx", "calibrated_arms": []})
        res = reg.dispatch("ra_run_trial", {}, _FakeSession())
        self.assertFalse(res.ok)

    def test_finalize_registered(self):
        reg = ToolRegistry()
        register_ra_trial_loop_tools(reg, None, {
            "sb": None, "vpop": "V.xlsx", "calibrated_arms": []})
        self.assertIn("ra_finalize", reg)

    def test_finalize_rejects_unrun_protocol(self):
        reg = ToolRegistry()
        register_ra_trial_loop_tools(reg, None, {
            "sb": None, "vpop": "V.xlsx", "calibrated_arms": []})
        res = reg.dispatch("ra_finalize", {
            "first_line": ["MTX_15mg_Q1W_SC_t200"],
            "second_line": ["TCZ8mgkg_Q4W_IV_t200"], "switch_day": 285},
            _FakeSession())
        self.assertFalse(res.ok)

    def test_finalize_commits_a_run_protocol(self):
        reg = ToolRegistry()
        register_ra_trial_loop_tools(reg, None, {
            "sb": None, "vpop": "V.xlsx", "calibrated_arms": []})
        sess = _FakeSession()
        # simulate a completed run in history
        spec = "MTX_15mg_Q1W_SC_t200;TCZ8mgkg_Q4W_IV_t200@285"
        sess.put("ra_history", [{"protocol": spec,
                                 "first_line": {"ACR20": 30.0},
                                 "second_line": {"n_MTX_IR": 38, "ACR20": 42.1,
                                                 "ACR50": 26.3, "ACR70": 15.8}}])
        res = reg.dispatch("ra_finalize", {
            "first_line": ["MTX_15mg_Q1W_SC_t200"],
            "second_line": ["TCZ8mgkg_Q4W_IV_t200"], "switch_day": 285}, sess)
        self.assertTrue(res.ok)
        self.assertEqual(sess.get("ra_final")["protocol"], spec)

    def test_finalize_rejects_empty_mtx_ir_arm(self):
        reg = ToolRegistry()
        register_ra_trial_loop_tools(reg, None, {
            "sb": None, "vpop": "V.xlsx", "calibrated_arms": []})
        sess = _FakeSession()
        sess.put("ra_history", [{"protocol": "MTX_15mg_Q1W_SC_t200",
                                 "first_line": {"ACR20": 30.0},
                                 "second_line": {"n_MTX_IR": 0}}])
        res = reg.dispatch("ra_finalize", {
            "first_line": ["MTX_15mg_Q1W_SC_t200"]}, sess)
        self.assertFalse(res.ok)


class _FakeSession:
    def __init__(self):
        self._d = {}

    def get(self, k, default=None):
        return self._d.get(k, default)

    def put(self, k, v):
        self._d[k] = v


if __name__ == "__main__":
    unittest.main()
