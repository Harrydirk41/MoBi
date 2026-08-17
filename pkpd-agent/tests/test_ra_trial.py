"""RA QSP virtual-trial summarizing/scoring and loop-tool registration.

All synthetic (no MATLAB/SimBiology): the pure-Python layer that turns a Vpop run
into arm response rates, builds the sequential-switch dose spec, and scores a
predicted flagship against a held-out target.
"""

import unittest

from pkpd_agent.engines import osp_ra_trial as R
from pkpd_agent.tools.registry import ToolRegistry
from pkpd_agent.tools.ra_trial_loop_tools import register_ra_trial_loop_tools


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


class _FakeSession:
    def __init__(self):
        self._d = {}

    def get(self, k, default=None):
        return self._d.get(k, default)

    def put(self, k, v):
        self._d[k] = v


if __name__ == "__main__":
    unittest.main()
