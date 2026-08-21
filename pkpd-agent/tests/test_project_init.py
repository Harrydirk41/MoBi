"""Project validation + conversational config builder (traditional-modeler onboarding).

validate_project catches the silent foot-guns (a role that is not a real model
parameter, an unfilled stub) as plain-English messages. build_tasks merges an LLM reply
(stubbed here) onto the drafter's pools. No API / MATLAB - stub call_fn and dict fixtures.
"""

import unittest

from pkpd_agent.engines import project_validate as PV
from pkpd_agent.engines import llm_config_build as LC


_NETWORK = {
    "name": "Test QSP",
    "species": [{"name": s} for s in ("ACR20", "MTX_NonResp", "TCZ_ACR20", "DAS28_CRP")],
    "parameters": [{"name": p} for p in ("F_IL6", "F_TNFa", "KD_TCZ", "kg_FLS_Baseline")],
}


def _good_tasks():
    return {
        "readout_states": ["ACR20", "MTX_NonResp", "TCZ_ACR20", "DAS28_CRP"],
        "run_columns": {"patient": "patient", "first_line": {"ACR20": "ACR20"},
                        "subgroup_flag": "MTX_NonResp", "second_line": {"ACR20": "TCZ_ACR20"},
                        "severity": {"baseline": "DAS28_base", "readout": "DAS28_read"}},
        "timeline": {"baseline_day": 200, "first_line_readout_day": 284,
                     "second_line_readout_day": 600},
        "drugs": {"TCZ": {"doses": ["TCZ8mg"]}},
        "vpop_drivers": {"F_IL6": {}, "F_TNFa": {}},
        "vpop_target": {"mean": 5.0, "sd": 1.0, "band": [3, 8]},
        "fit_params": {"KD_TCZ": {}},
        "design_targets": {"F_IL6": {}},
        "clinical_trials": {"TCZ": {}}, "refractory_target": {"ACR20": 50},
    }


_SPEC = {"readout_targets": ["DAS28_CRP"], "gsa_top": ["F_IL6"]}


class TestValidate(unittest.TestCase):
    def test_clean_config_passes(self):
        r = PV.validate_project(_good_tasks(), _SPEC, _NETWORK)
        self.assertEqual(r["errors"], [])
        self.assertEqual(r["warnings"], [])

    def test_bad_param_role_is_an_error(self):
        t = _good_tasks()
        t["vpop_drivers"] = {"F_IL6": {}, "F_NOT_A_PARAM": {}}
        r = PV.validate_project(t, _SPEC, _NETWORK)
        self.assertTrue(any("F_NOT_A_PARAM" in m for m in r["errors"]))

    def test_missing_run_columns_is_error(self):
        t = _good_tasks()
        del t["run_columns"]
        r = PV.validate_project(t, _SPEC, _NETWORK)
        self.assertTrue(any("run_columns" in m for m in r["errors"]))

    def test_stub_external_field_is_warning(self):
        t = _good_tasks()
        t["clinical_trials"] = "TODO: author must supply"
        r = PV.validate_project(t, _SPEC, _NETWORK)
        self.assertTrue(any("clinical_trials" in m for m in r["warnings"]))

    def test_readout_state_not_a_species_warns(self):
        t = _good_tasks()
        t["readout_states"] = t["readout_states"] + ["NotASpecies"]
        r = PV.validate_project(t, _SPEC, _NETWORK)
        self.assertTrue(any("NotASpecies" in m for m in r["warnings"]))

    def test_gsa_top_bad_param_warns(self):
        r = PV.validate_project(_good_tasks(), {"readout_targets": ["DAS28_CRP"],
                                                "gsa_top": ["ghost_param"]}, _NETWORK)
        self.assertTrue(any("ghost_param" in m for m in r["warnings"]))

    def test_format_report_ok(self):
        self.assertIn("OK", PV.format_report({"errors": [], "warnings": []}))


class TestBuildTasks(unittest.TestCase):
    def _draft(self):
        return {"readout_states": ["ACR20", "DAS28_CRP"],
                "vpop_drivers": {"F_IL6": {}, "F_TNFa": {}, "F_junk": {}},
                "design_targets": {"F_IL6": {}}, "fit_params": {"KD_TCZ": {}},
                "_review": "note", "readout_roles": {}, "drugs": "TODO"}

    def test_llm_output_wins_draft_fills_gaps(self):
        # LLM prunes vpop to F_IL6 and fills drugs; draft supplies readout_states
        call = lambda s, u: ('{"vpop_drivers":{"F_IL6":{"nominal":8}},'
                             '"drugs":{"TCZ":{"doses":["TCZ8mg"]}}}')
        out = LC.build_tasks(_NETWORK, "use IL-6 only", self._draft(), call)
        self.assertEqual(list(out["vpop_drivers"]), ["F_IL6"])   # LLM pruned
        self.assertIn("TCZ", out["drugs"])                        # LLM filled
        self.assertEqual(out["readout_states"], ["ACR20", "DAS28_CRP"])  # from draft
        # internal drafter markers are stripped
        self.assertNotIn("_review", out)
        self.assertNotIn("readout_roles", out)

    def test_bad_llm_reply_falls_back_to_draft(self):
        call = lambda s, u: "not json at all"
        out = LC.build_tasks(_NETWORK, "desc", self._draft(), call)
        self.assertIn("F_IL6", out["vpop_drivers"])              # draft preserved
        self.assertNotIn("_review", out)

    def test_build_spec_from_structure(self):
        spec = LC.build_spec(_NETWORK, "desc",
                             {"classification": {"readout": ["DAS28_CRP"]}})
        self.assertEqual(spec["readout_targets"], ["DAS28_CRP"])
        self.assertEqual(spec["name"], "Test QSP")


if __name__ == "__main__":
    unittest.main()
