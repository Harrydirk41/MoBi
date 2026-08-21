"""LLM tasks.json drafter: role classification, draft assembly, and the regression.

Uses a STUB call_fn (canned JSON) so the drafting/compare logic is tested without an
API. The real LLM is swapped in only in run_llm_draft_tasks. A perfect-drafter test
also regresses against the actual hand-written RA tasks.json via qsp_config.
"""

import unittest

from pkpd_agent.engines import llm_tasks as LT
from pkpd_agent.engines import qsp_config


class TestClassifyParameters(unittest.TestCase):
    def test_roles_parsed_and_filtered(self):
        call = lambda s, u: ('{"disease_drivers":["F_TNFa","F_IL6","kg_FLS_Baseline"],'
                             '"druggable":["F_TNFa","F_IL6"],'
                             '"calibratable":["KD_TCZ", 42]}')
        params = [{"name": n, "units": "x"} for n in
                  ("F_TNFa", "F_IL6", "kg_FLS_Baseline", "KD_TCZ")]
        r = LT.classify_parameters(params, ["rule"], call)
        self.assertEqual(r["disease_drivers"], ["F_TNFa", "F_IL6", "kg_FLS_Baseline"])
        self.assertEqual(r["druggable"], ["F_TNFa", "F_IL6"])
        self.assertEqual(r["calibratable"], ["KD_TCZ"])       # non-str 42 filtered


class TestClassifyReadout(unittest.TestCase):
    def test_pieces_and_assembly_order(self):
        call = lambda s, u: ('{"first_line_flags":["ACR20","ACR50","ACR70","Remission"],'
                             '"subgroup_flag":"MTX_NonResp",'
                             '"second_line_flags":["TCZ20","TCZ50","TCZ70","TCZR"],'
                             '"severity_states":["DAS28_CRP","DAS28_BL"]}')
        ro = LT.classify_readout_states(["ACR20"], ["r"], call)
        states = LT._readout_states(ro)
        self.assertEqual(states, ["ACR20", "ACR50", "ACR70", "Remission", "MTX_NonResp",
                                  "TCZ20", "TCZ50", "TCZ70", "TCZR",
                                  "DAS28_CRP", "DAS28_BL"])

    def test_no_subgroup_flag_skipped(self):
        call = lambda s, u: ('{"first_line_flags":["A"],"subgroup_flag":"",'
                             '"second_line_flags":[],"severity_states":["S"]}')
        ro = LT.classify_readout_states(["A"], ["r"], call)
        self.assertEqual(LT._readout_states(ro), ["A", "S"])


class TestDraftTasks(unittest.TestCase):
    def _call(self, s, u):
        if "role in a two-line trial" in u:
            return ('{"first_line_flags":["ACR20"],"subgroup_flag":"MTX_NonResp",'
                    '"second_line_flags":["TCZ20"],"severity_states":["DAS28_CRP"]}')
        return ('{"disease_drivers":["F_IL6"],"druggable":["F_IL6"],'
                '"calibratable":["KD_TCZ"]}')

    def test_draft_fills_derivable_and_stubs_external(self):
        data = {
            "species": [{"name": n} for n in
                        ("ACR20", "MTX_NonResp", "TCZ20", "DAS28_CRP")],
            "parameters": [{"name": "F_IL6", "units": "dimensionless"},
                           {"name": "KD_TCZ", "units": "M"}],
            "rules": [{"rule": "x=y"}]}
        d = LT.draft_tasks(data, self._call, name="Test model")
        self.assertEqual(d["readout_states"],
                         ["ACR20", "MTX_NonResp", "TCZ20", "DAS28_CRP"])
        self.assertIn("F_IL6", d["vpop_drivers"])
        self.assertIn("F_IL6", d["design_targets"])
        self.assertEqual(d["fit_params"]["KD_TCZ"]["unit"], "M")     # unit carried
        # external fields left as explicit stubs, not invented
        for f in ("drugs", "timeline", "vpop_target", "clinical_trials",
                  "refractory_target"):
            self.assertTrue(str(d[f]).startswith("TODO"), f)


class TestCompareTasks(unittest.TestCase):
    def test_perfect_draft_scores_1(self):
        ref = qsp_config.get("vantage_ra")
        reference = {"readout_states": ref.readout_states, "vpop_drivers": ref.vpop_drivers,
                     "design_targets": ref.design_targets, "fit_params": ref.fit_params}
        # a perfect drafter reproduces exactly the reference role sets
        draft = {
            "readout_states": list(ref.readout_states),
            "vpop_drivers": {k: {} for k in ref.vpop_drivers},
            "design_targets": {k: {} for k in ref.design_targets},
            "fit_params": {k: {} for k in ref.fit_params}}
        cmp = LT.compare_tasks(draft, reference)
        for field in ("readout_states", "vpop_drivers", "design_targets", "fit_params"):
            self.assertEqual(cmp[field]["f1"], 1.0, field)

    def test_partial_draft_reports_misses(self):
        reference = {"vpop_drivers": {"F_TNFa": {}, "F_IL6": {}, "kg_X": {}}}
        draft = {"vpop_drivers": {"F_TNFa": {}, "F_BOGUS": {}}}
        d = LT.compare_tasks(draft, reference)["vpop_drivers"]
        self.assertEqual(d["hit"], 1)
        self.assertIn("F_IL6", d["missed"])
        self.assertIn("F_BOGUS", d["extra"])


if __name__ == "__main__":
    unittest.main()
