"""Tests for the evaluation-report assembly (no PK-Sim, no API key).

Locks the honesty-critical behaviour: parameters are labelled by what actually
happened in the run (estimated vs fixed), fixed parameters appear as their own
rows, and no narrative section is ever shipped empty.
"""

import os
import unittest

from pkpd_agent.report import (assemble, llm_narrative, deterministic_narrative,
                               ReportData, _estimable_leftovers)
from pkpd_agent.config import AgentConfig
from pkpd_agent.state import ModelingSession, Observation

SNAP = os.path.join(os.path.dirname(__file__), "..", "..",
                    "OSP-PBPK-Model-Library", "Alfentanil", "benchmark",
                    "Alfentanil-Model.blanked.json")

INPUT = {
    "compound": "Alfentanil", "objective": "Build a model",
    "background": {"description": "bg", "literature_facts": ["CYP3A4"]},
    "given_data": {"clinical_observed_data": [
        {"dataset": "S1", "study": "S1", "route": "IV",
         "time_h": [0.5, 1, 2], "conc_mg_L": [0.1, 0.05, 0.02]}]},
}


def _session():
    s = ModelingSession(goal="g")
    s.transcript.append(Observation(call_id="1", tool="osp_optimize", ok=True,
                                    content={"message": "optimized",
                                             "optimized": {"Lipophilicity": 1.45}}))
    return s


@unittest.skipUnless(os.path.exists(SNAP), "benchmark snapshot not present")
class TestReportAssembly(unittest.TestCase):
    def setUp(self):
        # fraction unbound was moved INTO the fit; GFR fraction held fixed at 0
        self.best_edits = {
            "parameters": {"Lipophilicity": 1.4532, "Intrinsic clearance": 0.5623,
                           "Fraction unbound (plasma, reference value)": 0.1162},
            "fix": {"GFR fraction": 0.0},
            "calculation_methods": {"partition": "Rodgers and Rowland"}}
        self.d = assemble(_session(), AgentConfig(mock=True), None, INPUT, SNAP,
                          self.best_edits, run_models=False)

    def test_estimated_params_labelled_estimated(self):
        roles = {p["name"]: p["role"] for p in self.d.parameters}
        # fraction unbound was fitted this run -> must NOT be labelled "measured"
        self.assertEqual(roles["Fraction unbound (plasma, reference value)"],
                         "estimated")
        self.assertEqual(roles["Lipophilicity"], "estimated")

    def test_fixed_param_shown_as_fixed_row(self):
        roles = {p["name"]: p["role"] for p in self.d.parameters}
        self.assertEqual(roles.get("GFR fraction"), "fixed")

    def test_ode_section_present_and_model_aware(self):
        ode = self.d.odes
        self.assertTrue(ode.get("equations"))
        kidney = [f for n, f, _ in ode["equations"] if "Kidney" in n]
        note = [nt for n, _, nt in ode["equations"] if "Kidney" in n][0]
        self.assertIn("ZERO", note)   # GFR fraction = 0 -> renal sink off


class TestEstimableLeftovers(unittest.TestCase):
    """The parameter table must surface estimate-tier knobs left at a default
    (e.g. lipophilicity) WITHOUT re-listing derived siblings of a fitted
    specific-clearance (Specific clearance / Enzyme concentration), which are
    coupled outputs, not independent unfitted knobs."""

    COMP = {
        "Lipophilicity": [{"Parameters": [{"Name": "Lipophilicity", "Value": 0.016}]}],
        "Parameters": [{"Name": "Cl", "Value": 1.0},
                       {"Name": "Molecular weight", "Value": 408.9,
                        "ValueOrigin": {"Source": "Publication"}}],
        "Processes": [{"Molecule": "Hepatic-CYP", "Parameters": [
            {"Name": "Enzyme concentration", "Value": 1.0},
            {"Name": "Specific clearance", "Value": 0.0},
            {"Name": "CLspec/[Enzyme]", "Value": 0.223}]}],
    }

    def test_derived_siblings_hidden_real_knob_surfaced(self):
        rows = _estimable_leftovers(
            self.COMP,
            listed_names={"CLspec/[Enzyme]@Hepatic-CYP", "Lipophilicity"},
            ref_params={})
        names = {r["name"] for r in rows}
        # derived siblings of the fitted CLspec must NOT appear
        self.assertNotIn("Specific clearance@Hepatic-CYP", names)
        self.assertNotIn("Enzyme concentration@Hepatic-CYP", names)
        # a genuine unfitted knob (compound-level Cl) is surfaced
        self.assertIn("Cl", names)
        self.assertEqual(next(r["role"] for r in rows if r["name"] == "Cl"),
                         "held-at-default")
        # a literature constant (MW) is NOT surfaced (tier=constant, not estimate)
        self.assertNotIn("Molecular weight", names)


class TestNarrativeNeverEmpty(unittest.TestCase):
    def _data(self, gmfe, opt_ok):
        return ReportData(
            title="t", objective="o", background="", known_biology=[],
            data_overview={}, nca_rows=[],
            structure={"calculation_methods": [], "processes": ["CYP3A4"]},
            parameters=[{"name": "Lipophilicity", "value": 1.9, "unit": "",
                         "role": "estimated", "plausible_range": [-2, 7]}],
            fit={"gmfe": gmfe}, reference={}, profiles=[], narrative={},
            trajectory=[],
            diagnostics={"optimization_succeeded": opt_ok, "params_fitted": True,
                         "fit_verdict": "good" if gmfe and gmfe <= 1.5 else "poor",
                         "gmfe": gmfe})

    def test_deterministic_conclusion_non_empty(self):
        for gmfe, ok in [(1.45, True), (3.8, True), (None, False)]:
            nar = deterministic_narrative(self._data(gmfe, ok))
            for key in ("model_choice", "parameter_rationale", "conclusion"):
                self.assertTrue(str(nar.get(key) or "").strip(),
                                f"{key} empty for gmfe={gmfe}")

    def test_llm_narrative_falls_back_without_api(self):
        # no anthropic client available in tests -> deterministic, all sections set
        nar = llm_narrative(self._data(1.45, True), AgentConfig(mock=True))
        self.assertTrue(str(nar.get("conclusion") or "").strip())


class TestLiteratureAnchoredRationale(unittest.TestCase):
    """The rationale must COMPARE fitted values to literature anchors and flag
    departures / weak identifiability, not rubber-stamp everything as plausible."""

    def _data(self, params):
        from pkpd_agent.report import ReportData
        return ReportData(
            title="t", objective="o", background="", known_biology=[],
            data_overview={}, nca_rows=[],
            structure={"calculation_methods": ["Rodgers and Rowland"],
                       "processes": ["CYP3A4"]},
            parameters=params, fit={"gmfe": 1.48}, reference={}, profiles=[],
            narrative={}, trajectory=[],
            diagnostics={"optimization_succeeded": True, "params_fitted": True,
                         "fit_verdict": "good", "gmfe": 1.48},
            literature=[
                {"parameter": "logD (pH 7.4)", "reported_values": [2.1, 2.2]},
                {"parameter": "Fraction unbound in plasma",
                 "reported_range_percent": [8.6, 12.0]}])

    def test_flags_lipophilicity_departure(self):
        from pkpd_agent.report import deterministic_narrative
        d = self._data([{"name": "Lipophilicity", "value": 1.45, "unit": "",
                         "role": "estimated", "plausible_range": [-2, 7]}])
        pr = deterministic_narrative(d)["parameter_rationale"]
        self.assertIn("DEPARTS", pr)      # 1.45 vs literature 2.1-2.2

    def test_confirms_fu_consistent(self):
        from pkpd_agent.report import deterministic_narrative
        d = self._data([{"name": "Fraction unbound (plasma, reference value)",
                         "value": 0.116, "unit": "", "role": "estimated",
                         "plausible_range": [0.001, 1.0]}])
        pr = deterministic_narrative(d)["parameter_rationale"]
        self.assertIn("consistent with the literature", pr)

    def test_weak_identifiability_from_sensitivity_not_name(self):
        # driven by the optimizer's sensitivity number, not by the parameter name
        from pkpd_agent.report import deterministic_narrative
        d = self._data([{"name": "Permeability", "value": 1.4e-4, "unit": "cm/min",
                         "role": "estimated", "plausible_range": [1e-6, 1.0],
                         "sensitivity": 0.03}])
        pr = deterministic_narrative(d)["parameter_rationale"]
        self.assertIn("weakly identifiable", pr.lower())

    def test_high_sensitivity_reads_as_well_constrained(self):
        from pkpd_agent.report import deterministic_narrative
        d = self._data([{"name": "Intrinsic clearance", "value": 0.5, "unit": "l/min",
                         "role": "estimated", "plausible_range": [1e-3, 5.0],
                         "sensitivity": 0.9}])
        pr = deterministic_narrative(d)["parameter_rationale"]
        self.assertIn("well constrained", pr.lower())


class TestNarrativeScaffoldSanitizer(unittest.TestCase):
    """A model that echoes the write_sections tool-call format as text (dumping
    all three sections + tags into model_choice) must be recovered, not shipped
    with raw <parameter name="..."></invoke> tags in the report."""

    def test_recovers_leaked_sections(self):
        from pkpd_agent.report import _sanitize_sections
        leaked = {"model_choice": (
            "Structure rationale here.</model_choice>\n"
            "<parameter name=\"parameter_rationale\">Param rationale here.\n"
            "<parameter name=\"conclusion\">Conclusion here.</conclusion>\n</invoke>"),
            "parameter_rationale": "", "conclusion": ""}
        out = _sanitize_sections(leaked)
        self.assertIn("Structure rationale", out["model_choice"])
        self.assertIn("Param rationale", out["parameter_rationale"])
        self.assertIn("Conclusion here", out["conclusion"])
        for v in out.values():
            for tag in ("<parameter", "</invoke>", "</model_choice>", "</conclusion>"):
                self.assertNotIn(tag, v)

    def test_clean_input_untouched(self):
        from pkpd_agent.report import _sanitize_sections
        clean = {"model_choice": "A.", "parameter_rationale": "B.", "conclusion": "C."}
        self.assertEqual(_sanitize_sections(clean),
                         {"model_choice": "A.", "parameter_rationale": "B.",
                          "conclusion": "C."})


if __name__ == "__main__":
    unittest.main()


class TestGroundTruthComparison(unittest.TestCase):
    """The comparison section reports structure + per-parameter diff vs the
    ground truth and grades each on distance AND identifiability."""

    A = {"calculation_methods": ["Cellular partition coefficient method - Rodgers and Rowland",
                                  "Cellular permeability - PK-Sim Standard"],
         "processes": ["CYP3A4"]}

    def _cmp(self, params, ref=None):
        from pkpd_agent.report import _comparison_analysis
        return _comparison_analysis(self.A, ref or self.A, params,
                                    {"gmfe": 1.47}, {"gmfe": 1.45})

    def test_structure_match_detected(self):
        c = self._cmp([{"name": "Lipophilicity", "value": 1.9, "reference": 1.85,
                        "sensitivity": 1.0, "role": "estimated"}])
        self.assertTrue(all(s["match"] for s in c["structure"]))
        self.assertIn("MATCHES", c["summary"])

    def test_structure_difference_detected(self):
        from pkpd_agent.report import _comparison_analysis
        ref = {"calculation_methods": ["Cellular partition coefficient method - Schmitt"],
               "processes": ["CYP3A4", "GlomerularFiltration"]}
        c = _comparison_analysis(self.A, ref, [], {"gmfe": 1.5}, {"gmfe": 1.4})
        self.assertFalse(all(s["match"] for s in c["structure"]))
        self.assertIn("DIFFERS", c["summary"])

    def test_estimated_recovered_well_is_good(self):
        c = self._cmp([{"name": "CLint", "value": 0.5, "reference": 0.53,
                        "sensitivity": 0.6, "role": "estimated"}])
        self.assertEqual(c["parameters"][0]["grade"], "good")

    def test_estimated_far_and_influential_is_a_miss(self):
        c = self._cmp([{"name": "CLint", "value": 5.0, "reference": 0.5,
                        "sensitivity": 0.9, "role": "estimated"}])
        self.assertEqual(c["parameters"][0]["grade"], "bad")
        self.assertIn("miss", c["parameters"][0]["verdict"].lower())

    def test_estimated_far_but_weak_is_not_a_miss(self):
        c = self._cmp([{"name": "Permeability", "value": 0.001, "reference": 0.05,
                        "sensitivity": 0.03, "role": "estimated"}])
        self.assertEqual(c["parameters"][0]["grade"], "soft")

    def test_fixed_far_off_is_not_a_miss(self):
        # a FIXED parameter differing from the reference is a prior choice, not a fit error
        c = self._cmp([{"name": "Permeability", "value": 0.0015, "reference": 0.0069,
                        "sensitivity": None, "role": "fixed"}])
        self.assertEqual(c["parameters"][0]["grade"], "soft")
        self.assertIn("prior choice", c["parameters"][0]["verdict"])

    def test_collinear_far_off_is_not_a_miss(self):
        # two influential clearances that trade off: far from truth AND high
        # one-at-a-time sensitivity, but collinear -> only the combination is
        # identifiable, so the split is NOT graded a real miss.
        c = self._cmp([{"name": "CLspec/[Enzyme]@UGT2B7", "value": 0.08,
                        "reference": 0.0066, "sensitivity": 1.0,
                        "collinearity": 0.94, "collinear_with": "CLspec/[Enzyme]@UGT1A9",
                        "role": "estimated"}])
        self.assertEqual(c["parameters"][0]["grade"], "soft")
        self.assertIn("collinear", c["parameters"][0]["verdict"].lower())
        self.assertIn("collinear", c["summary"].lower())

    def test_held_at_default_far_off_is_not_a_miss(self):
        # an estimate-tier parameter left at a bare default is a prior choice,
        # not a fitting error (it was never estimated).
        c = self._cmp([{"name": "Lipophilicity", "value": 2.1, "reference": 0.7,
                        "sensitivity": None, "role": "held-at-default"}])
        self.assertEqual(c["parameters"][0]["grade"], "soft")
        self.assertIn("not estimated", c["parameters"][0]["verdict"])

    def test_no_reference_yields_empty(self):
        from pkpd_agent.report import _comparison_analysis
        c = _comparison_analysis(self.A, None,
                                 [{"name": "x", "value": 1.0, "reference": None}],
                                 {"gmfe": 1.5}, {})
        self.assertEqual(c, {})
