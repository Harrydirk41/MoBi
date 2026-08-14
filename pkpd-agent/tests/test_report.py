"""Tests for the evaluation-report assembly (no PK-Sim, no API key).

Locks the honesty-critical behaviour: parameters are labelled by what actually
happened in the run (estimated vs fixed), fixed parameters appear as their own
rows, and no narrative section is ever shipped empty.
"""

import os
import unittest

from pkpd_agent.report import assemble, llm_narrative, deterministic_narrative, ReportData
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

    def test_flags_permeability_weak_identifiability(self):
        from pkpd_agent.report import deterministic_narrative
        d = self._data([{"name": "Permeability", "value": 1.4e-4, "unit": "cm/min",
                         "role": "estimated", "plausible_range": [1e-6, 1.0]}])
        pr = deterministic_narrative(d)["parameter_rationale"]
        self.assertIn("WEAKLY IDENTIFIABLE", pr)


if __name__ == "__main__":
    unittest.main()
