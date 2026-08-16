"""The DDI LLM-loop tools: ddi_inspect (observe) and ddi_try_model (act).

Structure comes from the real Erythromycin DDI snapshot; the run is monkeypatched
so no PK-Sim is needed.
"""

import json
import os
import unittest

from pkpd_agent.config import AgentConfig
from pkpd_agent.engines import osp_ddi
from pkpd_agent.state import ModelingSession
from pkpd_agent.tools.registry import ToolRegistry
from pkpd_agent.tools.osp_ddi_loop_tools import register_osp_ddi_loop_tools

ERY = os.path.join(os.path.dirname(__file__), "..", "..",
                   "OSP-PBPK-Model-Library", "Erythromycin", "json",
                   "Erythromycin-Model.json")


@unittest.skipUnless(os.path.exists(ERY), "Erythromycin DDI snapshot not present")
class TestDDILoopTools(unittest.TestCase):
    def setUp(self):
        with open(ERY, encoding="utf-8") as fh:
            snap = json.load(fh)
        self.ddi = osp_ddi.analyze_ddi(snap)
        self.observed = [{"treatment": "Treatment_Olkkola1993_Oral",
                          "auc_ratio": 4.4, "observed_aucr": 4.4}]
        self.reg = ToolRegistry()
        register_osp_ddi_loop_tools(self.reg, AgentConfig(mock=False), {
            "cli": object(), "snapshot_path": ERY, "ddi": self.ddi,
            "victim": "Midazolam", "observed_ratios": self.observed,
            "input": {"objective": "Predict erythromycin -> midazolam."}})
        self.inspect = self.reg.get("ddi_inspect").handler
        self.try_model = self.reg.get("ddi_try_model").handler

    def test_inspect_exposes_mechanism_and_identifiability(self):
        r = self.inspect({}, ModelingSession(goal="g"))
        self.assertTrue(r.ok)
        self.assertEqual(r.data["victim"], "Midazolam")
        mechs = {m["internal_name"] for m in r.data["mechanisms"]}
        self.assertIn("IrreversibleInhibition", mechs)   # MBI
        # MBI is a 2-parameter trade-off -> flagged from a single observed ratio
        recs = r.data["identifiability_recommendations"]
        self.assertTrue(any(a["mechanism"] == "IrreversibleInhibition" for a in recs))

    def test_inspect_reports_current_parameters(self):
        r = self.inspect({}, ModelingSession(goal="g"))
        mbi = next(m for m in r.data["mechanisms"]
                   if m["internal_name"] == "IrreversibleInhibition")
        self.assertIn("kinact", mbi["current_parameters"])
        self.assertIn("K_kinact_half", mbi["current_parameters"])

    def test_try_model_runs_scores_and_tracks_best(self):
        # monkeypatch the engine run so no PK-Sim is needed
        orig = osp_ddi.run_ddi_prediction
        calls = {}

        def fake(cli, path, ddi, victim, edits=None, observed_ratios=None):
            calls["edits"] = edits
            return {"ok": True, "predicted_ratios": [
                        {"treatment": "Treatment_Olkkola1993_Oral", "auc_ratio": 4.0}],
                    "score": {"gmfe_aucr": 1.1, "within_2fold_pct": 100.0,
                              "per_arm": [{"treatment": "Treatment_Olkkola1993_Oral",
                                           "predicted_aucr": 4.0, "observed_aucr": 4.4,
                                           "fold_error": 0.91}]}}
        osp_ddi.run_ddi_prediction = fake
        try:
            sess = ModelingSession(goal="g")
            ip = [{"perpetrator": "Erythromycin", "internal_name": "IrreversibleInhibition",
                   "target": "CYP3A4", "parameters": {"kinact": 0.05, "K_kinact_half": 5.0}}]
            r = self.try_model({"interaction_parameters": ip}, sess)
            self.assertTrue(r.ok)
            self.assertEqual(r.data["gmfe_aucr"], 1.1)
            self.assertEqual(calls["edits"], {"interaction_parameters": ip})
            self.assertEqual(sess.get("ddi_best_gmfe"), 1.1)
            self.assertEqual(r.data["iteration"], 1)
        finally:
            osp_ddi.run_ddi_prediction = orig

    def test_try_model_surfaces_engine_error(self):
        orig = osp_ddi.run_ddi_prediction
        osp_ddi.run_ddi_prediction = lambda *a, **k: {"ok": False, "message": "boom"}
        try:
            r = self.try_model({"interaction_parameters": []}, ModelingSession(goal="g"))
            self.assertFalse(r.ok)
            self.assertIn("boom", r.message)
        finally:
            osp_ddi.run_ddi_prediction = orig


if __name__ == "__main__":
    unittest.main()
