"""DDI structure analysis + interaction-ratio metric.

Structure test runs against the real Erythromycin (Midazolam victim) library
snapshot; the metric tests are synthetic (no PK-Sim).
"""

import json
import math
import os
import unittest

from pkpd_agent.engines.osp_ddi import (analyze_ddi, interaction_ratios,
                                        score_ddi, _auc)

ERY = os.path.join(os.path.dirname(__file__), "..", "..",
                   "OSP-PBPK-Model-Library", "Erythromycin", "json",
                   "Erythromycin-Model.json")
ALF = os.path.join(os.path.dirname(__file__), "..", "..",
                   "OSP-PBPK-Model-Library", "Alfentanil", "json",
                   "Alfentanil-Model.json")


@unittest.skipUnless(os.path.exists(ERY), "Erythromycin DDI snapshot not present")
class TestAnalyzeDDI(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open(ERY, encoding="utf-8") as fh:
            cls.r = analyze_ddi(json.load(fh))

    def test_identifies_perpetrator_and_victim(self):
        self.assertTrue(self.r["is_ddi"])
        self.assertEqual([p["name"] for p in self.r["perpetrators"]], ["Erythromycin"])
        self.assertIn("Midazolam", self.r["victims"])

    def test_extracts_mbi_mechanism(self):
        mechs = self.r["perpetrators"][0]["mechanisms"]
        mbi = next(m for m in mechs if m["internal_name"] == "IrreversibleInhibition")
        self.assertEqual(mbi["target"], "CYP3A4")
        self.assertIn("kinact", mbi["parameters"])
        self.assertIn("K_kinact_half", mbi["parameters"])   # ground truth, not Ki

    def test_interaction_link(self):
        it = self.r["interactions"][0]
        self.assertEqual(it["molecule"], "CYP3A4")
        self.assertEqual(it["perpetrator"], "Erythromycin")

    def test_pairs_only_true_victim_arms(self):
        # exactly the two Olkkola1993 control/treatment pairs; erythromycin's own
        # PK arms are left unpaired, not falsely matched to the midazolam control
        treats = {p["treatment"] for p in self.r["pairs"]}
        self.assertEqual(treats, {"Treatment_Olkkola1993_IV",
                                  "Treatment_Olkkola1993_Oral"})
        self.assertGreater(len(self.r["unpaired_treatments"]), 10)

    def test_single_compound_is_not_ddi(self):
        if os.path.exists(ALF):
            with open(ALF, encoding="utf-8") as fh:
                self.assertEqual(analyze_ddi(json.load(fh)), {})


class TestDDIMetric(unittest.TestCase):
    def _prof(self, decline):
        # a simple mono-exponential-ish victim profile; smaller decline -> more AUC
        t = [0, 1, 2, 4, 8]
        c = [0.1 * math.exp(-decline * x) for x in t]
        return {"time_h": t, "conc": c}

    def test_ratio_greater_than_one_when_treatment_auc_higher(self):
        prof = {"C": self._prof(0.8), "T": self._prof(0.2)}   # T declines slower
        pairs = [{"control": "C", "treatment": "T", "route": "IV"}]
        r = interaction_ratios(pairs, prof)[0]
        self.assertGreater(r["auc_ratio"], 1.0)

    def test_score_fold_error_and_gmfe(self):
        pred = [{"treatment": "T", "route": "IV", "auc_ratio": 4.0, "cmax_ratio": 1.5}]
        obs = [{"treatment": "T", "auc_ratio": 4.0, "cmax_ratio": 1.5}]
        s = score_ddi(pred, obs)
        self.assertEqual(s["gmfe_aucr"], 1.0)            # perfect ratio match
        self.assertEqual(s["within_2fold_pct"], 100.0)

    def test_score_penalises_wrong_ratio(self):
        pred = [{"treatment": "T", "auc_ratio": 1.2}]    # missed the interaction
        obs = [{"treatment": "T", "auc_ratio": 4.4}]
        s = score_ddi(pred, obs)
        self.assertGreater(s["gmfe_aucr"], 3.0)          # ~3.7x off
        self.assertEqual(s["within_2fold_pct"], 0.0)

    def test_auc_trapezoid(self):
        self.assertAlmostEqual(_auc([0, 2], [1.0, 1.0]), 2.0, places=6)


if __name__ == "__main__":
    unittest.main()
