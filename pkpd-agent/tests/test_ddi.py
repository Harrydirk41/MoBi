"""DDI structure analysis + interaction-ratio metric.

Structure test runs against the real Erythromycin (Midazolam victim) library
snapshot; the metric tests are synthetic (no PK-Sim).
"""

import json
import math
import os
import unittest

from pkpd_agent.engines import osp_ddi
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

    def test_identifiability_flags_mbi_pair_not_single_ki(self):
        # MBI (kinact/K_kinact_half) is a 2-parameter trade-off -> flagged from a
        # single ratio; the P-gp competitive Ki (1 parameter) is identifiable ->
        # never flagged.
        acts = osp_ddi.ddi_identifiability(self.r, n_observed_ratios=1)
        mechs = {a["mechanism"] for a in acts}
        self.assertIn("IrreversibleInhibition", mechs)
        self.assertNotIn("CompetitiveInhibition", mechs)
        mbi = next(a for a in acts if a["mechanism"] == "IrreversibleInhibition")
        self.assertEqual(mbi["severity"], "high")
        self.assertIn("K_kinact_half", mbi["action"])

    def test_identifiability_relaxes_with_more_arms(self):
        # given enough independent ratios to separate both parameters, do not flag
        acts = osp_ddi.ddi_identifiability(self.r, n_observed_ratios=5)
        self.assertFalse(any(a["mechanism"] == "IrreversibleInhibition"
                             for a in acts))

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



class TestDDIEditingAndOrchestration(unittest.TestCase):
    @unittest.skipUnless(os.path.exists(ERY), "Erythromycin snapshot not present")
    def test_edit_interaction_parameters_on_perpetrator(self):
        from pkpd_agent.engines.snapshot_edit import apply_edits
        with open(ERY, encoding="utf-8") as fh:
            snap = json.load(fh)
        new, rep = apply_edits(snap, {"interaction_parameters": [
            {"perpetrator": "Erythromycin", "internal_name": "IrreversibleInhibition",
             "target": "CYP3A4", "parameters": {"kinact": 0.05, "K_kinact_half": 5.0}}]})
        ery = next(c for c in new["Compounds"] if c["Name"] == "Erythromycin")
        mbi = next(p for p in ery["Processes"]
                   if p.get("InternalName") == "IrreversibleInhibition")
        vals = {x["Name"]: x["Value"] for x in mbi["Parameters"]}
        self.assertEqual(vals["kinact"], 0.05)
        self.assertEqual(vals["K_kinact_half"], 5.0)
        self.assertIn("Erythromycin/IrreversibleInhibition/kinact",
                      rep["interaction_parameters"])

    def test_unknown_perpetrator_reported(self):
        from pkpd_agent.engines.snapshot_edit import apply_edits
        _, rep = apply_edits({"Compounds": [{"Name": "A"}]},
                             {"interaction_parameters": [
                                 {"perpetrator": "Z", "internal_name": "X",
                                  "parameters": {"Ki": 1}}]})
        self.assertTrue(any("Z" in m for m in rep["not_found"]))

    def test_orchestration_runs_pairs_and_scores(self):
        import math

        class P:
            def __init__(s, n, t, c): s.simulation, s.time_h, s.conc_mg_L = n, t, c

        class Cli:
            def build_and_run(s, snap, edits=None, simulations=None,
                              prune_simulations=False, target_molecule=None):
                assert target_molecule == "Midazolam"
                pr = []
                for nm in simulations:
                    d = 0.2 if nm.startswith("Treatment") else 0.8
                    pr.append(P(nm, [0, 1, 2, 4, 8],
                               [0.1 * math.exp(-d * x) for x in [0, 1, 2, 4, 8]]))
                return {"ok": True, "profiles": pr}

        ddi = {"pairs": [{"control": "Control_S_IV", "treatment": "Treatment_S_IV",
                          "route": "IV"}]}
        out = osp_ddi.run_ddi_prediction(
            Cli(), "x", ddi, "Midazolam",
            observed_ratios=[{"treatment": "Treatment_S_IV", "auc_ratio": 3.0}])
        self.assertTrue(out["ok"])
        self.assertEqual(out["n_ran"], 2)
        self.assertGreater(out["predicted_ratios"][0]["auc_ratio"], 1.0)
        self.assertIn("gmfe_aucr", out["score"])

class TestDDIBenchmarkGenerator(unittest.TestCase):
    @unittest.skipUnless(os.path.exists(ERY), "Erythromycin snapshot not present")
    def test_generates_leak_free_ddi_benchmark(self):
        import examples.build_ddi_benchmark as G
        res = G.build(ERY)
        self.assertFalse(res.get("skip"))
        # answer = the fitted interaction params (kinact + K_kinact_half)
        specs = res["answer_edits"]["interaction_parameters"]
        mbi = next(s for s in specs if s["internal_name"] == "IrreversibleInhibition")
        self.assertEqual(set(mbi["parameters"]), {"kinact", "K_kinact_half"})
        # blanked snapshot must NOT contain the fitted values (no leak)
        ery = next(c for c in res["blanked"]["Compounds"] if c["Name"] == "Erythromycin")
        p = next(x for x in ery["Processes"]
                 if x.get("InternalName") == "IrreversibleInhibition")
        blanked = {q["Name"]: q["Value"] for q in p["Parameters"]}
        self.assertNotAlmostEqual(blanked["kinact"], mbi["parameters"]["kinact"])
        # agent input must not contain any fitted interaction value
        import json as _j
        inp = _j.dumps(res["input"])
        for s in specs:
            for v in s["parameters"].values():
                self.assertNotIn(str(v), inp)

    @unittest.skipUnless(os.path.exists(ERY), "Erythromycin snapshot not present")
    def test_observed_ratio_direction_sensible(self):
        import examples.build_ddi_benchmark as G
        res = G.build(ERY)
        ratios = res["answer_key"]["observed_interaction_ratios"]
        self.assertTrue(ratios)
        # erythromycin inhibits CYP3A4 -> midazolam exposure INCREASES (AUCR > 1)
        self.assertGreater(ratios[0]["observed_aucr"], 1.0)


if __name__ == "__main__":
    unittest.main()
