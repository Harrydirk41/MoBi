"""Biologic (protein / mAb) structure analysis, biodistribution scoring, benchmark.

Structure/benchmark tests run against the real BAY794620, dAb2 and Tefibazumab
library snapshots; the metric and CSV-column tests are synthetic (no PK-Sim).
"""

import json
import math
import os
import unittest

from pkpd_agent.engines import osp_biologic as B
from pkpd_agent.engines.osp_cli import OSPCli

_LIB = os.path.join(os.path.dirname(__file__), "..", "..", "OSP-PBPK-Model-Library")
BAY = os.path.join(_LIB, "BAY794620", "json", "BAY794620.json")
DAB = os.path.join(_LIB, "dAb2", "json", "dAb2.json")
TEFI = os.path.join(_LIB, "Tefibazumab", "json", "Tefibazumab.json")
ALF = os.path.join(_LIB, "Alfentanil", "json", "Alfentanil-Model.json")


@unittest.skipUnless(os.path.exists(BAY), "BAY794620 biologic snapshot not present")
class TestAnalyzeBiologic(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open(BAY, encoding="utf-8") as fh:
            cls.snap = json.load(fh)
        cls.b = B.analyze_biologic(cls.snap)

    def test_detects_large_molecule(self):
        self.assertTrue(self.b["is_biologic"])
        self.assertEqual(self.b["molecule"], "BAY794620")

    def test_recovers_fcrn_kd_not_enzymes(self):
        names = {p["name"] for p in self.b["disposition_parameters"]}
        self.assertIn("Kd (FcRn) in endosomal space", names)
        self.assertFalse(self.b["has_metabolizing_processes"])   # no enzymes

    def test_biodistribution_matrices(self):
        # whole blood + many tissues (this is a biodistribution, not one plasma curve)
        organs = {m["organ"] for m in self.b["observed_matrices"]}
        self.assertIn("Blood", organs)
        self.assertGreater(len(self.b["observed_matrices"]), 5)

    def test_observed_keeps_tissue(self):
        obs = B.biologic_observed(self.snap, "BAY794620")
        matrices = {o["matrix"] for o in obs}
        # tissue matrices are kept (the small-molecule extractor would drop these)
        self.assertTrue(any("Tissue" in m for m in matrices))
        self.assertTrue(any("WholeBlood" in m for m in matrices))

    def test_matrix_specs_carry_organ_compartment(self):
        obs = B.biologic_observed(self.snap, "BAY794620")
        specs = B.matrix_specs(obs, "BAY794620")
        blood = next(s for s in specs
                     if s["organ"] == "Blood" and s["compartment"] == "WholeBlood")
        self.assertEqual(blood["molecule"], "BAY794620")
        self.assertEqual(blood["organ"], "Blood")
        self.assertEqual(blood["compartment"], "WholeBlood")

    def test_two_studies_are_separate_matrices(self):
        # BAY794620 has Autoradiography (500 ng) and TissueDissection (1.25 mg/kg)
        # studies; the same organ in each must be a SEPARATE matrix, so it is
        # scored against the simulation run at its own dose.
        obs = B.biologic_observed(self.snap, "BAY794620")
        blood_studies = {o["study"] for o in obs if o["organ"] == "Blood"}
        self.assertIn("autoradiography", blood_studies)
        self.assertIn("tissuedissection", blood_studies)
        # and the matrix keys differ by study
        keys = {o["matrix"] for o in obs if o["organ"] == "Blood"}
        self.assertGreaterEqual(len(keys), 2)

    def test_single_molecule_is_not_biologic(self):
        if os.path.exists(ALF):
            with open(ALF, encoding="utf-8") as fh:
                self.assertEqual(B.analyze_biologic(json.load(fh)), {})


@unittest.skipUnless(os.path.exists(DAB), "dAb2 snapshot not present")
class TestGfrBiologic(unittest.TestCase):
    def test_recovers_gfr_fraction(self):
        with open(DAB, encoding="utf-8") as fh:
            b = B.analyze_biologic(json.load(fh))
        names = {p["name"] for p in b["disposition_parameters"]}
        self.assertIn("GFR fraction", names)


class TestBiologicColumnAndScoring(unittest.TestCase):
    def test_token_column_is_space_insensitive(self):
        # observed organ 'Peripheral Venous Blood' must match the CSV's
        # 'PeripheralVenousBlood' (no spaces)
        header = ["Time [h]",
                  "Organism|PeripheralVenousBlood|Tefibazumab|Plasma [µmol/l]",
                  "Organism|Fat|Tefibazumab|Tissue [µmol/l]"]
        i = OSPCli._pick_column_by_tokens(
            header, ["Tefibazumab", "Peripheral Venous Blood", "Plasma"])
        self.assertEqual(i, 1)

    def test_token_column_picks_right_tissue(self):
        header = ["Time [h]",
                  "Organism|Fat|mAb|Tissue [µg/ml]",
                  "Organism|Liver|mAb|Tissue [µg/ml]"]
        self.assertEqual(OSPCli._pick_column_by_tokens(header, ["mAb", "Liver"]), 2)
        self.assertIsNone(OSPCli._pick_column_by_tokens(header, ["mAb", "Spleen"]))

    # -- the segment matcher, against the REAL PK-Sim mAb column format --
    _REAL = ["﻿\"IndividualId\"", "Time [min]",
             "Organism|PeripheralVenousBlood|BAY794620|Plasma (Peripheral Venous Blood) [µmol/l]",
             "Organism|VenousBlood|BAY794620|Whole Blood [µmol/l]",
             "Organism|Fat|BAY794620|Whole Organ incl. FcRn_Complex-BAY794620 [µmol/l]",
             "Organism|Gonads|BAY794620|Whole Organ incl. FcRn_Complex-BAY794620 [µmol/l]",
             "Organism|SmallIntestine|BAY794620|Whole Organ (Small Intestine) incl. FcRn_Complex-BAY794620 [µmol/l]"]

    def test_segment_match_real_tissue(self):
        self.assertEqual(OSPCli._pick_matrix_column(self._REAL, "BAY794620", "Fat"), 4)

    def test_ovaries_aliases_to_gonads(self):
        # observed 'Ovaries' must map to the model's 'Gonads' column
        self.assertEqual(OSPCli._pick_matrix_column(self._REAL, "BAY794620", "Ovaries"), 5)

    def test_intestine_does_not_match_smallintestine(self):
        # the substring 'Intestine' must NOT hijack the 'SmallIntestine' column
        self.assertIsNone(OSPCli._pick_matrix_column(self._REAL, "BAY794620", "Intestine"))
        self.assertEqual(
            OSPCli._pick_matrix_column(self._REAL, "BAY794620", "SmallIntestine"), 6)

    def test_whole_blood_vs_plasma_disambiguated_by_compartment(self):
        # 'Blood'+'WholeBlood' -> the Whole Blood column, not the plasma one
        self.assertEqual(OSPCli._pick_matrix_column(
            self._REAL, "BAY794620", "Blood", "WholeBlood"), 3)

    def test_tumor_has_no_model_column(self):
        self.assertIsNone(OSPCli._pick_matrix_column(self._REAL, "BAY794620", "Tumor"))

    def _obs(self, matrix, scale=1.0):
        t = [0, 24, 72, 168]
        return {"dataset": matrix, "matrix": matrix, "molecule": "mAb",
                "time_h": t, "conc_mg_L": [scale * math.exp(-0.01 * x) for x in t]}

    def _pred(self, scale=1.0):
        class P:
            pass
        p = P()
        p.time_h = [0, 24, 72, 168]
        p.conc_mg_L = [scale * math.exp(-0.01 * x) for x in p.time_h]
        return p

    def test_perfect_biodistribution_scores_one(self):
        obs = [self._obs("Blood|WholeBlood"), self._obs("Liver|Tissue", 0.5)]
        pbm = {"Blood|WholeBlood": [self._pred()], "Liver|Tissue": [self._pred(0.5)]}
        s = B.score_biologic(obs, pbm)
        self.assertAlmostEqual(s["overall"]["gmfe"], 1.0, places=2)
        self.assertEqual(s["n_matched"], 2)

    def test_unmatched_tissue_flagged(self):
        obs = [self._obs("Blood|WholeBlood"), self._obs("Spleen|Tissue")]
        pbm = {"Blood|WholeBlood": [self._pred()]}     # spleen column absent
        s = B.score_biologic(obs, pbm)
        self.assertFalse(s["per_matrix"]["Spleen|Tissue"]["matched"])
        self.assertEqual(s["n_matched"], 1)

    def test_scores_each_study_against_its_own_simulation(self):
        # two studies at two doses -> two simulations. The low-dose observed must
        # be scored against the LOW-dose sim, not whichever sorts first, else a
        # dose mismatch shows up as a systematic bias.
        class P:
            def __init__(s, sim, scale):
                s.simulation = sim
                s.time_h = [0, 24, 72, 168]
                s.conc_mg_L = [scale * math.exp(-0.01 * x) for x in s.time_h]
        lo = {"dataset": "d", "matrix": "lo|Blood|WholeBlood", "molecule": "mAb",
              "study": "lo", "time_h": [0, 24, 72, 168],
              "conc_mg_L": [1.0 * math.exp(-0.01 * x) for x in [0, 24, 72, 168]]}
        hi = {"dataset": "d", "matrix": "hi|Blood|WholeBlood", "molecule": "mAb",
              "study": "hi", "time_h": [0, 24, 72, 168],
              "conc_mg_L": [50.0 * math.exp(-0.01 * x) for x in [0, 24, 72, 168]]}
        # each matrix carries BOTH sims' profiles; the scorer must pick by study
        pbm = {"lo|Blood|WholeBlood": [P("mAb_hi", 50.0), P("mAb_lo", 1.0)],
               "hi|Blood|WholeBlood": [P("mAb_hi", 50.0), P("mAb_lo", 1.0)]}
        s = B.score_biologic([lo, hi], pbm)
        # both match their OWN dose sim -> GMFE ~1, not ~50
        self.assertAlmostEqual(s["per_matrix"]["lo|Blood|WholeBlood"]["gmfe"], 1.0, places=1)
        self.assertAlmostEqual(s["per_matrix"]["hi|Blood|WholeBlood"]["gmfe"], 1.0, places=1)


class TestBiologicOrchestration(unittest.TestCase):
    def test_run_scores_matrices_from_one_run(self):
        class P:
            def __init__(s, t, c):
                s.simulation, s.study, s.route, s.dose = "sim", None, None, None
                s.time_h, s.conc_mg_L = t, c

        class Cli:
            def build_and_run(s, path, edits=None, simulations=None,
                              prune_simulations=False, target_matrices=None):
                assert target_matrices
                t = [0, 24, 72, 168]
                return {"ok": True, "profiles_by_matrix": {
                    m["key"]: [P(t, [math.exp(-0.01 * x) for x in t])]
                    for m in target_matrices}}

        snap = {"ObservedData": [
            {"Name": "blood", "ExtendedProperties": [
                {"Name": "Molecule", "Value": "mAb"},
                {"Name": "Organ", "Value": "Blood"},
                {"Name": "Compartment", "Value": "WholeBlood"}],
             "BaseGrid": {"Values": [0, 24, 72, 168], "Unit": "h"},
             "Columns": [{"Unit": "mg/l", "DataInfo": {"MolWeight": 150000.0},
                          "Values": [math.exp(-0.01 * x) for x in [0, 24, 72, 168]]}]}]}
        bstruct = {"molecule": "mAb", "disposition_parameters": []}
        out = B.run_biologic_prediction(Cli(), "x.json", bstruct, snapshot=snap)
        self.assertTrue(out["ok"])
        self.assertAlmostEqual(out["score"]["overall"]["gmfe"], 1.0, places=1)


class TestBiologicBenchmarkGenerator(unittest.TestCase):
    @unittest.skipUnless(os.path.exists(BAY), "BAY794620 snapshot not present")
    def test_generates_leak_free_biologic_benchmark(self):
        import examples.build_biologic_benchmark as G
        res = G.build(BAY)
        self.assertFalse(res.get("skip"))
        fitted = res["answer_edits"]["disposition_parameters"][0]["parameters"]
        kd = next(p for p in fitted if p["name"] == "Kd (FcRn) in endosomal space")
        # the fitted value is blanked in the snapshot
        for c in res["blanked"]["Compounds"]:
            for p in c.get("Parameters") or []:
                if p.get("Name") == "Kd (FcRn) in endosomal space":
                    self.assertNotAlmostEqual(p["Value"], kd["value"])
                    self.assertEqual((p.get("ValueOrigin") or {}).get("Source"), "Unknown")
        # and never leaks into the agent input
        self.assertNotIn(str(kd["value"]), json.dumps(res["input"]))

    @unittest.skipUnless(os.path.exists(BAY), "BAY794620 snapshot not present")
    def test_input_forbids_enzyme_clearance(self):
        import examples.build_biologic_benchmark as G
        res = G.build(BAY)
        guidance = res["input"]["unknowns_guidance"].lower()
        self.assertIn("not a small molecule", guidance)
        self.assertIn("do not add enzyme", guidance)


if __name__ == "__main__":
    unittest.main()
