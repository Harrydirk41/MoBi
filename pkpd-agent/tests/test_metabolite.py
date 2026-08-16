"""Metabolite-cascade structure analysis, per-molecule scoring, and benchmark.

Structure/benchmark tests run against the real Itraconazole cascade snapshot
(Itraconazole -> Hydroxy- -> Keto- -> N-desalkyl-Itraconazole); the metric and
CLI-column tests are synthetic (no PK-Sim).
"""

import json
import math
import os
import unittest

from pkpd_agent.engines import osp_metabolite as M
from pkpd_agent.engines.osp_cli import OSPCli

ITRA = os.path.join(os.path.dirname(__file__), "..", "..",
                    "OSP-PBPK-Model-Library", "Itraconazole", "json",
                    "Itraconazole-Model.json")
ALF = os.path.join(os.path.dirname(__file__), "..", "..",
                   "OSP-PBPK-Model-Library", "Alfentanil", "json",
                   "Alfentanil-Model.json")


@unittest.skipUnless(os.path.exists(ITRA), "Itraconazole cascade snapshot not present")
class TestAnalyzeMetabolites(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open(ITRA, encoding="utf-8") as fh:
            cls.snap = json.load(fh)
        cls.m = M.analyze_metabolites(cls.snap)

    def test_detects_cascade_and_root(self):
        self.assertTrue(self.m["is_metabolite"])
        self.assertEqual(self.m["root"], "Itraconazole")

    def test_chain_order(self):
        self.assertEqual(self.m["chain"],
                         ["Itraconazole", "Hydroxy-Itraconazole",
                          "Keto-Itraconazole", "N-desalkyl-Itraconazole"])

    def test_edges_name_the_product(self):
        first = self.m["edges"][0]
        self.assertEqual(first["parent"], "Itraconazole")
        self.assertEqual(first["metabolite"], "Hydroxy-Itraconazole")
        self.assertEqual(first["enzyme"], "CYP3A4")

    def test_metabolites_have_scorable_plasma(self):
        # the parent AND at least the first metabolite carry measured plasma
        self.assertIn("Itraconazole", self.m["scorable_molecules"])
        self.assertIn("Hydroxy-Itraconazole", self.m["scorable_molecules"])

    def test_plasma_observed_is_per_molecule(self):
        parent = M.plasma_observed(self.snap, "Itraconazole")
        daughter = M.plasma_observed(self.snap, "Hydroxy-Itraconazole")
        self.assertTrue(parent)
        self.assertTrue(daughter)
        # each dataset belongs to its own molecule (no cross-contamination)
        self.assertTrue(all(o["molecule"] == "Itraconazole" for o in parent))
        self.assertTrue(all(o["molecule"] == "Hydroxy-Itraconazole" for o in daughter))

    def test_single_compound_is_not_cascade(self):
        if os.path.exists(ALF):
            with open(ALF, encoding="utf-8") as fh:
                self.assertEqual(M.analyze_metabolites(json.load(fh)), {})


_OME = os.path.join(os.path.dirname(__file__), "..", "..",
                    "OSP-PBPK-Model-Library", "Omeprazole", "json",
                    "Omeprazole-Model.json")
_CIM = os.path.join(os.path.dirname(__file__), "..", "..",
                    "OSP-PBPK-Model-Library", "Cimetidine", "json",
                    "Cimetidine-Model.json")


class TestMulticompound(unittest.TestCase):
    @unittest.skipUnless(os.path.exists(_OME), "Omeprazole snapshot not present")
    def test_parallel_enantiomers_detected(self):
        # Omeprazole = R/S enantiomers + racemate, each with plasma + its own
        # clearance, NO metabolite cascade edge -> a PARALLEL multicompound task.
        with open(_OME, encoding="utf-8") as fh:
            snap = json.load(fh)
        mc = M.analyze_multicompound(snap)
        self.assertEqual(mc["kind"], "parallel")
        self.assertFalse(mc["edges"])            # no cascade
        self.assertGreaterEqual(len(mc["scorable_molecules"]), 2)
        # the cascade-only analyzer must NOT claim it
        self.assertEqual(M.analyze_metabolites(snap), {})

    @unittest.skipUnless(os.path.exists(_CIM), "Cimetidine snapshot not present")
    def test_victim_ddi_excluded_from_multicompound(self):
        # a perpetrator->victim DDI is the DDI task, not a joint parallel fit
        with open(_CIM, encoding="utf-8") as fh:
            self.assertEqual(M.analyze_multicompound(json.load(fh)), {})

    @unittest.skipUnless(os.path.exists(ITRA), "Itraconazole snapshot not present")
    def test_cascade_still_detected_by_multicompound(self):
        with open(ITRA, encoding="utf-8") as fh:
            mc = M.analyze_multicompound(json.load(fh))
        self.assertEqual(mc["kind"], "cascade")
        self.assertTrue(mc["edges"])

    @unittest.skipUnless(os.path.exists(_OME), "Omeprazole snapshot not present")
    def test_parallel_benchmark_has_no_hard_mode(self):
        import examples.build_metabolite_benchmark as G
        res = G.build(_OME)
        self.assertFalse(res.get("skip"))
        self.assertEqual(res["kind"], "parallel")
        self.assertFalse(res["is_cascade"])
        self.assertIsNone(res["hard_input"])          # no cascade to strip
        self.assertIn("modelled_compounds", res["input"])


class TestMetaboliteScoring(unittest.TestCase):
    def _prof(self, name, decline, scale=1.0):
        class P:
            pass
        p = P()
        p.simulation, p.study, p.route, p.dose = name, "S2000", "PO", "100 mg"
        p.time_h = [0, 1, 2, 4, 8]
        p.conc_mg_L = [scale * 0.1 * math.exp(-decline * t) for t in p.time_h]
        return p

    def _obs(self, decline, scale=1.0):
        return [{"dataset": "S2000, PO", "study": "S2000", "route": "PO",
                 "dose": "100 mg", "molecule": "X",
                 "time_h": [0, 1, 2, 4, 8],
                 "conc_mg_L": [scale * 0.1 * math.exp(-decline * t)
                               for t in [0, 1, 2, 4, 8]]}]

    def test_perfect_cascade_scores_near_one(self):
        obs = {"Parent": self._obs(0.3), "Meta": self._obs(0.15)}
        prof = {"Parent": [self._prof("S2000, PO", 0.3)],
                "Meta": [self._prof("S2000, PO", 0.15)]}
        s = M.score_metabolites(obs, prof)
        self.assertAlmostEqual(s["per_molecule"]["Parent"]["overall"]["gmfe"], 1.0, places=2)
        self.assertAlmostEqual(s["per_molecule"]["Meta"]["overall"]["gmfe"], 1.0, places=2)
        self.assertAlmostEqual(s["cascade"]["gmfe"], 1.0, places=2)

    def test_metabolite_miss_penalises_cascade(self):
        # parent fits perfectly but the metabolite is 3x under-predicted: a
        # parent-only score would look great; the cascade score must not.
        obs = {"Parent": self._obs(0.3), "Meta": self._obs(0.15, scale=3.0)}
        prof = {"Parent": [self._prof("S2000, PO", 0.3)],
                "Meta": [self._prof("S2000, PO", 0.15, scale=1.0)]}
        s = M.score_metabolites(obs, prof)
        self.assertAlmostEqual(s["per_molecule"]["Parent"]["overall"]["gmfe"], 1.0, places=2)
        self.assertGreater(s["per_molecule"]["Meta"]["overall"]["gmfe"], 2.5)
        self.assertGreater(s["cascade"]["gmfe"], 1.3)   # cascade caught the miss

    def test_identifiability_flags_dataless_metabolite(self):
        mstruct = {"unscorable_metabolites": ["Keto-X"]}
        acts = M.metabolite_identifiability(mstruct)
        self.assertEqual(acts[0]["molecule"], "Keto-X")
        self.assertEqual(acts[0]["severity"], "medium")


class TestMultiMoleculeColumn(unittest.TestCase):
    def test_strict_match_picks_own_column_not_parent(self):
        # a cascade CSV has a plasma column per molecule; strict matching must
        # return the daughter's OWN column and NOT fall back to the parent's.
        header = ["Time [h]",
                  "Organism|PeripheralVenousBlood|Itraconazole|Plasma [µmol/l]",
                  "Organism|PeripheralVenousBlood|Hydroxy-Itraconazole|Plasma [µmol/l]"]
        i_parent = OSPCli._pick_conc_column(header, "Itraconazole", True)
        i_daughter = OSPCli._pick_conc_column(header, "Hydroxy-Itraconazole", True)
        self.assertEqual(i_parent, 1)
        self.assertEqual(i_daughter, 2)

    def test_strict_match_returns_none_when_absent(self):
        header = ["Time [h]",
                  "Organism|PeripheralVenousBlood|Itraconazole|Plasma [µmol/l]"]
        # a molecule with no column must be None under strict matching, NOT the
        # parent's column (which would score it against the wrong compound)
        self.assertIsNone(OSPCli._pick_conc_column(header, "Keto-Itraconazole", True))

    def test_nonstrict_falls_back_to_plasma(self):
        header = ["Time [h]",
                  "Organism|PeripheralVenousBlood|Itraconazole|Plasma [µmol/l]"]
        self.assertEqual(OSPCli._pick_conc_column(header, "Keto-Itraconazole", False), 1)


class TestMetaboliteOrchestration(unittest.TestCase):
    def test_run_scores_every_molecule_from_one_run(self):
        class P:
            def __init__(s, n, t, c):
                s.simulation, s.study, s.route, s.dose = n, "S2000", "PO", "100 mg"
                s.time_h, s.conc_mg_L = t, c

        class Cli:
            def build_and_run(s, path, edits=None, simulations=None,
                              prune_simulations=False, target_molecules=None):
                assert target_molecules == ["Parent", "Meta"]
                t = [0, 1, 2, 4, 8]
                return {"ok": True, "profiles_by_molecule": {
                    "Parent": [P("S2000, PO", t, [0.1 * math.exp(-0.3 * x) for x in t])],
                    "Meta": [P("S2000, PO", t, [0.1 * math.exp(-0.15 * x) for x in t])]}}

        snap = {"ObservedData": [
            {"Name": "S2000, PO",
             "ExtendedProperties": [{"Name": "Molecule", "Value": mol},
                                    {"Name": "Compartment", "Value": "Plasma"},
                                    {"Name": "Study Id", "Value": "S2000"},
                                    {"Name": "Route", "Value": "PO"}],
             "BaseGrid": {"Values": [0, 1, 2, 4, 8], "Unit": "h"},
             "Columns": [{"Unit": "mg/l", "DataInfo": {"MolWeight": 700.0},
                          "Values": [0.1 * math.exp(-d * x) for x in [0, 1, 2, 4, 8]]}]}
            for mol, d in (("Parent", 0.3), ("Meta", 0.15))]}
        mstruct = {"root": "Parent", "scorable_molecules": ["Parent", "Meta"]}
        out = M.run_metabolite_prediction(Cli(), "x.json", mstruct, snapshot=snap)
        self.assertTrue(out["ok"])
        self.assertIn("cascade", out["score"])
        self.assertAlmostEqual(out["score"]["cascade"]["gmfe"], 1.0, places=1)


class TestMetaboliteBenchmarkGenerator(unittest.TestCase):
    @unittest.skipUnless(os.path.exists(ITRA), "Itraconazole snapshot not present")
    def test_generates_leak_free_cascade_benchmark(self):
        import examples.build_metabolite_benchmark as G
        res = G.build(ITRA)
        self.assertFalse(res.get("skip"))
        # answer = fitted disposition params for the cascade compounds
        specs = res["answer_edits"]["cascade_parameters"]
        itra = next(s for s in specs if s["compound"] == "Itraconazole")
        names = {p["name"] for p in itra["parameters"]}
        self.assertIn("Lipophilicity", names)
        # each fitted value is blanked in the snapshot (exact value gone from its
        # own parameter slot)
        blanked = {c["Name"]: c for c in res["blanked"]["Compounds"]}
        for s in specs:
            comp = blanked[s["compound"]]
            fitted_now = []

            def w(o):
                if isinstance(o, dict):
                    if isinstance(o.get("Name"), str) and \
                            (o.get("ValueOrigin") or {}).get("Source") == "ParameterIdentification":
                        fitted_now.append(o["Name"])
                    for v in o.values():
                        w(v)
                elif isinstance(o, list):
                    for v in o:
                        w(v)
            w(comp)
            # nothing in the blanked snapshot still carries a ParameterIdentification
            # origin (all reset to Unknown / naive prior)
            self.assertEqual(fitted_now, [], f"{s['compound']} leaked a fitted origin")

    @unittest.skipUnless(os.path.exists(ITRA), "Itraconazole snapshot not present")
    def test_hard_mode_removes_cascade_structure(self):
        import examples.build_metabolite_benchmark as G
        res = G.build(ITRA)
        # normal keeps the Metabolite edges; hard strips them
        norm = any(p.get("Metabolite") for c in res["blanked"]["Compounds"]
                   for p in (c.get("Processes") or []))
        hard = any(p.get("Metabolite") for c in res["hard_blanked"]["Compounds"]
                   for p in (c.get("Processes") or []))
        self.assertTrue(norm)
        self.assertFalse(hard)
        self.assertEqual(len(res["removed_processes"]), 3)
        # the cascade topology (the answer) must NOT be echoed in the hard input
        self.assertNotIn("cascade_structure", res["hard_input"])
        # but the biology pathway facts remain (the only structure source)
        facts = res["hard_input"]["background"]["metabolic_pathway"]
        self.assertTrue(any("metabolized by CYP3A4 to Hydroxy-Itraconazole" in f
                            for f in facts))


if __name__ == "__main__":
    unittest.main()
