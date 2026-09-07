"""HARD (mechanism-discovery) benchmark invariants + the two subtle no-leak fixes.

The HARD variant must (a) strip the model STRUCTURE and (b) WITHHOLD the metabolizing
enzyme identity, so the agent builds the model and discovers the mechanism from a
candidate pool. Two regressions are guarded here specifically:

  * ``_blank`` (normal mode too) must neutralize the partition/permeability methods on
    BOTH the compound and its simulation copies - else the reference's method leaks.
  * the candidate pool must be DEDUPED by molecule (the expression system lists a
    molecule once per tissue), else a single-enzyme model looks like a big pool.
"""

import copy
import os
import sys
import unittest

_HERE = os.path.dirname(os.path.abspath(__file__))
_LIB = os.path.abspath(os.path.join(_HERE, "..", "..", "OSP-PBPK-Model-Library"))
if _LIB not in sys.path:
    sys.path.insert(0, _LIB)

import build_benchmark as G                                        # noqa: E402

_ALFENTANIL = os.path.join(_LIB, "Alfentanil", "json", "Alfentanil-Model.json")


@unittest.skipUnless(os.path.exists(_ALFENTANIL), "Alfentanil snapshot not present")
class TestHardDiscovery(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.res = G.build_files(_ALFENTANIL)
        cls.ans = sorted({p["molecule"] for p in
                          cls.res["answer_key"]["structural_choices"]["metabolizing_processes"]
                          if p.get("molecule")})

    def test_structure_stripped(self):
        hb = self.res["hard_blanked"]
        for c in hb.get("Compounds") or []:
            self.assertEqual(c.get("Processes"), [], "compound processes must be stripped")
        for s in hb.get("Simulations") or []:
            for sc in s.get("Compounds") or []:
                self.assertFalse(sc.get("Processes"), "sim process refs must be cleared")

    def test_enzyme_identity_withheld(self):
        facts = " ".join(self.res["hard_input"]["background"]["literature_facts"])
        for mol in self.ans:
            self.assertNotIn(mol, facts, f"HARD facts must not NAME the answer enzyme {mol}")
        self.assertIn("metabolism", facts.lower(),
                      "HARD facts must still state that metabolism is present")

    def test_candidate_pool_honest(self):
        pool = [m["molecule"] for m in
                self.res["hard_input"]["background"]["candidate_clearance_molecules"]]
        self.assertEqual(len(pool), len(set(pool)), "candidate pool must be deduped by molecule")
        for mol in self.ans:
            self.assertIn(mol, pool, "pool must contain the answer (task must be solvable)")
        self.assertGreater(len(set(pool)) - len(self.ans), 0,
                           "Alfentanil is a genuine discovery test: pool must have distractors")

    def test_normal_blank_neutralizes_methods_everywhere(self):
        """The reference partition/permeability method must not survive in the NORMAL
        blanked snapshot - on the compound OR its simulation copies."""
        blanked = self.res["blanked"]
        comps = list(blanked.get("Compounds") or [])
        for s in blanked.get("Simulations") or []:
            comps += list(s.get("Compounds") or [])
        for c in comps:
            for m in c.get("CalculationMethods") or []:
                low = str(m).lower()
                if "partition" in low or "permeability" in low:
                    self.assertTrue(str(m).endswith("PK-Sim Standard"),
                                    f"non-neutral method leaked: {m}")


class TestCandidatePoolDedup(unittest.TestCase):
    def test_pool_dedups_per_tissue_entries(self):
        """A molecule expressed in many tissues (many ExpressionProfiles entries) must
        appear ONCE in the candidate pool."""
        base = {
            "objective": "x", "compound": "D",
            "background": {"description": "d", "literature_facts": []},
            "given_data": {},
        }
        comp = {"Name": "D", "Processes": [
            {"Molecule": "CYP3A4", "InternalName": "MetabolizationIntrinsic_FirstOrder"}]}
        eps = [{"Molecule": "CYP3A4", "Type": "Enzyme"} for _ in range(22)]  # 22 tissues
        eps += [{"Molecule": "CYP1A2", "Type": "Enzyme"} for _ in range(5)]  # a distractor
        hi = G._hard_input(copy.deepcopy(base), comp, eps)
        pool = [m["molecule"] for m in hi["background"]["candidate_clearance_molecules"]]
        self.assertEqual(sorted(pool), ["CYP1A2", "CYP3A4"])
        self.assertNotIn("CYP3A4", " ".join(hi["background"]["literature_facts"]))


if __name__ == "__main__":
    unittest.main()
