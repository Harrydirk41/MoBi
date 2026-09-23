"""Reference library (library-assisted mode) invariants.

The library is the OTHER compounds' finished models, offered to the agent as a
human modeler's model library. The inviolable rule is LEAVE-ONE-OUT: the target
compound's own reference must never appear.
"""
import os
import unittest

from pkpd_agent.engines import reference_library as rl

_LIB = os.path.join(os.path.dirname(__file__), "..", "..",
                    "OSP-PBPK-Model-Library")
_TRIA = os.path.join(_LIB, "Triazolam", "benchmark", "Triazolam-Model.blanked.json")


@unittest.skipUnless(os.path.isfile(_TRIA), "OSP library not present")
class TestReferenceLibrary(unittest.TestCase):
    def test_leave_one_out_excludes_target(self):
        for st in (False, True):
            lib = rl.library_for_snapshot(_TRIA, same_type=st, full=False)
            names = [m["compound"] for m in lib["models"]]
            self.assertNotIn("Triazolam", names, "target leaked into its own library")
            self.assertEqual(len(names), len(set(names)), "compound variants not deduped")
            self.assertGreater(len(names), 0)

    def test_same_type_shares_target_enzyme(self):
        # Triazolam is a CYP3A4 substrate; same-type must all involve CYP3A4.
        self.assertEqual(rl.target_clearing_molecules(_TRIA), ["CYP3A4"])
        lib = rl.library_for_snapshot(_TRIA, same_type=True, full=False)
        for m in lib["models"]:
            mols = {p["molecule"] for p in m["processes"]}
            mols |= {p.split("@", 1)[1] for p in m["estimated_parameters"]
                     if isinstance(p, str) and "@" in p}
            self.assertIn("CYP3A4", mols, f"{m['compound']} is not a CYP3A4 analogue")
        # same-type is a subset of all
        allm = rl.library_for_snapshot(_TRIA, same_type=False, full=False)
        self.assertLessEqual(len(lib["models"]), len(allm["models"]))

    def test_structural_hides_values_full_shows_them(self):
        struct = rl.library_for_snapshot(_TRIA, same_type=True, full=False)
        full = rl.library_for_snapshot(_TRIA, same_type=True, full=True)
        # structural: estimated_parameters are plain names (strings)
        ep_s = struct["models"][0]["estimated_parameters"]
        self.assertTrue(all(isinstance(x, str) for x in ep_s))
        # full: estimated_parameters carry values (dicts)
        ep_f = next(m["estimated_parameters"] for m in full["models"]
                    if m["estimated_parameters"])
        self.assertTrue(all(isinstance(x, dict) and "value" in x for x in ep_f))


if __name__ == "__main__":
    unittest.main()
