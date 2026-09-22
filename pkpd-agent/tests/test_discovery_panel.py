"""Hard-mode (mechanism-discovery) anti-leak invariants.

The clearing enzyme/transporter identity is withheld in hard mode; the candidate
pool must be a FIXED, case-independent panel with many decoys, never the
reference model's own molecules (which would hand over the answer). These tests
lock that in so the leak cannot silently regress.
"""
import glob
import json
import os
import unittest

from pkpd_agent.engines import osp_catalog

_LIB = os.path.join(os.path.dirname(__file__), "..", "..",
                    "OSP-PBPK-Model-Library")
MIN_DISTRACTORS = 8


class TestStandardPanel(unittest.TestCase):
    def test_panel_is_broad_and_fixed(self):
        panel = osp_catalog.STANDARD_DISCOVERY_PANEL
        mols = [m["molecule"] for m in panel]
        self.assertGreaterEqual(len(panel), 20)
        self.assertEqual(len(mols), len(set(mols)), "panel has duplicates")
        # covers all three mechanism classes
        types = {m["type"] for m in panel}
        self.assertTrue({"Enzyme", "Transporter", "OtherProtein"} <= types)

    def test_every_panel_molecule_has_a_real_profile(self):
        # each candidate must be materializable, or a decoy would be unattachable
        # (and unattachability itself leaks which molecule is real).
        for m in osp_catalog.STANDARD_DISCOVERY_PANEL:
            self.assertIsNotNone(
                osp_catalog.panel_expression_profile(m["molecule"]),
                f"no harvested expression profile for {m['molecule']}")

    def test_pool_unions_needed_without_shrinking(self):
        pool = osp_catalog.hard_candidate_pool(
            [{"molecule": "Cleavage-protein", "type": "Enzyme"}])
        mols = [m["molecule"] for m in pool]
        self.assertGreaterEqual(len(pool), len(osp_catalog.STANDARD_DISCOVERY_PANEL))
        self.assertIn("Cleavage-protein", mols)   # unusual answer stays reachable
        self.assertIn("CYP3A4", mols)              # panel decoys still present


class TestLibraryHasNoDiscoveryLeak(unittest.TestCase):
    """Every built hard task in the library must carry a broad pool."""

    def _answer_molecules(self, hard_input_path):
        ak = hard_input_path.replace("json_input", "answer_key").replace(
            ".hard.input.json", ".answer_key.json")
        used = set()
        if os.path.exists(ak):
            data = json.load(open(ak, encoding="utf-8"))
            for p in data.get("estimated_parameters", []):
                if "@" in p["parameter"]:
                    used.add(p["parameter"].split("@")[1])
        return used

    def test_all_hard_pools_have_enough_decoys(self):
        files = sorted(glob.glob(os.path.join(_LIB, "*", "json_input",
                                              "*.hard.input.json")))
        if not files:
            self.skipTest("OSP library not present")
        offenders = []
        for f in files:
            bg = (json.load(open(f, encoding="utf-8")).get("background") or {})
            pool = {m.get("molecule")
                    for m in bg.get("candidate_clearance_molecules") or []}
            ans = self._answer_molecules(f)
            decoys = pool - ans
            if pool and len(decoys) < MIN_DISTRACTORS:
                offenders.append((os.path.basename(f), len(decoys)))
        self.assertEqual(offenders, [], f"discovery leak (too few decoys): {offenders}")


if __name__ == "__main__":
    unittest.main()
