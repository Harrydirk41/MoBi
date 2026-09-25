"""Pre-flight validation: catch model configs PK-Sim can't run (the opaque
'no result CSVs parsed' crash) and turn them into actionable messages."""

import unittest
from pkpd_agent.tools.osp_loop_tools import _structure_problems


class TestStructureProblems(unittest.TestCase):
    def test_duplicate_molecule_processes_flagged(self):
        struct = {"add_processes": [
            {"type": "metabolization_mm", "molecule": "CYP3A4"},
            {"type": "metabolization_mm", "molecule": "CYP3A4"}]}
        probs = _structure_problems(struct, {"kcat@CYP3A4": [0.1, 5]})
        self.assertEqual(len(probs), 1)
        self.assertIn("same molecule 'CYP3A4'", probs[0])
        self.assertIn("lumped", probs[0])

    def test_distinct_molecules_ok(self):
        struct = {"add_processes": [
            {"type": "metabolization_mm", "molecule": "CYP3A4"},
            {"type": "metabolization_mm", "molecule": "UGT1A4"}]}
        self.assertEqual(_structure_problems(struct, {}), [])

    def test_degenerate_link_scale_flagged(self):
        probs = _structure_problems({}, {}, link_groups=[["kcat@CYP3A4", "kcat@CYP3A4"]])
        self.assertTrue(any("fewer than two distinct" in p for p in probs))

    def test_valid_link_scale_ok(self):
        self.assertEqual(
            _structure_problems({}, {}, link_groups=[["CLspec@CYP3A4", "CLspec@UGT1A4"]]), [])

    def test_permeability_under_computed_method_flagged(self):
        struct = {"calculation_methods": {"permeability": "Charge dependent Schmitt"}}
        probs = _structure_problems(struct, {"Permeability": [1e-4, 1]})
        self.assertTrue(any("Permeability" in p and "PK-Sim Standard" in p for p in probs))

    def test_permeability_under_pksim_standard_ok(self):
        struct = {"calculation_methods": {"permeability": "PK-Sim Standard"}}
        self.assertEqual(_structure_problems(struct, {"Permeability": [1e-4, 1]}), [])


if __name__ == "__main__":
    unittest.main()


class TestResolveSims(unittest.TestCase):
    """fit_simulations from the agent (loose labels) resolve to real simulation
    names, so staged IV-first fits don't fail with 'Simulation not found'."""
    def setUp(self):
        from pkpd_agent.engines.osp_optimize import _resolve_sims
        self._r = _resolve_sims
        self.sims = ["Kroboth 1988 - IV - 0.5 mg", "Kroboth 1988 - IV - 1.0 mg",
                     "Kroboth 1988 - IV - 2.0 mg", "Smith 1984 - IV - 1 mg"]

    def test_loose_label_resolves(self):
        got, un = self._r(["Kroboth 1988 - 1.0 mg"], self.sims)
        self.assertEqual(got, ["Kroboth 1988 - IV - 1.0 mg"]); self.assertEqual(un, [])

    def test_study_only_matches_all_its_sims(self):
        got, _ = self._r(["Kroboth 1988"], self.sims)
        self.assertEqual(len(got), 3)

    def test_unmatched_reported(self):
        got, un = self._r(["Nope 2050"], self.sims)
        self.assertEqual(got, []); self.assertEqual(un, ["Nope 2050"])
