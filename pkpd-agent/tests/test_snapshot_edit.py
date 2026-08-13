"""Tests for the snapshot edit engine (the agent's PBPK action space)."""

import unittest

from pkpd_agent.engines.snapshot_edit import apply_edits

SNAP = {
    "Compounds": [{
        "Name": "Drug",
        "Lipophilicity": [{"Name": "Fit", "Parameters": [
            {"Name": "Lipophilicity", "Value": 2.0,
             "ValueOrigin": {"Source": "ParameterIdentification"}}]}],
        "CalculationMethods": [
            "Cellular partition coefficient method - Rodgers and Rowland",
            "Cellular permeability - PK-Sim Standard"],
        "Processes": [
            {"InternalName": "MetabolizationIntrinsic_FirstOrder", "Molecule": "CYP3A4",
             "Parameters": [{"Name": "Intrinsic clearance", "Value": 0.5}]},
            {"InternalName": "GlomerularFiltration", "DataSource": "GFR",
             "Parameters": [{"Name": "GFR fraction", "Value": 0.06}]}],
        "Parameters": [{"Name": "Molecular weight", "Value": 400.0}],
    }],
}


class TestApplyEdits(unittest.TestCase):
    def test_parameter_set_and_origin_cleared(self):
        out, rep = apply_edits(SNAP, {"parameters": {"Lipophilicity": 3.3}})
        lipo = out["Compounds"][0]["Lipophilicity"][0]["Parameters"][0]
        self.assertEqual(lipo["Value"], 3.3)
        self.assertNotIn("ValueOrigin", lipo)          # user value, not a fit result
        self.assertEqual(rep["parameters"]["Lipophilicity"], 3.3)

    def test_input_not_mutated(self):
        apply_edits(SNAP, {"parameters": {"Lipophilicity": 9.9}})
        self.assertEqual(
            SNAP["Compounds"][0]["Lipophilicity"][0]["Parameters"][0]["Value"], 2.0)

    def test_partition_method_swap(self):
        out, rep = apply_edits(SNAP, {"calculation_methods":
                                      {"partition": "PK-Sim Standard"}})
        cms = out["Compounds"][0]["CalculationMethods"]
        self.assertIn("Cellular partition coefficient method - PK-Sim Standard", cms)
        # the permeability entry is untouched
        self.assertIn("Cellular permeability - PK-Sim Standard", cms)
        self.assertEqual(len(cms), 2)
        self.assertEqual(rep["calculation_methods"]["partition"], "PK-Sim Standard")

    def test_process_disable(self):
        out, rep = apply_edits(SNAP, {"processes": {"CYP3A4": False}})
        mols = [p.get("Molecule") or p.get("InternalName")
                for p in out["Compounds"][0]["Processes"]]
        self.assertNotIn("CYP3A4", mols)
        self.assertIn("GlomerularFiltration", mols)
        self.assertEqual(rep["processes"]["CYP3A4"], "disabled")

    def test_gfr_alias_disable(self):
        out, _ = apply_edits(SNAP, {"processes": {"GFR": False}})
        internals = [p.get("InternalName") for p in out["Compounds"][0]["Processes"]]
        self.assertNotIn("GlomerularFiltration", internals)

    def test_not_found_reported(self):
        _, rep = apply_edits(SNAP, {"parameters": {"No Such Param": 1.0}})
        self.assertIn("parameter:No Such Param", rep["not_found"])

    def test_empty_edits_noop(self):
        out, rep = apply_edits(SNAP, None)
        self.assertEqual(out["Compounds"][0]["Lipophilicity"][0]["Parameters"][0]["Value"], 2.0)
        self.assertEqual(rep["not_found"], [])


if __name__ == "__main__":
    unittest.main()
