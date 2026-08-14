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


class TestAddProcess(unittest.TestCase):
    def setUp(self):
        self.snap = {
            "ExpressionProfiles": [
                {"Type": "Enzyme", "Molecule": "CYP3A4"},
                {"Type": "Enzyme", "Molecule": "AADAC"}],
            "Compounds": [{"Name": "Drug", "Processes": [
                {"InternalName": "MetabolizationIntrinsic_FirstOrder",
                 "Molecule": "CYP3A4", "Parameters": [
                     {"Name": "Intrinsic clearance", "Value": 0.5}]}]}],
        }

    def test_add_process_on_expressed_enzyme(self):
        out, rep = apply_edits(self.snap, {"add_processes": [
            {"type": "metabolization_first_order", "molecule": "AADAC",
             "parameters": {"Intrinsic clearance": 0.05}}]})
        mols = [p.get("Molecule") for p in out["Compounds"][0]["Processes"]]
        self.assertIn("AADAC", mols)
        self.assertEqual(rep["processes"]["add:AADAC"], "added")

    def test_reject_unexpressed_molecule(self):
        _, rep = apply_edits(self.snap, {"add_processes": [
            {"type": "metabolization_first_order", "molecule": "CYP2D6"}]})
        self.assertTrue(any("CYP2D6" in n for n in rep["not_found"]))

    def test_add_duplicate_is_noop(self):
        _, rep = apply_edits(self.snap, {"add_processes": [
            {"type": "metabolization_first_order", "molecule": "CYP3A4"}]})
        self.assertEqual(rep["processes"]["add:CYP3A4"], "already_present")


class TestSimulationProcessMirroring(unittest.TestCase):
    """A process lives on the compound AND on every simulation that uses it;
    editing one without the other leaves a dangling reference that makes
    PK-Sim's snapshot mapper fail. These lock in that the two stay in sync."""

    def setUp(self):
        # two simulations, each mirroring the compound's two processes
        self.snap = {
            "ExpressionProfiles": [
                {"Type": "Enzyme", "Molecule": "CYP3A4"},
                {"Type": "Transporter", "Molecule": "ABCB1"}],
            "Compounds": [{"Name": "Drug", "Processes": [
                {"InternalName": "MetabolizationIntrinsic_FirstOrder",
                 "DataSource": "1st order CL", "Molecule": "CYP3A4"},
                {"InternalName": "GlomerularFiltration", "DataSource": "GFR",
                 "Molecule": None}]}],
            "Simulations": [
                {"Name": "sim1", "Compounds": [{"Name": "Drug", "Processes": [
                    {"Name": "CYP3A4-1st order CL", "MoleculeName": "CYP3A4"},
                    {"Name": "Glomerular Filtration-GFR",
                     "SystemicProcessType": "GFR"}]}]},
                {"Name": "sim2", "Compounds": [{"Name": "Drug", "Processes": [
                    {"Name": "CYP3A4-1st order CL", "MoleculeName": "CYP3A4"},
                    {"Name": "Glomerular Filtration-GFR",
                     "SystemicProcessType": "GFR"}]}]},
            ],
        }

    def _dangling(self, snap):
        comp = snap["Compounds"][0]
        mols = {p.get("Molecule") for p in comp["Processes"]}
        ds = {(p.get("DataSource") or "").lower() for p in comp["Processes"]}
        bad = []
        for sim in snap["Simulations"]:
            for sc in sim.get("Compounds") or []:
                for p in sc.get("Processes") or []:
                    mn = p.get("MoleculeName")
                    spt = (p.get("SystemicProcessType") or "").lower()
                    if mn is not None and mn not in mols:
                        bad.append(p.get("Name"))
                    elif mn is None and spt and spt not in ds:
                        bad.append(p.get("Name"))
        return bad

    def test_original_has_no_dangling(self):
        self.assertEqual(self._dangling(self.snap), [])

    def test_disable_systemic_removes_sim_refs(self):
        out, rep = apply_edits(self.snap, {"processes": {"GFR": False}})
        # compound lost GFR
        self.assertNotIn("GlomerularFiltration",
                         [p["InternalName"] for p in out["Compounds"][0]["Processes"]])
        # both simulations lost their GFR reference -> no dangling
        self.assertEqual(self._dangling(out), [])
        self.assertEqual(rep["processes"]["GFR"], "disabled")
        self.assertEqual(rep["processes"]["_simulation_refs_removed"], 2)

    def test_disable_enzyme_removes_sim_refs(self):
        out, _ = apply_edits(self.snap, {"processes": {"CYP3A4": False}})
        self.assertEqual(self._dangling(out), [])
        names = [p.get("Name") for sim in out["Simulations"]
                 for sc in sim["Compounds"] for p in sc["Processes"]]
        self.assertNotIn("CYP3A4-1st order CL", names)

    def test_re_add_systemic_restores_sim_refs(self):
        disabled, _ = apply_edits(self.snap, {"processes": {"GFR": False}})
        out, rep = apply_edits(disabled, {"add_processes": [
            {"type": "glomerular_filtration", "parameters": {"GFR fraction": 0.3}}]})
        self.assertEqual(rep["processes"]["add:glomerular_filtration"], "added")
        self.assertEqual(self._dangling(out), [])
        # mirrored to BOTH simulations
        gfr_refs = [p for sim in out["Simulations"] for sc in sim["Compounds"]
                    for p in sc["Processes"]
                    if p.get("SystemicProcessType") == "GFR"]
        self.assertEqual(len(gfr_refs), 2)

    def test_add_enzyme_mirrors_to_sims(self):
        out, rep = apply_edits(self.snap, {"add_processes": [
            {"type": "active_transport_mm", "molecule": "ABCB1"}]})
        self.assertEqual(rep["processes"]["add:ABCB1"], "added")
        self.assertEqual(self._dangling(out), [])
        refs = [p.get("Name") for sim in out["Simulations"]
                for sc in sim["Compounds"] for p in sc["Processes"]
                if p.get("MoleculeName") == "ABCB1"]
        self.assertEqual(len(refs), 2)


if __name__ == "__main__":
    unittest.main()


class TestQualifiedParameters(unittest.TestCase):
    """A parameter name that lives on several processes (per-enzyme clearance)
    can be targeted with '<Name>@<Molecule>' to set only that process."""

    SNAP = {"Compounds": [{"Name": "Drug", "Lipophilicity": [
        {"Parameters": [{"Name": "Lipophilicity", "Value": 2.0}]}], "Processes": [
        {"InternalName": "MetabolizationSpecific_FirstOrder", "Molecule": "UGT1A9",
         "Parameters": [{"Name": "CLspec/[Enzyme]", "Value": 0.1}]},
        {"InternalName": "MetabolizationSpecific_FirstOrder", "Molecule": "UGT2B7",
         "Parameters": [{"Name": "CLspec/[Enzyme]", "Value": 0.1}]}]}]}

    def _cl(self, snap, mol):
        p = next(x for x in snap["Compounds"][0]["Processes"] if x["Molecule"] == mol)
        return next(q["Value"] for q in p["Parameters"] if q["Name"] == "CLspec/[Enzyme]")

    def test_qualified_targets_one_process(self):
        out, rep = apply_edits(self.SNAP, {"parameters": {
            "CLspec/[Enzyme]@UGT1A9": 0.4, "CLspec/[Enzyme]@UGT2B7": 0.007}})
        self.assertEqual(self._cl(out, "UGT1A9"), 0.4)
        self.assertEqual(self._cl(out, "UGT2B7"), 0.007)
        self.assertEqual(rep["not_found"], [])

    def test_unqualified_still_sets_all(self):
        # backward compatible: a plain name sets every matching parameter
        out, _ = apply_edits(self.SNAP, {"parameters": {"CLspec/[Enzyme]": 0.5}})
        self.assertEqual(self._cl(out, "UGT1A9"), 0.5)
        self.assertEqual(self._cl(out, "UGT2B7"), 0.5)

    def test_qualified_not_found_reported(self):
        _, rep = apply_edits(self.SNAP,
                             {"parameters": {"CLspec/[Enzyme]@NOPE": 1.0}})
        self.assertTrue(any("NOPE" in m for m in rep["not_found"]))
