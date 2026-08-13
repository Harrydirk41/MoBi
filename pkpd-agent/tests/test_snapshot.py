import json
import os
import tempfile
import unittest

from pkpd_agent.engines.osp_snapshot import OSPSnapshot

# A minimal snapshot mimicking the OSP structure (µmol/l observed data).
FIXTURE = {
    "Version": 1,
    "Compounds": [{
        "Name": "TestDrug",
        "Lipophilicity": [{"Name": "Measurement"}],
        "IntestinalPermeability": [{"Name": "fit"}],
        "CalculationMethods": ["Cellular partition coefficient method - Rodgers and Rowland"],
        "Processes": [{"Molecule": "CYP3A4", "DataSource": "fitted"}],
        "Parameters": [{"Name": "Molecular weight", "Value": 400.0, "Unit": "g/mol"}],
        "Lipo": {"Name": "Lipophilicity", "Value": 2.0, "Unit": "Log Units"},
    }],
    "Simulations": [{"Name": "sim1", "Model": "4Comp"}],
    "Protocols": [{"Name": "p"}],
    "ObservedData": [{
        "Name": "Study A, 100 mg po",
        "Columns": [{
            "Name": "Concentration",
            "DataInfo": {"MolWeight": 400.0},
            "Values": [0.0, 2.5, 1.25, 0.6],
            "Unit": "µmol/l",
        }],
        "BaseGrid": {"Name": "Time", "Values": [0.0, 1.0, 4.0, 8.0], "Unit": "h"},
    }],
}


class TestSnapshotExtract(unittest.TestCase):
    def setUp(self):
        self.snap = OSPSnapshot(FIXTURE)

    def test_summary(self):
        s = self.snap.summary()
        self.assertEqual(s["compound"], "TestDrug")
        self.assertEqual(s["n_observed_datasets"], 1)

    def test_observed_profile_and_unit_conversion(self):
        p = self.snap.observed_profiles()[0]
        self.assertEqual(p.time_h, [0.0, 1.0, 4.0, 8.0])
        # µmol/L -> mg/L uses MW/1000: 2.5 µmol/L * 400/1000 = 1.0 mg/L
        self.assertAlmostEqual(p.conc_mg_L()[1], 1.0, places=6)

    def test_nca(self):
        row = self.snap.nca_table()[0]
        self.assertEqual(row["study"], "Study A, 100 mg po")
        self.assertAlmostEqual(row["c_max_mg_L"], 1.0, places=6)  # peak 2.5 µmol/L
        self.assertGreater(row["auc_mg_h_L"], 0)

    def test_modeling_choices_answer_key(self):
        mc = self.snap.modeling_choices()
        self.assertIn("Rodgers and Rowland", mc["calculation_methods"][0])
        self.assertEqual(mc["fit_vs_measured"]["IntestinalPermeability"], "fit")
        self.assertEqual(mc["fit_vs_measured"]["Lipophilicity"], "Measurement")
        self.assertEqual(mc["processes"][0]["molecule"], "CYP3A4")
        self.assertEqual(mc["model_type"], "4Comp")

    def test_compound_parameters(self):
        params = {p["parameter"]: p["value"] for p in self.snap.compound_parameters()}
        self.assertIn("Molecular weight", params)
        self.assertEqual(params["Molecular weight"], 400.0)

    def test_write_csvs(self):
        with tempfile.TemporaryDirectory() as d:
            paths = self.snap.write_csvs(d)
            self.assertTrue(os.path.exists(paths["observed_csv"]))
            self.assertTrue(os.path.exists(paths["params_csv"]))
            with open(paths["observed_csv"]) as fh:
                self.assertIn("conc_mg_L", fh.readline())

    def test_tool_registered_and_dispatches(self):
        from pkpd_agent.config import AgentConfig
        from pkpd_agent.tools import build_default_registry
        from pkpd_agent.state import ModelingSession
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "snap.json")
            json.dump(FIXTURE, open(p, "w"))
            reg = build_default_registry(AgentConfig())
            self.assertIn("snapshot_extract", reg)
            r = reg.dispatch("snapshot_extract", {"path": p}, ModelingSession(goal="t"))
            self.assertTrue(r.ok)
            self.assertEqual(r.data["compound"], "TestDrug")


if __name__ == "__main__":
    unittest.main()
