"""Template A - the cell life-cycle engine: discovery, base-rate fit, 3-flux assembly, and
calibration self-consistency (the assembled cell holds at its target). Pure, synthetic fixture."""

import os
import tempfile
import unittest

from pkpd_agent.engines import cell_lifecycle as CL, model_assembly as MA
from pkpd_agent.engines.sbml_import import sbml_to_network
from examples.run_qsp_end_to_end import integrate


def _prov(pairs):
    return {n: {"value_from_reference": v, "from_literature": lit} for n, v, lit in pairs}


# one cell with influx (non-marginal) + one without (marginal), two cytokine regulators
PROV = _prov([
    ("kd_K1_Baseline", 0.10, True), ("kIn_K1_Baseline", 1000.0, True),
    ("K1Prolif_MaxbyC1", 2.0, True), ("K1Apop_MaxbyC2", 0.5, True),
    ("kd_K2_Baseline", 0.20, True),                      # K2: no kIn -> marginal
    ("K2Prolif_MaxbyC1", 3.0, True),
])
TARGETS = [{"model_species": "K1", "kind": "cell", "target_model_unit": 1e6},
           {"model_species": "K2", "kind": "cell", "target_model_unit": 5e5},
           {"model_species": "C1", "kind": "cytokine", "target_model_unit": 40.0},
           {"model_species": "C2", "kind": "cytokine", "target_model_unit": 8.0}]
LEVELS = {"C1": 40.0, "C2": 8.0}


class TestDiscover(unittest.TestCase):
    def setUp(self):
        self.cells = CL.discover_cells(PROV, TARGETS, {})

    def test_flux_sets_and_baselines(self):
        k1 = self.cells["K1"]
        self.assertEqual(set(k1["prolif"]), {"C1"})
        self.assertEqual(set(k1["apop"]), {"C2"})
        self.assertEqual(k1["influx"], {})
        self.assertEqual(k1["kd_val"], 0.10)
        self.assertEqual(k1["kin_val"], 1000.0)
        self.assertEqual(k1["target"], 1e6)

    def test_marginal_vs_influx(self):
        _, m1 = CL.fit_base_prolif(self.cells["K1"], LEVELS)
        _, m2 = CL.fit_base_prolif(self.cells["K2"], LEVELS)
        self.assertFalse(m1)          # K1 has influx -> level-pinned
        self.assertTrue(m2)           # K2 has none -> marginal


class TestFit(unittest.TestCase):
    def test_marginal_balances_birth_and_death(self):
        cells = CL.discover_cells(PROV, TARGETS, {})
        info = cells["K2"]
        kp, marg = CL.fit_base_prolif(info, LEVELS)
        # marginal: kprolif = kd * apop_eff / prolif_eff, with prolif_eff = 1+(3-1)*0.5 = 2, apop=1
        self.assertTrue(marg)
        self.assertAlmostEqual(kp, 0.20 * 1.0 / 2.0, places=9)


class TestAssembleHoldsTarget(unittest.TestCase):
    def _build_and_hold(self, cell):
        cells = CL.discover_cells(PROV, TARGETS, {})
        info = cells[cell]
        kpp = f"kprolif_{cell}"
        kp, _ = CL.fit_base_prolif(info, LEVELS)
        rxns, vals = CL.cell_reactions(cell, info, LEVELS, kpp)
        vals[kpp] = kp
        regs = sorted(c for c in CL.all_regulators(info) if c in LEVELS)
        species = [{"name": cell, "initial": info["target"]}] + \
                  [{"name": c, "initial": LEVELS[c], "boundary": True} for c in regs]
        spec = {"name": cell, "species": species,
                "parameters": [{"name": k, "value": v} for k, v in vals.items()],
                "reactions": rxns, "rules": []}
        xml = os.path.join(tempfile.gettempdir(), f"tc_{cell}.xml")
        with open(xml, "w", encoding="utf-8") as f:
            f.write(MA.to_sbml(spec))
        net = sbml_to_network(xml)
        clamp = {c: LEVELS[c] for c in regs}
        return integrate(net, clamp, t_end=50.0, dt=5e-3)[cell], info["target"]

    def test_influx_cell_holds(self):
        got, target = self._build_and_hold("K1")
        self.assertLess(abs(got - target) / target, 1e-3)

    def test_marginal_cell_holds(self):
        got, target = self._build_and_hold("K2")
        self.assertLess(abs(got - target) / target, 1e-3)


class TestRoundTripInitial(unittest.TestCase):
    def test_sbml_preserves_initial_and_boundary(self):
        spec = {"name": "m", "species": [{"name": "X", "initial": 1234.0},
                {"name": "Y", "initial": 7.0, "boundary": True}],
                "parameters": [{"name": "k", "value": 1.0}], "reactions": [], "rules": []}
        xml = os.path.join(tempfile.gettempdir(), "tc_rt.xml")
        with open(xml, "w", encoding="utf-8") as f:
            f.write(MA.to_sbml(spec))
        net = sbml_to_network(xml)
        byname = {s["name"]: s for s in net["species"]}
        self.assertEqual(byname["X"]["initial"], 1234.0)
        self.assertTrue(byname["Y"].get("boundary"))


if __name__ == "__main__":
    unittest.main()
