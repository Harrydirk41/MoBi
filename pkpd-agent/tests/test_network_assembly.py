"""Full coupled network: secretion discovery, clearance-name tolerance, joint-steady-state
calibration self-consistency, and single-driver coupling direction. Pure, synthetic fixture."""

import unittest

from pkpd_agent.engines import cell_lifecycle as CL, network_assembly as NA


def _prov(pairs):
    return {n: {"value_from_reference": v, "from_literature": lit} for n, v, lit in pairs}


# a real 2-cytokine / 1-cell loop: K1 makes C1 (modulated by C2) and C2 (modulated by C1);
# C1 drives K1 proliferation. Every free rate is pinned by its own target.
PROV = _prov([
    ("kd_K1_Baseline", 0.10, True), ("kIn_K1_Baseline", 1000.0, True),
    ("K1Prolif_MaxbyC1", 2.0, True),
    ("kcl_C1", 2.0, True), ("kcl_C2", 1.0, True),
    ("C1SecK1_MaxbyC2", 3.0, True), ("C2SecK1_MaxbyC1", 1.5, True),
])
TARGETS = [{"model_species": "K1", "kind": "cell", "target_model_unit": 1e6},
           {"model_species": "C1", "kind": "cytokine", "target_model_unit": 50.0},
           {"model_species": "C2", "kind": "cytokine", "target_model_unit": 10.0}]
LEVELS = {"C1": 50.0, "C2": 10.0}


class TestKclTolerance(unittest.TestCase):
    def test_name_variants(self):
        prov = _prov([("kcl_IL-6", 5.5, True)])
        self.assertEqual(NA._kcl_of(prov, "IL6"), ("kcl_IL-6", 5.5))


class TestSecretion(unittest.TestCase):
    def test_discovers_and_maps_cells(self):
        sec = NA.discover_secretion(PROV, NA.cell_token_map({}))
        self.assertIn("C1", sec)
        self.assertIn("K1", sec["C1"])
        self.assertEqual(sec["C1"]["K1"]["C2"], (3.0, True))

    def test_cell_token_map_inverts_aliases(self):
        self.assertEqual(NA.cell_token_map({"Macrophages": ["Macro", "Macrophase"]}),
                         {"Macro": "Macrophages", "Macrophase": "Macrophages"})


class TestAssembleAndCalibrate(unittest.TestCase):
    def setUp(self):
        cells = CL.discover_cells(PROV, TARGETS, {})
        self.spec, self.meta = NA.assemble_network(PROV, LEVELS, cells, {})
        self.targ = {s["name"]: s["initial"] for s in self.spec["species"]}

    def test_species_present(self):
        names = set(self.targ)
        self.assertEqual(names, {"K1", "C1", "C2"})

    def test_joint_steady_state_self_consistent(self):
        ss = NA.integrate_network(self.spec, t_end=5.0, dt=2e-3)
        for k in self.targ:
            self.assertLess(abs(ss[k] - self.targ[k]) / self.targ[k], 1e-3, f"{k} drifted")

    def test_free_rates_one_per_species(self):
        self.assertEqual(set(self.meta["free_ksec"]), {"ksec_C1", "ksec_C2"})
        self.assertIn("K1", self.meta["free_kprolif"])


class TestCouplingDirection(unittest.TestCase):
    def test_knockdown_drops_downstream_secretion(self):
        cells = CL.discover_cells(PROV, TARGETS, {})
        spec, meta = NA.assemble_network(PROV, LEVELS, cells, {})
        targ = {s["name"]: s["initial"] for s in spec["species"]}
        # C1 secretion is up-modulated by C2 (Max 3.0); knock C2 low -> C1 must fall
        base = NA.integrate_network(spec, clamp={"K1": targ["K1"]}, t_end=10.0, dt=2e-3)
        knock = NA.integrate_network(spec, clamp={"K1": targ["K1"], "C2": targ["C2"] * 0.1},
                                     t_end=10.0, dt=2e-3)
        self.assertLess(knock["C1"], base["C1"] * 0.95)


if __name__ == "__main__":
    unittest.main()
