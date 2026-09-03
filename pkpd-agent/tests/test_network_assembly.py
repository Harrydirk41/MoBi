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


class TestApplyStructure(unittest.TestCase):
    def test_model_structure_roundtrips(self):
        cells = CL.discover_cells(PROV, TARGETS, {})
        model_sec = NA.discover_secretion(PROV, NA.cell_token_map({}))
        sec_struct = {cyt: {cell: list(mods) for cell, mods in per.items()}
                      for cyt, per in model_sec.items()}
        cell_struct = {c: {f: list(cells[c][f]) for f in ("prolif", "influx", "apop")}
                       for c in cells}
        sec2, cells2 = NA.apply_structure(model_sec, cells, sec_struct, cell_struct)
        spec, _ = NA.assemble_network(PROV, LEVELS, cells2, {}, sec_override=sec2)
        targ = {s["name"]: s["initial"] for s in spec["species"]}
        ss = NA.integrate_network(spec, t_end=5.0, dt=2e-3)
        for k in targ:
            self.assertLess(abs(ss[k] - targ[k]) / targ[k], 1e-3)

    def test_uncited_agent_edge_uses_prior(self):
        cells = CL.discover_cells(PROV, TARGETS, {})
        model_sec = NA.discover_secretion(PROV, NA.cell_token_map({}))
        # agent adds a NON-model modulator C-fake to C1<-K1: value must be the prior flag (None)
        sec2, _ = NA.apply_structure(model_sec, cells, {"C1": {"K1": ["Cfake"]}}, {})
        self.assertEqual(sec2["C1"]["K1"]["Cfake"], (None, False))


class TestLivePathWithMock(unittest.TestCase):
    def test_full_agent_pipeline_runs_and_calibrates(self):
        # a mock call that echoes every offered candidate back as a regulator (recall high),
        # exercising propose_structure -> apply_structure -> assemble -> integrate end to end
        import json as _json

        def call(system, user):
            if "Candidate cell types:" in user:            # secreting-cells question
                after = user.split("Candidate cell types:", 1)[-1]
                cand = [c.strip() for c in after.split("\n", 1)[0].split(",") if c.strip()]
                return _json.dumps({"cells": cand})
            after = user.split("Available cytokine nodes:", 1)[-1]
            cand = [c.strip() for c in after.split("\n", 1)[0].split(",") if c.strip()]
            return _json.dumps({"regulators": [{"cytokine": c, "direction": "up",
                                "confidence": "high", "basis": "x"} for c in cand]})

        cells = CL.discover_cells(PROV, TARGETS, {})
        sec_struct, cell_struct, scores = NA.propose_structure(PROV, LEVELS, cells, {}, call)
        self.assertIn("secreting_cells", scores)
        sec2, cells2 = NA.apply_structure(
            NA.discover_secretion(PROV, NA.cell_token_map({})), cells, sec_struct, cell_struct)
        spec, _ = NA.assemble_network(PROV, LEVELS, cells2, {}, sec_override=sec2)
        targ = {s["name"]: s["initial"] for s in spec["species"]}
        ss = NA.integrate_network(spec, t_end=5.0, dt=2e-3)      # any structure stays self-consistent
        for k in targ:
            self.assertLess(abs(ss[k] - targ[k]) / targ[k], 1e-3)


class TestPruneStructure(unittest.TestCase):
    def test_drops_low_conf_uncited_keeps_cited_and_high(self):
        cells = CL.discover_cells(PROV, TARGETS, {})
        # agent structure: the real edges + spurious low-confidence uncited extras
        sec_struct = {"C1": {"K1": ["C2", "Cfake"]}}       # C2 cited (keep), Cfake uncited
        cell_struct = {"K1": {"prolif": ["C1", "C2"], "influx": [], "apop": []}}  # C1 cited, C2 not
        conf = {"sec_cell": {"C1": {"K1": "high"}},
                "sec_mod": {"C1": {"C2": "high", "Cfake": "low"}},
                "flux": {"K1": {"prolif": {"C1": "high", "C2": "low"}}}}
        sec2, cell2, dropped = NA.prune_structure(sec_struct, cell_struct, conf, PROV, {}, cells)
        self.assertEqual(sec2["C1"]["K1"], ["C2"])          # Cfake (low+uncited) dropped
        self.assertIn("C1.Cfake", dropped["sec_mod"])
        self.assertEqual(cell2["K1"]["prolif"], ["C1"])     # C2 (low+uncited) dropped
        self.assertIn("K1.prolif.C2", dropped["flux"])

    def test_high_confidence_uncited_edge_is_kept(self):
        cells = CL.discover_cells(PROV, TARGETS, {})
        # an uncited edge the agent is HIGH-confidence about is kept (biology the sparse model pruned)
        cell_struct = {"K1": {"prolif": ["C2"], "influx": [], "apop": []}}
        conf = {"sec_cell": {}, "sec_mod": {}, "flux": {"K1": {"prolif": {"C2": "high"}}}}
        _, cell2, dropped = NA.prune_structure({}, cell_struct, conf, PROV, {}, cells)
        self.assertEqual(cell2["K1"]["prolif"], ["C2"])     # kept on high confidence
        self.assertEqual(dropped["flux"], [])


class TestDivergenceGuard(unittest.TestCase):
    def test_runaway_flags_diverged(self):
        # autocatalysis with no clearance: dX/dt = X -> unbounded -> must flag, not return garbage
        spec = {"name": "u", "species": [{"name": "X", "initial": 1.0}],
                "parameters": [{"name": "k", "value": 1.0}],
                "reactions": [{"id": "g", "reactants": [], "products": ["X"], "rate": "k * X"}],
                "rules": []}
        out = NA.integrate_network(spec, t_end=100.0, dt=1e-2, diverge_fold=1e6)
        self.assertTrue(out.get("__diverged__"))

    def test_stable_has_no_flag(self):
        cells = CL.discover_cells(PROV, TARGETS, {})
        spec, _ = NA.assemble_network(PROV, LEVELS, cells, {})
        out = NA.integrate_network(spec, t_end=5.0, dt=2e-3)
        self.assertNotIn("__diverged__", out)


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
