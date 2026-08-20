"""Readout-mapping benchmark: driver extraction, scoring, loop tools (no LLM/MATLAB)."""

import json
import os
import tempfile
import unittest

from pkpd_agent.engines import ra_readout as RO
from pkpd_agent.tools.registry import ToolRegistry
from pkpd_agent.tools.ra_readout_loop_tools import register_ra_readout_loop_tools


def _write(rules, species):
    fd, path = tempfile.mkstemp(suffix=".json")
    os.close(fd)
    json.dump({"rules": [{"type": "repeatedAssignment", "rule": r} for r in rules],
               "species": [{"name": s} for s in species]}, open(path, "w"))
    return path


class TestDriverExtraction(unittest.TestCase):
    def test_walks_algebraic_graph_to_species(self):
        # DAS28 depends on SCD (a cell sum) and a CRP proxy (off IL6)
        path = _write(
            rules=["DAS28_CRP = 0.5*SCD + 0.3*CRPproxy",
                   "SCD = FLS + Macrophages + Th1",
                   "CRPproxy = 2*IL6"],
            species=["FLS", "Macrophages", "Th1", "IL6", "TNFa"])
        try:
            ext = RO.readout_drivers(path)
            self.assertIn("DAS28_CRP", ext["targets_found"])
            self.assertEqual(set(ext["drivers"]), {"FLS", "Macro", "Th1", "IL6"})
            self.assertNotIn("TNFa", ext["drivers"])   # not referenced by the readout
        finally:
            os.remove(path)

    def test_real_das28_rule_with_compartment_prefix(self):
        # the actual model rule: DAS28 is a Hill-sum of the 9 cell densities (Treg negative)
        das = ("Synovium.DAS28_CRP = 2*FLS^2.5/(1.3E7^2.5+FLS^2.5)+0.5*Endothelial^2.5/"
               "(4.2e7^2.5+Endothelial^2.5)+1.5*Th1^2.5/(4e6^2.5+Th1^2.5)+0.5*Th17^2.5/"
               "(1e5^2.5+Th17^2.5)-0.5*Treg^2.5/(2E6^2.5+Treg^2.5)+0.5*CTL^2.5/"
               "(1.3E5^2.5+CTL^2.5)+1.5*BCells^2.5/(3E6^2.5+BCells^2.5)+1*PlasmaCells^2.5/"
               "(1.8E6^2.5+PlasmaCells^2.5)+2.5*Macrophages^2.5/(2.2e7^2.5+Macrophages^2.5)")
        path = _write(
            rules=["ACR_Perc = 100*delta_DAS28_CRP/DAS28_BL",
                   "delta_DAS28_CRP = DAS28_BL-DAS28_CRP", das],
            species=["FLS", "Endothelial", "Th1", "Th17", "Treg", "CTL", "BCells",
                     "PlasmaCells", "Macrophages", "TNFa", "IL6"])
        try:
            ext = RO.readout_drivers(path)
            self.assertEqual(set(ext["drivers"]),
                             {"FLS", "Endo", "Th1", "Th17", "Treg", "CTL", "BCell",
                              "PlasmaCell", "Macro"})     # the 9 cells, no cytokines
        finally:
            os.remove(path)

    def test_missing_target_reports_rule_names(self):
        path = _write(rules=["SomethingElse = FLS"], species=["FLS"])
        try:
            ext = RO.readout_drivers(path)
            self.assertEqual(ext["targets_found"], [])
            self.assertIn("SomethingElse", ext["all_rule_names"])
        finally:
            os.remove(path)


class TestScore(unittest.TestCase):
    def test_perfect(self):
        s = RO.score_readout(["FLS", "Macrophages", "IL6"], ["FLS", "Macro", "IL6"])
        self.assertEqual(s["recall"], 1.0)
        self.assertEqual(s["precision"], 1.0)

    def test_partial_and_extra(self):
        s = RO.score_readout(["FLS", "TNFa"], ["FLS", "Macro", "IL6"])
        self.assertEqual(s["hit"], 1)
        self.assertEqual(s["recall"], round(1 / 3, 3))
        self.assertEqual(s["precision"], 0.5)          # 1 hit of 2 picks
        self.assertIn("TNFa", s["extra"])              # valid node, but not a driver

    def test_junk_in_extra(self):
        s = RO.score_readout(["aspirinXYZ"], ["FLS"])
        self.assertIn("ASPIRINXYZ", s["extra"])


class _FakeSession:
    def __init__(self):
        self._d = {}

    def get(self, k, default=None):
        return self._d.get(k, default)

    def put(self, k, v):
        self._d[k] = v


class TestLoopTools(unittest.TestCase):
    def _reg(self, drivers):
        reg = ToolRegistry()
        register_ra_readout_loop_tools(reg, None, {"drivers": drivers})
        return reg

    def test_registers_and_hides_key(self):
        reg = self._reg(["FLS", "IL6"])
        for t in ("readout_inspect", "readout_propose", "readout_finalize"):
            self.assertIn(t, reg)
        res = reg.dispatch("readout_inspect", {}, _FakeSession())
        self.assertNotIn("fls", str(res.data).lower())

    def test_propose_and_finalize(self):
        reg = self._reg(["FLS", "Macro", "IL6"])
        sess = _FakeSession()
        reg.dispatch("readout_propose", {"nodes": ["FLS", "macrophage", "IL-6"]}, sess)
        res = reg.dispatch("readout_finalize", {}, sess)
        self.assertTrue(res.ok)
        self.assertEqual(sess.get("readout_final")["recall"], 1.0)

    def test_finalize_requires_proposal(self):
        self.assertFalse(self._reg(["FLS"]).dispatch(
            "readout_finalize", {}, _FakeSession()).ok)


if __name__ == "__main__":
    unittest.main()
