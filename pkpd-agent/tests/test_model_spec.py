"""Layered model spec: build_spec fills the four layers with provenance tags; to_subsystem
round-trips to runnable SBML. Pure."""

import unittest

from pkpd_agent.engines import model_spec as MS
from pkpd_agent.engines.sbml_import import sbml_to_network
from pkpd_agent.engines import model_assembly as MA


PROV = {
    "IL6SecFLS_MaxbyIL1b": {"value_from_reference": 69.0, "from_literature": True, "reference": "Fig1A"},
    "IL6SecFLS_MaxbyIL17": {"value_from_reference": 4.0, "from_literature": True, "reference": "X"},
    "kcl_IL6": {"value_from_reference": 5.5, "from_literature": True, "reference": "Y"},
}
LEVELS = {"IL6": 289.0, "IL1b": 1.2, "IL17": 0.007, "VEGF": 118.0}


class TestBuildSpec(unittest.TestCase):
    def setUp(self):
        self.chosen = [{"cytokine": "IL1b", "direction": "up", "basis": "b"},
                       {"cytokine": "VEGF", "direction": "up"}]        # VEGF = over-inclusion
        self.motif = {"proliferation_order": "zeroth", "combination": "product", "cap": None}
        self.spec = MS.build_spec({"objective": "o"}, "IL6", self.chosen, self.motif,
                                  PROV, LEVELS, truth_regulators={"IL1b", "IL17"})

    def test_edges_tagged_and_verified(self):
        e = {x["src"]: x for x in self.spec["edges"]}
        self.assertEqual(e["IL1b"]["source"], MS.AGENT)
        self.assertTrue(e["IL1b"]["verify"]["in_model"])       # real regulator
        self.assertFalse(e["VEGF"]["verify"]["in_model"])      # over-inclusion flagged

    def test_constants_provenance(self):
        c = {x["param"]: x for x in self.spec["constants"]}
        self.assertEqual(c["M_IL1b"]["provenance"], MS.LITERATURE)  # cited
        self.assertEqual(c["M_VEGF"]["provenance"], MS.PRIOR)       # no citation -> prior
        self.assertEqual(c["K_IL1b"]["provenance"], MS.FIT)         # half-effect always fit
        self.assertFalse(c["K_IL1b"]["identifiable"])
        self.assertEqual(c["K_IL1b"]["needs_data"], "dose_response")
        self.assertTrue(c["kg_IL6"]["identifiable"])               # baseline pinned by steady state

    def test_form_capped_sum_unresponsive_flagged(self):
        # a binding cap (IL1b excess huge) makes the hub unresponsive - verify catches it
        spec = MS.build_spec({}, "IL6", [{"cytokine": "IL1b", "direction": "up"}],
                             {"proliferation_order": "zeroth", "combination": "capped_sum",
                              "cap": 2}, PROV, LEVELS, truth_regulators={"IL1b"})
        self.assertFalse(spec["forms"][0]["verify"]["responds_to_single_knockdown"])

    def test_rollup_counts(self):
        r = MS.provenance_rollup(self.spec)
        self.assertEqual(r["edges_matching_model"], 1)             # only IL1b
        self.assertEqual(r["edges_total"], 2)
        self.assertGreaterEqual(r["constants_by_source"][MS.FIT], 2)

    def test_to_subsystem_roundtrips_and_hits_target(self):
        sub = MS.to_subsystem(self.spec, PROV, LEVELS)
        net = sbml_to_network_str(MA.to_sbml(sub))
        self.assertIn("IL6", {s["name"] for s in net["species"]})
        prolif = next(r for r in net["reactions"] if r.get("products") == ["IL6"])
        self.assertIn("IL1b", prolif["rate"])


def sbml_to_network_str(xml):
    import os, tempfile
    p = os.path.join(tempfile.gettempdir(), "spec_test.xml")
    open(p, "w", encoding="utf-8").write(xml)
    return sbml_to_network(p)


if __name__ == "__main__":
    unittest.main()
