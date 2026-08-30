"""Stage-1 topology reconstruction: ground-truth edge extraction, drafting, scoring.

Pure functions only - a STUB call_fn stands in for the drafter LLM, and a small hand-built
network.json dict stands in for the model dump. No MATLAB, no API. This proves the answer-key
extraction and the precision/recall scoring; the real LLM is swapped in only in the runner.
"""

import unittest

from pkpd_agent.engines import llm_topology as LT


# a tiny fake network.json: TNFa drives FLS proliferation via a regulatory intermediate,
# IL6 modifies the same reaction directly, and FLS mass-converts into an activated form.
_NETWORK = {
    "species": [{"name": n} for n in ["TNFa", "IL6", "FLS", "FLS_active", "Drug"]],
    "rules": [
        {"rule": "Pro_FLSProlif_byTNFa_effect = MM(TNFa, k1)"},
    ],
    "reactions": [
        # FLS proliferation: rate uses the TNFa-defined intermediate + IL6 directly
        {"reaction": "-> FLS", "rate": "Pro_FLSProlif_byTNFa_effect * IL6 * kbase",
         "reactants": [], "products": ["FLS"]},
        # activation: mass flow FLS -> FLS_active, modified by Drug (non-species-free)
        {"reaction": "FLS -> FLS_active", "rate": "kact * FLS",
         "reactants": ["FLS"], "products": ["FLS_active"]},
    ],
    "parameters": [{"name": "k1"}, {"name": "kbase"}, {"name": "kact"}],
}


class TestGroundTruthEdges(unittest.TestCase):
    def setUp(self):
        self.edges = LT.ground_truth_edges(_NETWORK)

    def test_intermediate_param_expands_to_source_species(self):
        # TNFa influences FLS via Pro_FLSProlif_byTNFa_effect
        self.assertIn(("TNFa", "FLS"), self.edges)

    def test_direct_modifier_species_in_rate(self):
        self.assertIn(("IL6", "FLS"), self.edges)

    def test_mass_flow_reactant_to_product(self):
        self.assertIn(("FLS", "FLS_active"), self.edges)

    def test_no_self_edges(self):
        self.assertFalse(any(s == d for s, d in self.edges))

    def test_drug_not_a_species_is_ignored(self):
        # Drug appears in no rate here; even if it did, only species-list names count
        self.assertNotIn("Drug", {s for s, _ in self.edges} | {d for _, d in self.edges})


class TestDraftTopology(unittest.TestCase):
    def test_filters_to_real_node_pairs(self):
        call = lambda s, u: (
            '{"edges": ['
            '{"src": "TNFa", "dst": "FLS", "sign": "activate", "basis": "ref 1"},'
            '{"src": "TNFa", "dst": "Nonexistent", "sign": "activate", "basis": "x"},'
            '{"src": "FLS", "dst": "FLS", "sign": "activate", "basis": "self"}]}'
        )
        out = LT.draft_topology(["TNFa", "FLS", "IL6"], "refs...", call)
        self.assertEqual(len(out), 1)
        self.assertEqual(out[0]["src"], "TNFa")
        self.assertEqual(out[0]["dst"], "FLS")
        self.assertEqual(out[0]["sign"], "activate")


class TestCompareTopology(unittest.TestCase):
    def test_perfect_draft(self):
        truth = {("TNFa", "FLS"), ("IL6", "FLS")}
        draft = [{"src": "TNFa", "dst": "FLS"}, {"src": "IL6", "dst": "FLS"}]
        r = LT.compare_topology(draft, truth)
        self.assertEqual(r["recall"], 1.0)
        self.assertEqual(r["precision"], 1.0)
        self.assertEqual(r["f1"], 1.0)

    def test_partial_with_extra(self):
        truth = {("TNFa", "FLS"), ("IL6", "FLS")}
        draft = [{"src": "TNFa", "dst": "FLS"}, {"src": "IL6", "dst": "TNFa"}]
        r = LT.compare_topology(draft, truth)
        self.assertEqual(r["hit"], 1)
        self.assertEqual(r["recall"], 0.5)
        self.assertEqual(r["precision"], 0.5)
        self.assertIn(("IL6", "TNFa"), r["extra"])
        self.assertIn(("IL6", "FLS"), r["missed"])

    def test_empty_draft(self):
        r = LT.compare_topology([], {("A", "B")})
        self.assertEqual(r["precision"], 0.0)
        self.assertEqual(r["recall"], 0.0)
        self.assertEqual(r["f1"], 0.0)


if __name__ == "__main__":
    unittest.main()
