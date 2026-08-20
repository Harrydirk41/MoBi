"""Stage-1 network-reconstruction: truth parsing, scoring, and loop-tool registration.

All synthetic (no MATLAB, no LLM): the pure-Python layer that parses the model's
regulatory edges from its naming convention, scores a proposed network against the
truth, and drives the propose/finalize loop.
"""

import os
import unittest

from pkpd_agent.engines import ra_network as N
from pkpd_agent.tools.registry import ToolRegistry
from pkpd_agent.tools.ra_network_loop_tools import register_ra_network_loop_tools

_SBPROJ = os.path.abspath(os.path.join(
    os.path.dirname(__file__), "..", "..", "RA-QSP-Model",
    "Vantage RA QSP Model v1.0.sbproj"))


class TestNodeCanon(unittest.TestCase):
    def test_aliases(self):
        self.assertEqual(N.canon_node("macrophage"), "Macro")
        self.assertEqual(N.canon_node("TNF-a"), "TNFa")
        self.assertEqual(N.canon_node("IFNgamma"), "IFNg")
        self.assertEqual(N.canon_node("FLS"), "FLS")

    def test_unknown_is_none(self):
        self.assertIsNone(N.canon_node("Aspirin"))
        self.assertIsNone(N.canon_node(""))


class TestEdgeParsing(unittest.TestCase):
    def test_pro_edge(self):
        e = N.edges_from_names(["Pro_IL6Sec_byMacro_effect"])
        self.assertEqual(len(e), 1)
        self.assertEqual(e[0].signed(), ("Macro", 1, "IL6"))

    def test_anti_edge_is_negative(self):
        e = N.edges_from_names(["Anti_TNFaSec_byFLS_effect"])
        self.assertEqual(e[0].signed(), ("FLS", -1, "TNFa"))

    def test_hill_edge_strips_producer_cell(self):
        # GMCSF secreted by Macro, DRIVEN by TNFa -> regulator is TNFa
        e = N.edges_from_names(["Hill_GMCSFSecMacro_byTNFa"])
        self.assertEqual(e[0].signed(), ("TNFa", 1, "GMCSF"))

    def test_unknown_nodes_dropped(self):
        self.assertEqual(N.edges_from_names(["Pro_XSec_byAspirin"]), [])

    def test_self_loop_dropped(self):
        self.assertEqual(N.edges_from_names(["Pro_IL6Sec_byIL6"]), [])

    def test_dedupe(self):
        e = N.edges_from_names(["Pro_IL6Sec_byMacro", "Pro_IL6Sec_byMacro_effect"])
        self.assertEqual(len(e), 1)


class TestRuleParsing(unittest.TestCase):
    def test_mm_terms_become_edges(self):
        rules = [{"rule": "Pro_FLSProlif_effect = min(10,MM(TNFa,a,b,c)+MM(IL6,a,b,c))"}]
        e = {x.signed() for x in N.edges_from_rules(rules)}
        self.assertIn(("TNFa", 1, "FLS"), e)
        self.assertIn(("IL6", 1, "FLS"), e)

    def test_anti_rule_negative(self):
        rules = [{"rule": "Anti_EndoInflux_effect = min(0.9,MM(TGFb,a,b,c))"}]
        self.assertEqual(N.edges_from_rules(rules)[0].signed(), ("TGFb", -1, "Endo"))

    def test_non_regulatory_rule_ignored(self):
        rules = [{"rule": "SCD = FLS+Endothelial+Macrophages+BCells"}]
        self.assertEqual(N.edges_from_rules(rules), [])

    def test_new_nodes_recognized(self):
        self.assertEqual(N.canon_node("TGFb"), "TGFb")
        self.assertEqual(N.canon_node("IL10"), "IL10")
        self.assertEqual(N.canon_node("Macrophages"), "Macro")
        self.assertEqual(N.canon_node("PlasmaCells"), "PlasmaCell")


class TestScoring(unittest.TestCase):
    def setUp(self):
        self.truth = [N.Edge("Macro", 1, "IL6"), N.Edge("Macro", 1, "TNFa"),
                      N.Edge("Treg", -1, "TNFa"), N.Edge("Th17", 1, "IL17")]

    def test_perfect(self):
        s = N.score_network(self.truth, self.truth)
        self.assertEqual((s["precision"], s["recall"], s["f1"]), (1.0, 1.0, 1.0))

    def test_partial_and_extra(self):
        prop = [N.Edge("Macro", 1, "IL6"), N.Edge("Th1", 1, "IFNg")]  # 1 hit, 1 extra
        s = N.score_network(prop, self.truth)
        self.assertEqual(s["hit"], 1)
        self.assertEqual(s["extra"], 1)
        self.assertEqual(s["missed"], 3)
        self.assertAlmostEqual(s["precision"], 0.5)
        self.assertAlmostEqual(s["recall"], 0.25)

    def test_sign_matters_when_sign_aware(self):
        prop = [N.Edge("Treg", 1, "TNFa")]      # right pair, wrong sign
        self.assertEqual(N.score_network(prop, self.truth, sign_aware=True)["hit"], 0)
        self.assertEqual(N.score_network(prop, self.truth, sign_aware=False)["hit"], 1)

    def test_empty_proposal(self):
        s = N.score_network([], self.truth)
        self.assertEqual((s["precision"], s["recall"]), (0.0, 0.0))


class TestProposalParsing(unittest.TestCase):
    def test_string_signs(self):
        items = [{"source": "macrophage", "target": "IL6", "sign": "promote"},
                 {"source": "Treg", "target": "TNFa", "sign": "inhibit"}]
        e = {x.signed() for x in N.edges_from_proposal(items)}
        self.assertIn(("Macro", 1, "IL6"), e)
        self.assertIn(("Treg", -1, "TNFa"), e)

    def test_bad_nodes_skipped(self):
        self.assertEqual(N.edges_from_proposal([{"source": "x", "target": "y"}]), [])


class TestDiagramTruth(unittest.TestCase):
    @unittest.skipUnless(os.path.exists(_SBPROJ), "sbproj not present")
    def test_diagram_yields_edges(self):
        truth = N.parse_truth_from_diagram(_SBPROJ)
        self.assertGreater(len(truth), 15)
        # every edge is over the known cast
        for e in truth:
            self.assertIn(e.source, N.NODES)
            self.assertIn(e.target, N.NODES)


class TestSignScoring(unittest.TestCase):
    def setUp(self):
        self.truth = [N.Edge("Macro", 1, "IL6"), N.Edge("Macro", 1, "TNFa"),
                      N.Edge("Treg", -1, "TNFa"), N.Edge("TGFb", -1, "BCell")]

    def test_perfect(self):
        pred = {e.pair(): e.sign for e in self.truth}
        s = N.score_signs(pred, self.truth)
        self.assertEqual(s["accuracy"], 1.0)

    def test_majority_baseline(self):
        # 2 of 4 positive -> majority baseline 0.5
        s = N.score_signs({}, self.truth)
        self.assertEqual(s["majority_baseline"], 0.5)
        self.assertEqual(s["correct"], 0)

    def test_all_positive_guess(self):
        pred = {e.pair(): 1 for e in self.truth}   # guess all activate
        s = N.score_signs(pred, self.truth)
        self.assertEqual(s["correct"], 2)          # the 2 real positives
        self.assertFalse(s["beats_majority"])      # equals majority, not beats


class _FakeSession:
    def __init__(self):
        self._d = {}

    def get(self, k, default=None):
        return self._d.get(k, default)

    def put(self, k, v):
        self._d[k] = v


class TestLoopTools(unittest.TestCase):
    def _reg(self, truth):
        reg = ToolRegistry()
        register_ra_network_loop_tools(reg, None, {"truth": truth})
        return reg

    def test_registers(self):
        reg = self._reg([])
        for t in ("network_inspect", "network_propose", "network_finalize"):
            self.assertIn(t, reg)

    def test_inspect_hides_key(self):
        reg = self._reg([N.Edge("Macro", 1, "IL6")])
        res = reg.dispatch("network_inspect", {}, _FakeSession())
        self.assertTrue(res.ok)
        self.assertIn("cells", res.data)
        self.assertIn("cytokines", res.data)
        # the answer key must not leak
        blob = str(res.data).lower()
        self.assertNotIn("truth", blob)

    def test_propose_accumulates_and_dedupes(self):
        reg = self._reg([N.Edge("Macro", 1, "IL6")])
        sess = _FakeSession()
        reg.dispatch("network_propose", {"edges": [
            {"source": "Macro", "target": "IL6", "sign": 1}]}, sess)
        res = reg.dispatch("network_propose", {"edges": [
            {"source": "Macro", "target": "IL6", "sign": 1},   # dup
            {"source": "Th17", "target": "IL17", "sign": 1}]}, sess)
        self.assertEqual(res.data["total_edges"], 2)

    def test_finalize_requires_edges(self):
        reg = self._reg([N.Edge("Macro", 1, "IL6")])
        res = reg.dispatch("network_finalize", {}, _FakeSession())
        self.assertFalse(res.ok)

    def test_finalize_scores_and_commits(self):
        truth = [N.Edge("Macro", 1, "IL6"), N.Edge("Th17", 1, "IL17")]
        reg = self._reg(truth)
        sess = _FakeSession()
        reg.dispatch("network_propose", {"edges": [
            {"source": "Macro", "target": "IL6", "sign": 1}]}, sess)
        res = reg.dispatch("network_finalize", {}, sess)
        self.assertTrue(res.ok)
        self.assertEqual(res.data["topology"]["hit"], 1)
        self.assertEqual(sess.get("net_final")["topology"]["n_truth"], 2)


if __name__ == "__main__":
    unittest.main()
