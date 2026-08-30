"""Model assembly: infix->MathML, combined-effect motif, SBML emit round-trips through the
importer, subsystem build. Pure - validated by re-parsing the emitted SBML with sbml_import."""

import unittest
import xml.etree.ElementTree as ET

from pkpd_agent.engines import model_assembly as MA
from pkpd_agent.engines.sbml_import import mathml_to_infix, sbml_to_network


def _mathml_roundtrip(expr):
    ml = MA.infix_to_mathml(expr)
    return mathml_to_infix(ET.fromstring(
        f'<math xmlns="http://www.w3.org/1998/Math/MathML">{ml}</math>'))


class TestInfixMathml(unittest.TestCase):
    def test_precedence_and_functions(self):
        self.assertEqual(_mathml_roundtrip("a + b * c ^ 2"), "(a + (b * (c ^ 2)))")
        self.assertEqual(_mathml_roundtrip("MM(TNFa, k)"), "MM(TNFa,k)")

    def test_division_and_parens(self):
        self.assertEqual(_mathml_roundtrip("(Max - 1) * X / (K + X)"),
                         "(((Max - 1) * X) / (K + X))")


class TestCombinedEffect(unittest.TestCase):
    def test_motif_multiplies_fold_changes(self):
        r = MA.combined_effect("kbase", [{"species": "IL6", "max_param": "M6", "k_param": "K6"},
                                         {"species": "TNFa", "max_param": "MT", "k_param": "KT"}])
        self.assertIn("kbase", r)
        self.assertIn("(1 + (M6 - 1) * IL6 / (K6 + IL6))", r)
        self.assertIn("(1 + (MT - 1) * TNFa / (KT + TNFa))", r)


class TestAssembleRoundTrip(unittest.TestCase):
    def test_subsystem_emits_valid_roundtripping_sbml(self):
        spec = MA.build_subsystem(
            "FLS", "kg", "kd",
            [{"species": "IL6", "max_param": "M6", "k_param": "K6"}],
            values={"kg": 0.5, "kd": 0.1, "M6": 2.0, "K6": 50, "FLS_init": 1e6},
            clamp={"IL6": 100.0})
        net = sbml_to_network_from_str(MA.to_sbml(spec))
        names = {s["name"] for s in net["species"]}
        self.assertEqual(names, {"FLS", "IL6"})
        self.assertEqual(len(net["parameters"]), 4)     # kg, kd, M6, K6 (FLS_init excluded)
        rxn = {(r.get("reaction") or r.get("name")): r for r in net["reactions"]}
        # proliferation produces FLS and its rate carries the IL6 fold-change term
        prolif = next(r for r in net["reactions"] if r.get("products") == ["FLS"])
        self.assertIn("M6", prolif["rate"])
        self.assertIn("IL6", prolif["rate"])


class TestMotif(unittest.TestCase):
    REGS = [{"species": "IL6", "max_param": "M6", "k_param": "K6"},
            {"species": "TNFa", "max_param": "MT", "k_param": "KT"}]

    def test_propose_motif_parses(self):
        call = lambda s, u: ('{"proliferation_order":"zeroth","combination":"capped_sum",'
                             '"cap":10,"per_regulator":"hill","reason":"saturating sum"}')
        m = MA.propose_motif("FLS", self.REGS, "ref rate", call)
        self.assertEqual(m["proliferation_order"], "zeroth")
        self.assertEqual(m["combination"], "capped_sum")
        self.assertEqual(m["cap"], 10)

    def test_rate_capped_sum_form(self):
        m = {"proliferation_order": "zeroth", "combination": "capped_sum", "cap": 10}
        r = MA.rate_from_motif(m, "kg", self.REGS, "FLS")
        self.assertIn("min(10, 1 +", r)
        self.assertIn("(M6 - 1) * IL6 / (K6 + IL6)", r)
        self.assertNotIn("* FLS", r)                   # zeroth order

    def test_rate_first_order_product(self):
        m = {"proliferation_order": "first", "combination": "product", "cap": None}
        r = MA.rate_from_motif(m, "kg", self.REGS, "FLS")
        self.assertTrue(r.rstrip().endswith("* FLS"))  # first order
        self.assertIn(") * (1 +", r)                   # product of folds


def sbml_to_network_from_str(xml: str):
    import os, tempfile
    p = os.path.join(tempfile.gettempdir(), "asm_test.xml")
    with open(p, "w", encoding="utf-8") as fh:
        fh.write(xml)
    return sbml_to_network(p)


if __name__ == "__main__":
    unittest.main()
