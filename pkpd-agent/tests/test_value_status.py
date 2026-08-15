"""Parameter value-status classification (given vs placeholder vs structural).

The agent must be able to tell a real measured input from a blanked benchmark
placeholder, so it does not trust a naive default (e.g. solubility=100) as if it
were measured. This is problem-setup information, never the answer value.
"""

import os
import unittest

from pkpd_agent.tools.osp_loop_tools import _value_status, _current_model


class TestValueStatus(unittest.TestCase):
    def test_published_is_given(self):
        self.assertEqual(_value_status(
            {"Source": "Publication", "Description": "Smith 2014"}), "given")
        self.assertEqual(_value_status({"Source": "Database"}), "given")
        self.assertEqual(_value_status({"Source": "In Vitro"}), "given")

    def test_blanked_marker_is_placeholder(self):
        self.assertEqual(_value_status(
            {"Source": "Unknown",
             "Description": "benchmark naive prior (blanked)"}), "placeholder")

    def test_untagged_constant_is_structural_not_placeholder(self):
        # molecular weight with no ValueOrigin is a real constant, NOT a blanked
        # unknown - it must not be flagged "determine this".
        self.assertEqual(_value_status(None), "structural")
        self.assertEqual(_value_status({}), "structural")

    def test_unknown_without_blank_marker_is_structural(self):
        # an 'Unknown' source that is not the benchmark blank marker is structural
        self.assertEqual(_value_status(
            {"Source": "Unknown", "Description": "model default"}), "structural")


BLANKED = os.path.join(
    os.path.dirname(__file__), "..", "..", "OSP-PBPK-Model-Library",
    "Dapagliflozin", "benchmark", "Dapagliflozin-Model.blanked.json")


@unittest.skipUnless(os.path.exists(BLANKED), "Dapagliflozin blanked not present")
class TestCurrentModelStatus(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.m = _current_model(BLANKED)
        cls.by = {p["name"]: p.get("value_status") for p in cls.m["parameters"]}

    def test_measured_inputs_are_given(self):
        self.assertEqual(self.by["Fraction unbound (plasma, reference value)"], "given")
        self.assertEqual(self.by["Reference pH"], "given")

    def test_blanked_unknowns_are_placeholders(self):
        for n in ("Lipophilicity", "Solubility at reference pH",
                  "CLspec/[Enzyme]@UGT1A9"):
            self.assertEqual(self.by[n], "placeholder", n)

    def test_molecular_weight_is_not_a_placeholder(self):
        # a real physical constant present as-is - must not be "determine this"
        self.assertEqual(self.by["Molecular weight"], "structural")


if __name__ == "__main__":
    unittest.main()
