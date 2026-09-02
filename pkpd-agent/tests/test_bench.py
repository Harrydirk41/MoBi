"""Topology+form benchmark scoring: recall/precision/high-conf-extra and form matching. Pure."""

import unittest

from examples import run_qsp_bench as B

PROV = {
    "IL6SecFLS_MaxbyIL1b": {}, "IL6SecFLS_MaxbyIL17": {}, "IL6SecMacro_MaxbyTNFa": {},
    "RANTESSecFLS_MaxbyIL1b": {}, "RANTESSecFLS_byIL1b_MaxbyIFNg": {},   # nested modifier
    "TNFaSecMacro_MaxbyIL17": {}, "TNFaSecMacro_MaxbyAutoAb": {},        # non-cytokine regulator
}


class TestTruth(unittest.TestCase):
    def test_full_regulators_catches_nested(self):
        self.assertEqual(B.full_regulators(PROV, "RANTES"), {"IL1b", "IFNg"})   # nested IFNg caught

    def test_full_regulators_includes_noncandidate(self):
        self.assertEqual(B.full_regulators(PROV, "TNFa"), {"IL17", "AutoAb"})   # AutoAb present here


class TestScoreTopology(unittest.TestCase):
    def test_recall_precision_and_high_conf_extra(self):
        truth = {"IL1b", "IL17", "TNFa"}
        chosen = [{"cytokine": "IL1b", "confidence": "high"},
                  {"cytokine": "TNFa", "confidence": "high"},
                  {"cytokine": "IL17", "confidence": "high"},
                  {"cytokine": "IL6", "confidence": "high"},        # high-conf over-inclusion
                  {"cytokine": "IL10", "confidence": "low"}]
        s = B.score_topology(chosen, truth)
        self.assertEqual(s["recall"], 1.0)
        self.assertEqual(s["precision"], round(3 / 5, 3))
        self.assertEqual(s["missed"], [])
        self.assertEqual(s["high_conf_extra"], ["IL6"])            # low-conf IL10 not flagged

    def test_missed_edges(self):
        s = B.score_topology([{"cytokine": "IL1b"}], {"IL1b", "IL17"})
        self.assertEqual(s["recall"], 0.5)
        self.assertEqual(s["missed"], ["IL17"])


class TestScoreForm(unittest.TestCase):
    def test_order_and_combination_match(self):
        self.assertEqual(B.score_form({"proliferation_order": "zeroth", "combination": "capped_sum"},
                                      "zeroth", "capped_sum"),
                         {"order": "zeroth", "combination": "capped_sum",
                          "order_match": True, "comb_match": True})
        f = B.score_form({"proliferation_order": "zeroth", "combination": "product"},
                         "zeroth", "capped_sum")
        self.assertTrue(f["order_match"])
        self.assertFalse(f["comb_match"])


if __name__ == "__main__":
    unittest.main()
