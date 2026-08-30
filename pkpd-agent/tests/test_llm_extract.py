"""Data-acquisition benchmark: prompt building, extraction parsing, grading. Pure - a stub
call stands in for the web-enabled LLM."""

import unittest

from pkpd_agent.engines import llm_extract as EX
from examples.extract_ra_provenance import extract_from_rows


class TestExtractValue(unittest.TestCase):
    PROV = {"name": "FLSProlif_MaxbyIL6", "units": "dimensionless",
            "reference": "Hashizume et al, 2018", "figure": "Fig 2B",
            "experiment": "FLS proliferation under IL-6", "description": "fold change"}

    def test_prompt_mentions_reference_and_figure(self):
        p = EX.build_prompt(self.PROV)
        self.assertIn("Hashizume et al, 2018", p)
        self.assertIn("Fig 2B", p)

    def test_extract_parses(self):
        call = lambda s, u: ('{"found_paper": true, "value": 2.0, "in_figure_only": true, '
                             '"note": "read from Fig 2B"}')
        r = EX.extract_value(self.PROV, call)
        self.assertTrue(r["found_paper"])
        self.assertEqual(r["value"], 2.0)
        self.assertTrue(r["in_figure_only"])

    def test_extract_handles_null_value(self):
        call = lambda s, u: '{"found_paper": false, "value": null, "note": "paywalled"}'
        r = EX.extract_value(self.PROV, call)
        self.assertIsNone(r["value"])


class TestGrade(unittest.TestCase):
    def test_hit_within_tol(self):
        self.assertTrue(EX.grade(2.1, 2.0, tol=0.25)["hit"])

    def test_miss_outside_tol(self):
        g = EX.grade(1.0, 2.0, tol=0.25)
        self.assertFalse(g["hit"])
        self.assertTrue(g["extracted"])
        self.assertAlmostEqual(g["rel_err"], 0.5)

    def test_no_value_not_extracted(self):
        g = EX.grade(None, 2.0)
        self.assertFalse(g["extracted"])
        self.assertFalse(g["hit"])


class TestProvenanceExtract(unittest.TestCase):
    def test_flags_override_and_literature(self):
        rows = [
            ["FLS parameters"],
            ["name", "Units", "Value in the model", "Value from reference", "Description",
             "H/D", "Reference", "Fig", "Exp", "Calc", "Comments"],
            ["FLSProlif_MaxbyIL6", "dimensionless", 2.0, 2.0, "d", "RA", "Hashizume 2018",
             "Fig 2B", "prolif", "ratio", ""],
            ["IL6SecFLS_MaxbyIL1b", "ng", 1.0, 69.0, "d", "RA", "Georganas 2000",
             "Fig 1A", "secretion", "", ""],
        ]
        out = extract_from_rows(rows)
        by = {p["name"]: p for p in out}
        self.assertTrue(by["FLSProlif_MaxbyIL6"]["from_literature"])
        self.assertFalse(by["FLSProlif_MaxbyIL6"]["overridden"])
        self.assertTrue(by["IL6SecFLS_MaxbyIL1b"]["overridden"])   # 1 vs 69


if __name__ == "__main__":
    unittest.main()
