"""Route-A reading-material assembly: citation parsing, title extraction, abstract enrich.

Pure functions + a STUB fetcher (no network). The live PubMed path (pubmed_fetch) is
exercised only in the runner. The fixture mimics the ragged multi-line numbering real
papers print (number alone on a line; 10+ leading the text).
"""

import unittest

from pkpd_agent.engines import llm_references as R

_REFS = """References
1.
Almutairi, K. A.-O. et al. The Prevalence of Rheumatoid Arthritis: A
Systematic Review. J. Rheumatol. 48, 669-676 (2021).
2.
Smolen, J. S. et al. Rheumatoid arthritis. Nat. Rev. Dis. Prim. 4,
18001 (2018).
10. Ray, T. et al. QSP Model of Rheumatoid Arthritis. PAGE 2019.
"""


class TestParse(unittest.TestCase):
    def setUp(self):
        self.c = R.parse_references(_REFS)

    def test_count_and_numbers(self):
        # only the contiguous 1,2 survive; the out-of-sequence '10' is not #3, so dropped
        self.assertEqual([x["n"] for x in self.c], [1, 2])

    def test_lines_collapsed(self):
        self.assertIn("Systematic Review", self.c[0]["text"])
        self.assertNotIn("\n", self.c[0]["text"])

    def test_stray_year_marker_not_a_citation(self):
        # '18001 (2018)' contains no leading number marker; nothing spurious added
        self.assertTrue(all(x["text"] for x in self.c))


class TestTitle(unittest.TestCase):
    def test_et_al_title(self):
        t = R._title_of("Smolen, J. S. et al. Rheumatoid arthritis. Nat. Rev. Dis. 4, 1 (2018).")
        self.assertEqual(t, "Rheumatoid arthritis")

    def test_multi_author_no_etal(self):
        t = R._title_of("Frey, N., Grange, S. & Woodworth, T. Population PK of tocilizumab. "
                        "J Clin Pharmacol. 50, 754 (2010).")
        self.assertEqual(t, "Population PK of tocilizumab")


class TestFetchAbstracts(unittest.TestCase):
    def test_pluggable_fetch_fills_abstracts(self):
        calls = []
        def fake(title):
            calls.append(title)
            return f"ABSTRACT for {title}"
        c = R.parse_references(_REFS)
        R.fetch_abstracts(c, fetch=fake, sleep=0)
        self.assertEqual(len(calls), 2)
        self.assertTrue(all(x["abstract"].startswith("ABSTRACT") for x in c))

    def test_format_includes_abstract(self):
        c = [{"n": 1, "text": "Foo et al. Bar.", "abstract": "an abstract"}]
        out = R.format_references(c)
        self.assertIn("1. Foo et al. Bar.", out)
        self.assertIn("ABSTRACT: an abstract", out)

    def test_fetch_miss_is_empty_not_error(self):
        c = R.parse_references(_REFS)
        R.fetch_abstracts(c, fetch=lambda t: "", sleep=0)
        self.assertEqual([x["abstract"] for x in c], ["", ""])


if __name__ == "__main__":
    unittest.main()
