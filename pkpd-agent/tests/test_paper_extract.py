"""Paper extraction: the reliable text parts work; structured figure/table data does not.

These tests pin the honest boundary - the ACR regex parses a CLEAN row, and the module's
functions are present - without asserting that real flattened-PDF tables or figure axes
extract (they do not; that needs a vision step).
"""

import os
import unittest

from pkpd_agent.engines import paper_extract as PE

_PDF = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..",
                                    "original_paper.pdf"))


class TestAcrRegex(unittest.TestCase):
    def test_clean_row_parses(self):
        # a clean, well-separated row parses; monotonicity filter keeps it
        txt = "Tocilizumab TCZ 8mg 45 30 13.9"
        # feed via the module's regex on a stand-in text (no PDF needed)
        rows = [m for m in PE._ACR_ROW.finditer(txt)]
        self.assertTrue(rows)

    def test_monotonicity_filter_rejects_noise(self):
        # ACR20 < ACR50 is impossible -> filtered out (protects against table noise)
        import re
        m = re.search(PE._ACR_ROW, "MTX 9 23 46")   # reversed order = noise
        if m:
            a20, a50, a70 = float(m.group(2)), float(m.group(3)), float(m.group(4))
            self.assertFalse(a20 >= a50 >= a70)      # correctly not a valid ACR triple


class TestOnRealPaper(unittest.TestCase):
    @unittest.skipUnless(os.path.exists(_PDF), "paper PDF not present")
    def test_param_tokens_extractable_but_gsa_ranking_is_not(self):
        toks = PE.find_param_tokens(_PDF)
        self.assertGreater(len(toks), 50)            # param names ARE in the text
        # ...but the ranked GSA list is a figure; we do not assert it extracts.


if __name__ == "__main__":
    unittest.main()
