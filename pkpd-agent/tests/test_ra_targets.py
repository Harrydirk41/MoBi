"""Steady-state target extraction: the dependency-free xlsx reader and the row parser.

read_sheet is tested on a minimal in-memory .xlsx (shared strings + numbers) so no external
file is needed; extract_from_rows is tested on hand-built rows mimicking MOESM1's layout
(section headers, primary rows, secondary source rows, unit-keyed converted values).
"""

import io
import unittest
import zipfile

from pkpd_agent.engines import xlsx_read as X
from examples.extract_ra_targets import extract_from_rows


def _minimal_xlsx() -> bytes:
    """Build the smallest valid .xlsx with one sheet 'S': a string cell and two numbers."""
    ct = ('<?xml version="1.0"?><Types xmlns="http://schemas.openxmlformats.org/'
          'package/2006/content-types">'
          '<Default Extension="xml" ContentType="application/xml"/>'
          '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-'
          'officedocument.spreadsheetml.sheet.main+xml"/></Types>')
    rels = ('<?xml version="1.0"?><Relationships xmlns="http://schemas.openxmlformats.org/'
            'package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.'
            'openxmlformats.org/officeDocument/2006/relationships/officeDocument" '
            'Target="xl/workbook.xml"/></Relationships>')
    wb = ('<?xml version="1.0"?><workbook xmlns:r="http://schemas.openxmlformats.org/'
          'officeDocument/2006/relationships"><sheets><sheet name="S" sheetId="1" '
          'r:id="rId1"/></sheets></workbook>')
    wbrels = ('<?xml version="1.0"?><Relationships xmlns="http://schemas.openxmlformats.org/'
              'package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.'
              'openxmlformats.org/officeDocument/2006/relationships/worksheet" '
              'Target="worksheets/sheet1.xml"/></Relationships>')
    sst = ('<?xml version="1.0"?><sst xmlns="http://schemas.openxmlformats.org/'
           'spreadsheetml/2006/main" count="1" uniqueCount="1"><si><t>FLS</t></si></sst>')
    sheet = ('<?xml version="1.0"?><worksheet xmlns="http://schemas.openxmlformats.org/'
             'spreadsheetml/2006/main"><sheetData>'
             '<row r="1"><c r="A1" t="s"><v>0</v></c><c r="B1"><v>441</v></c>'
             '<c r="C1"><v>2.28E7</v></c></row></sheetData></worksheet>')
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        z.writestr("[Content_Types].xml", ct)
        z.writestr("_rels/.rels", rels)
        z.writestr("xl/workbook.xml", wb)
        z.writestr("xl/_rels/workbook.xml.rels", wbrels)
        z.writestr("xl/sharedStrings.xml", sst)
        z.writestr("xl/worksheets/sheet1.xml", sheet)
    return buf.getvalue()


class TestXlsxRead(unittest.TestCase):
    def setUp(self):
        self.path = "/tmp/_mini_test.xlsx"
        with open(self.path, "wb") as fh:
            fh.write(_minimal_xlsx())

    def test_sheet_names(self):
        self.assertEqual(X.sheet_names(self.path), ["S"])

    def test_reads_string_and_numbers(self):
        rows = X.read_sheet(self.path, "S")
        self.assertEqual(rows[0][0], "FLS")            # shared string resolved
        self.assertEqual(rows[0][1], 441)              # integer
        self.assertAlmostEqual(rows[0][2], 2.28e7)     # scientific notation

    def test_missing_sheet_raises(self):
        with self.assertRaises(KeyError):
            X.read_sheet(self.path, "Nope")


class TestExtractRows(unittest.TestCase):
    ROWS = [
        ["Cell numbers"],
        ["serial no", "Cell type", "Mean/Median", "SD/IQR", "Lower bound",
         "Upper bound", "Units", "Converted to cells/mL*"],
        [1, "Fibrocytes Like Synoviocytes (FLS)", 441, "35-2405", 35, 2405,
         "cells/mm2", "2.28E+07"],
        [28, 479, 334, 145, 1147, "cells/mm2", 35, "secondary source row"],
        [10, "Macrophages", 2190, "NR", 192, 6765, "cells/mm2", "4.94E+07"],
        ["Cytokine numbers"],
        [2, "IL-6", 289, 238, 51, 714, "ng/mL", 289],
        [14, "upper bound", "Units"],                  # a header leaked as a row
    ]

    def setUp(self):
        self.t = extract_from_rows(self.ROWS)

    def test_primary_cell_row_parsed_with_model_unit(self):
        fls = next(x for x in self.t if x["model_species"] == "FLS")
        self.assertEqual(fls["kind"], "cell")
        self.assertAlmostEqual(fls["target_model_unit"], 2.28e7)

    def test_secondary_source_row_skipped(self):
        # the '28 | 479 | ...' row has a numeric second col -> not a target
        self.assertFalse(any(x["name"] == "479" for x in self.t))

    def test_header_leak_skipped(self):
        self.assertFalse(any(x["name"].lower() == "upper bound" for x in self.t))

    def test_cytokine_section_and_alias(self):
        il6 = next(x for x in self.t if x["model_species"] == "IL6")
        self.assertEqual(il6["kind"], "cytokine")
        self.assertEqual(il6["target_units"], "ng/mL")


if __name__ == "__main__":
    unittest.main()
