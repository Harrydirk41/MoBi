"""Tiny dependency-free .xlsx reader (stdlib only: zipfile + regex over the OOXML parts).

The project must run without pandas/openpyxl (the pipeline already avoids openpyxl for Vpop
I/O). This reads a worksheet into rows of strings/numbers - enough to pull the paper's
supplementary tables (steady-state targets, parameters) into JSON. Not a general xlsx
library: it handles shared strings, inline numbers, and blank cells; it ignores styles,
formulas' cached vs live values (returns the cached <v>), merged-cell geometry, and dates.
"""

from __future__ import annotations

import re
import zipfile


def sheet_names(path: str) -> list[str]:
    with zipfile.ZipFile(path) as z:
        wb = z.read("xl/workbook.xml").decode("utf-8", "replace")
    return re.findall(r'<sheet[^>]*name="([^"]+)"', wb)


def read_sheet(path: str, sheet: str) -> list[list]:
    """Return the worksheet ``sheet`` as a list of rows (each a list of cell values, strings
    or floats; blank cells are ''). Rows are trimmed of trailing blanks."""
    with zipfile.ZipFile(path) as z:
        wb = z.read("xl/workbook.xml").decode("utf-8", "replace")
        pairs = re.findall(r'<sheet[^>]*name="([^"]+)"[^>]*r:id="(rId\d+)"', wb)
        rels = z.read("xl/_rels/workbook.xml.rels").decode("utf-8", "replace")
        idmap = dict(re.findall(r'Id="(rId\d+)"[^>]*Target="([^"]+)"', rels))
        by_name = {nm: idmap[rid] for nm, rid in pairs if rid in idmap}
        if sheet not in by_name:
            raise KeyError(f"sheet {sheet!r} not found; have {sorted(by_name)}")
        try:
            ss_xml = z.read("xl/sharedStrings.xml").decode("utf-8", "replace")
            sst = [re.sub(r"<[^>]+>", "", s)
                   for s in re.findall(r"<si>(.*?)</si>", ss_xml, re.S)]
        except KeyError:
            sst = []
        data = z.read("xl/" + by_name[sheet].lstrip("/")).decode("utf-8", "replace")

    rows: list[list] = []
    for row in re.findall(r"<row[^>]*>(.*?)</row>", data, re.S):
        cells: list = []
        for cell in re.findall(r"<c\b(.*?)</c>", row, re.S):
            head = cell.split(">", 1)[0]
            m = re.search(r"<v>(.*?)</v>", cell, re.S)
            if not m:
                cells.append("")
                continue
            v = m.group(1)
            if 't="s"' in head and v.isdigit():          # shared-string index
                cells.append(sst[int(v)] if int(v) < len(sst) else "")
            elif 't="str"' in head or 't="inlineStr"' in head:
                cells.append(v)
            else:                                          # numeric
                try:
                    f = float(v)
                    cells.append(int(f) if f.is_integer() else f)
                except ValueError:
                    cells.append(v)
        while cells and cells[-1] == "":
            cells.pop()
        rows.append(cells)
    return rows
