r"""Extract each OSP model's PUBLISHED goodness-of-fit GMFE from its evaluation report.

Every OSP evaluation report states the model's own GMFE in a small markdown table
(a "Group | GMFE" table), split into model-building, model-validation/verification,
and an "All" row. That published number is the TRUE yardstick: a faithful harness
re-running the reference model must reproduce roughly this GMFE. If it does not, the
harness is not modelling the reference correctly (food/formulation/dissolution
handling, observed->simulation mapping) and every agent-vs-reference comparison is
built on sand.

This scans every ``*_evaluation_report.md`` and writes a machine-readable
``published_gmfe.json`` mapping model -> {all, building, validation, groups{...}}.

    python -m examples.published_gmfe            # print the table, write the json
"""

from __future__ import annotations

import glob
import json
import os
import re

_HERE = os.path.dirname(os.path.abspath(__file__))
_LIB = os.path.abspath(os.path.join(_HERE, "..", "..", "OSP-PBPK-Model-Library"))


def _tables_with_gmfe(text: str) -> list[list[tuple[str, float]]]:
    """Return every markdown table that has a GMFE column, as a list of
    (row-label, gmfe) pairs. Robust to column order and extra columns."""
    out = []
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        if "|" in line and "gmfe" in line.lower():
            header = [c.strip().lower() for c in line.strip().strip("|").split("|")]
            try:
                gcol = next(k for k, c in enumerate(header) if c == "gmfe")
            except StopIteration:
                i += 1
                continue
            # label column = first non-gmfe column
            lcol = next((k for k in range(len(header)) if k != gcol), 0)
            rows = []
            j = i + 1
            if j < len(lines) and set(lines[j].replace("|", "").strip()) <= set("-: "):
                j += 1                                   # skip the |---|---| separator
            while j < len(lines) and "|" in lines[j]:
                cells = [c.strip() for c in lines[j].strip().strip("|").split("|")]
                if len(cells) > max(gcol, lcol):
                    m = re.search(r"[-+]?\d+\.?\d*", cells[gcol])
                    if m:
                        rows.append((cells[lcol], float(m.group())))
                j += 1
            if rows:
                out.append(rows)
            i = j
        else:
            i += 1
    return out


def _classify(rows: list[tuple[str, float]]) -> dict:
    """Collapse a GMFE table's rows into all / building / validation + raw groups."""
    rec: dict = {"groups": {lbl: v for lbl, v in rows}}
    for lbl, v in rows:
        l = lbl.lower()
        if l.strip() in ("all", "overall") or l.strip().startswith("all "):
            rec["all"] = v
        if "build" in l:
            rec.setdefault("building", v)
        if "valid" in l or "verif" in l or "evaluat" in l or "test" in l:
            rec.setdefault("validation", v)
    # no explicit "All" row: use the max group as a conservative overall estimate
    if "all" not in rec and rows:
        rec["all"] = max(v for _, v in rows)
        rec["all_is_estimated"] = True
    return rec


def parse_report(path: str) -> dict:
    with open(path, encoding="utf-8") as fh:
        text = fh.read()
    tables = _tables_with_gmfe(text)
    if not tables:
        return {"gmfe_found": False}
    # merge all GMFE tables (a metabolite-cascade report has one per molecule);
    # keep the primary compound's "All" as the headline, but retain every group.
    merged: list[tuple[str, float]] = []
    for t in tables:
        merged += t
    rec = _classify(merged)
    rec["gmfe_found"] = True
    rec["n_gmfe_tables"] = len(tables)
    return rec


def main() -> None:
    reports = sorted(glob.glob(os.path.join(_LIB, "*", "*_evaluation_report.md")))
    result = {}
    print(f"{'model':22} {'all':>6} {'build':>6} {'valid':>6}  n_tables")
    for r in reports:
        model = os.path.basename(os.path.dirname(r))
        rec = parse_report(r)
        result[model] = rec
        if rec.get("gmfe_found"):
            print(f"{model:22} {rec.get('all','-'):>6} {rec.get('building','-'):>6} "
                  f"{rec.get('validation','-'):>6}  {rec.get('n_gmfe_tables')}")
        else:
            print(f"{model:22} {'(no GMFE table found)':>6}")
    out = os.path.join(_LIB, "published_gmfe.json")
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(result, fh, indent=1, ensure_ascii=False)
    n = sum(1 for v in result.values() if v.get("gmfe_found"))
    print(f"\nparsed a published GMFE for {n}/{len(result)} models -> {out}")


if __name__ == "__main__":
    main()
