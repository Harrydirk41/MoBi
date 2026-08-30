r"""Extract the paper's steady-state calibration TARGETS from Supplementary Data 1 (MOESM1)
into a machine-readable JSON - the top-down cell densities and cytokine concentrations the
model was calibrated to reproduce. This is the fitting TARGET data (the answer key for a
Stage-0 calibration), the counterpart to ra_params_esm2.json (the fitted parameter values).

    python -m examples.extract_ra_targets \
        --xlsx ..\41540_2024_454_MOESM1_ESM.xlsx \
        --out projects\vantage_ra\data\steady_state_targets.json

Robust to the sheet's shifting sub-layouts by keying each row off its UNITS cell: the value
immediately after the units is the paper's own conversion to the model's unit (cells/mL for
cells, ng/mL for cytokines). Names are mapped to the model's species where an alias is known.
"""

from __future__ import annotations

import argparse
import json
import os
import re

from pkpd_agent.engines import xlsx_read as X

# MOESM1 name -> model species (only the ones the model represents as a node)
_ALIAS = {
    "fibrocytes like synoviocytes (fls)": "FLS", "fls": "FLS",
    "b cells": "BCells", "plasma cells": "PlasmaCells", "macrophages": "Macrophages",
    "th1": "Th1", "th17": "Th17", "treg": "Treg", "cd8 cells": "CTL",
    "endothelial cells( sum of mature and immature vessels)": "Endothelial",
    "il-10": "IL10", "tgf-beta@ time 0": "TGFb", "tnf-alpha": "TNFa", "il-6": "IL6",
    "il-1 beta": "IL1b", "il-17": "IL17", "ifn-gamma": "IFNg", "il-12": "IL12",
    "vegf": "VEGF", "rf-igm": "AutoAb", "gm-csf": "GMCSF", "rantes": "RANTES",
    "baff": "BAFF", "mcp-1": "MCP1", "il23": "IL23", "mip3a": "MIP3",
}
_UNIT = re.compile(r"(cells?/m|vessels/m|%\s*of|pg/m|ng/m|pM|U/m|protein)", re.I)


def _num(x):
    if isinstance(x, (int, float)):
        return float(x)
    s = str(x).strip().replace(",", "")
    m = re.match(r"^[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?$", s)
    return float(s) if m else None


def _units_idx(row):
    for i, c in enumerate(row):
        if isinstance(c, str) and _UNIT.search(c):
            return i
    return None


def extract(xlsx: str) -> list[dict]:
    return extract_from_rows(X.read_sheet(xlsx, "Steady_state_cells_cytokines"))


def extract_from_rows(rows: list) -> list[dict]:
    out, kind = [], None
    for row in rows:
        first = str(row[0]).strip().lower() if row else ""
        if first.startswith("cell numbers"):
            kind = "cell"; continue
        if first.startswith("cytokine numbers"):
            kind = "cytokine"; continue
        if kind is None or len(row) < 3:
            continue
        # a primary target row: col0 a small serial int, col1 a non-numeric NAME
        serial = _num(row[0])
        name = row[1] if len(row) > 1 else ""
        if serial is None or serial > 30 or _num(name) is not None or not str(name).strip():
            continue                                   # secondary source row / header / blank
        if str(name).strip().lower() in {"upper bound", "lower bound", "mean", "units",
                                          "cell type", "mean/median"}:
            continue                                   # a sub-section header leaked as a row
        ui = _units_idx(row)
        if ui is None or ui + 1 >= len(row):
            continue
        units = str(row[ui]).strip()
        converted = _num(row[ui + 1])
        mean = _num(row[2])
        key = str(name).strip().lower()
        out.append({
            "name": str(name).strip(),
            "model_species": _ALIAS.get(key),
            "kind": kind,
            "mean_reported": mean,
            "units_reported": units,
            "target_model_unit": converted,
            "target_units": "cells/mL" if kind == "cell" else "ng/mL",
        })
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xlsx", required=True, help="path to 41540_2024_454_MOESM1_ESM.xlsx")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    targets = extract(args.xlsx)
    mapped = [t for t in targets if t["model_species"] and t["target_model_unit"] is not None]
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as fh:
        json.dump(targets, fh, indent=2)
    print(f"extracted {len(targets)} targets ({len(mapped)} mapped to a model species "
          "with a model-unit value) -> " + args.out)
    for t in targets:
        ms = t["model_species"] or "-"
        v = t["target_model_unit"]
        print(f"  {t['kind']:8} {t['name'][:38]:38} -> {ms:12} "
              f"{v if v is not None else '?'} {t['target_units']}")


if __name__ == "__main__":
    main()
