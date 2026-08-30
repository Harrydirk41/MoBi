r"""Extract the per-parameter PROVENANCE from Supplementary Data 2 (MOESM2) into JSON - the
answer key for the data-acquisition benchmark (A): for each parameter, which reference /
figure it came from, the value read from that reference, the value used in the model, and the
experiment description. This lets us grade whether an LLM can READ the cited paper and EXTRACT
the value the modellers used.

    python -m examples.extract_ra_provenance \
        --xlsx ..\41540_2024_454_MOESM2_ESM.xlsx \
        --out projects\vantage_ra\data\param_provenance.json
"""

from __future__ import annotations

import argparse
import json
import os
import re


def _num(x):
    try:
        return float(str(x).strip())
    except (TypeError, ValueError):
        return None


def extract_from_rows(rows: list) -> list[dict]:
    out = []
    for r in rows:
        if not r or not isinstance(r[0], str):
            continue
        name = r[0].strip()
        if not name or name.lower().endswith("parameters") or name.lower() == "units":
            continue                                   # section header / blank
        def col(i):
            return (r[i] if len(r) > i and str(r[i]).strip() not in ("", "NR") else None)
        vm, vr = _num(col(2)), _num(col(3))
        if vm is None and vr is None:                  # a continuation/source row w/o a param
            continue
        out.append({
            "name": name,
            "units": col(1),
            "value_in_model": vm,
            "value_from_reference": vr,
            "reference": col(6),
            "figure": col(7),
            "experiment": col(8),
            "calculation": col(9),
            "from_literature": vr is not None,
            "overridden": (vm is not None and vr is not None and vm != 0
                           and not (0.5 <= vr / vm <= 2.0)),
        })
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xlsx", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    from pkpd_agent.engines import xlsx_read as X
    prov = extract_from_rows(X.read_sheet(args.xlsx, "Model_parameters"))
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as fh:
        json.dump(prov, fh, indent=2)
    lit = [p for p in prov if p["from_literature"]]
    withref = [p for p in lit if p["reference"]]
    print(f"{len(prov)} parameters; {len(lit)} have a literature value; "
          f"{len(withref)} name a reference; "
          f"{sum(1 for p in lit if p['overridden'])} were overridden by the model")
    print(f"-> {args.out}")


if __name__ == "__main__":
    main()
