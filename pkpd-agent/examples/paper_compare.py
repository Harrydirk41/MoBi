r"""Diff my pipeline's model output against the paper's published model predictions.

Reads the three per-patient CSVs written by ``sb_paper_compare.m`` (full Vpop, no
subsampling), rebuilds the paper's Fig 4 cascade with its EXACT inadequate-responder
rules, and prints a three-column table for each arm:

    my model (raw sim, full pop)  |  paper's model (Fig 5/6 bars)  |  real clinical data

The paper's model bars were digitized by eye from Fig 5 / Fig 6 (+/-2-3pp) and the
clinical numbers are the raw Table 1 / ESM1 values. The apples-to-apples comparison
that answers "do I reproduce their model" is column 1 vs column 2 (both are RAW model
outputs). Fig 5's red data bars are placebo-corrected; column 3 here is RAW data, so a
raw-data gap is expected and is not the model-vs-model question.

    python -m examples.paper_compare --dir path\to\csvs

No MATLAB required - this is pure post-processing of the CSVs.
"""

from __future__ import annotations

import argparse
import csv
import os


# --- paper's own numbers (fixed references) ------------------------------- #
# Model predictions: digitized from the BLACK (Simulation) bars, +/-2-3pp.
PAPER_MODEL = {
    "Fig5A MTX naive (Wk12)":        {"ACR20": 27, "ACR50": 16, "ACR70": 6},
    "Fig5C TCZ on MTX-IR (Wk24)":    {"ACR20": 30, "ACR50": 22, "ACR70": 12},
    "Fig6 TCZ on MTX-IR&ADA-IR":     {"ACR20": 34, "ACR50": 24, "ACR70": 14},
}
# Raw clinical data (Table 1 / ESM1), NOT placebo-corrected.
REAL_RAW = {
    "Fig5A MTX naive (Wk12)":        {"ACR20": 46, "ACR50": 23, "ACR70": 9},    # Strand 1999
    "Fig5C TCZ on MTX-IR (Wk24)":    {"ACR20": 45, "ACR50": 30.1, "ACR70": 13.9},  # ROSE
    "Fig6 TCZ on MTX-IR&ADA-IR":     {"ACR20": 50, "ACR50": 28.8, "ACR70": 12.4},  # Emery/RADIATE
}
PAPER_N = {"Fig5A MTX naive (Wk12)": 300, "Fig5C TCZ on MTX-IR (Wk24)": 251,
           "Fig6 TCZ on MTX-IR&ADA-IR": 216}


def _read(path: str) -> dict[int, dict]:
    """patient-index -> row of floats (blank -> None)."""
    out = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            rec = {}
            for k, v in row.items():
                try:
                    rec[k] = float(v)
                except (TypeError, ValueError):
                    rec[k] = None
            out[int(rec["patient"])] = rec
    return out


def _rate(rows: list[dict], key: str) -> float | None:
    vals = [r[key] for r in rows if r.get(key) is not None]
    if not vals:
        return None
    return 100.0 * sum(1 for v in vals if v >= 0.5) / len(vals)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True, help="folder with arm_mtx/ada/tcz.csv")
    args = ap.parse_args()

    mtx = _read(os.path.join(args.dir, "arm_mtx.csv"))
    ada = _read(os.path.join(args.dir, "arm_ada.csv"))
    tcz = _read(os.path.join(args.dir, "arm_tcz.csv"))
    pats = sorted(set(mtx) & set(ada) & set(tcz))
    print(f"patients common to all three arms: {len(pats)}\n")

    # -- the paper's exact IR rules (Fig 4) -------------------------------- #
    #   MTX-IR : ACR<50 at Wk12  (ACR50 flag == 0)
    #   ADA-IR : DAS28-CRP>3.2 AND ACR<50 at Wk24
    mtx_ir = {p for p in pats if (mtx[p].get("ACR50") or 0) < 0.5}
    ada_ir = {p for p in pats
              if (ada[p].get("ACR50") or 0) < 0.5
              and (ada[p].get("DAS28_read") or 0) > 3.2}
    valid = mtx_ir & ada_ir
    print(f"MTX-IR (ACR<50 @Wk12)         n={len(mtx_ir):3d}  (paper 251)")
    print(f"ADA-IR (DAS>3.2 & ACR<50)     n={len(ada_ir):3d}  (paper 217)")
    print(f"validation MTX-IR & ADA-IR    n={len(valid):3d}  (paper 216)\n")

    arms = {
        "Fig5A MTX naive (Wk12)":     [mtx[p] for p in pats],
        "Fig5C TCZ on MTX-IR (Wk24)": [tcz[p] for p in mtx_ir],
        "Fig6 TCZ on MTX-IR&ADA-IR":  [tcz[p] for p in valid],
    }

    hdr = f"{'arm / ACR':32s} {'mine':>7s} {'paper':>7s} {'d(m-p)':>7s} {'real':>7s}"
    for arm, rows in arms.items():
        print("=" * len(hdr))
        print(f"{arm}   (mine n={len(rows)}, paper n={PAPER_N[arm]})")
        print(hdr)
        for acr in ("ACR20", "ACR50", "ACR70"):
            mine = _rate(rows, acr)
            paper = PAPER_MODEL[arm][acr]
            real = REAL_RAW[arm][acr]
            dm = f"{mine - paper:+.1f}" if mine is not None else "  n/a"
            ms = f"{mine:.1f}" if mine is not None else " n/a"
            print(f"  {acr:30s} {ms:>7s} {paper:>7.0f} {dm:>7s} {real:>7.1f}")
    print("=" * len(hdr))
    print("\ncol 'mine'  = my pipeline, raw model output, full Vpop")
    print("col 'paper' = paper's model bars (Fig 5/6, digitized by eye +/-2-3pp)")
    print("col 'real'  = raw clinical data (Table 1 / ESM1), NOT placebo-corrected")
    print("The model-vs-model question is  mine vs paper  (both raw sim). A large")
    print("|d(m-p)| means my pipeline does NOT reproduce their model on that endpoint.")


if __name__ == "__main__":
    main()
