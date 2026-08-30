r"""L1 calibration: calibrate the FLS subsystem using ONLY the real curated data (MOESM1
steady-state targets + MOESM2 literature parameter values) - NO answer-model oracle, no
synthesized isolating experiments. Tests the un-scaffolded question: with the data a modeller
actually has, which parameters can be pinned and which cannot?

Finding this exposes: MOESM2 supplies the bottom-up parameters (the Max fold-changes, the
apoptosis rate), but NOT the half-effect K's, Hill slopes, or the baseline proliferation rate -
the ones the paper fitted top-down. So from the real data the baseline rate is identifiable
from the steady-state target (a 1-D fit), but the shape parameters (K's) are NOT: any K's
reproduce the one target once the baseline rate re-absorbs the difference. The data that WOULD
identify them - per-cytokine dose-responses - lives in the primary-paper figures (benchmark A,
the hard access problem), not the supplement.

    python -m examples.run_qsp_l1_calibrate

Pure - reads only projects/vantage_ra/data/{steady_state_targets,param_provenance}.json.
No MATLAB, no answer model.
"""

from __future__ import annotations

import json
import os

_DATA = os.path.join(os.path.dirname(__file__), "..", "projects", "vantage_ra", "data")


def _load(name):
    with open(os.path.join(_DATA, name), encoding="utf-8") as fh:
        return json.load(fh)


def effect(maxes, levels, ks, cap=10.0):
    """capped-sum motif effect: 1 + min(cap, sum_i (Max_i-1) * X_i/(K_i+X_i))."""
    s = sum((maxes[c] - 1) * levels[c] / (ks[c] + levels[c]) for c in maxes)
    return 1.0 + min(cap, s)


def main() -> None:
    prov = {p["name"]: p for p in _load("param_provenance.json")}
    targets = {t["model_species"]: t for t in _load("steady_state_targets.json")
               if t.get("model_species") and t.get("target_model_unit") is not None}

    cyts = ["TNFa", "IL1b", "IL17", "TGFb", "IL6"]
    maxes, have, missing = {}, [], []
    for c in cyts:
        p = prov.get(f"FLSProlif_Maxby{c}")
        if p and p.get("from_literature"):
            maxes[c] = float(p["value_from_reference"])
            have.append(f"FLSProlif_Maxby{c}={maxes[c]:g} ({p['reference']})")
    kd = prov.get("kd_FLS_Baseline")
    kd_v = float(kd["value_from_reference"]) if kd and kd.get("from_literature") else None

    print("== what the REAL curated data (MOESM1/MOESM2) gives, with NO answer model ==")
    print("  bottom-up parameters WITH a literature value (MOESM2):")
    for h in have:
        print(f"    {h}")
    print(f"    kd_FLS_Baseline={kd_v} ({kd['reference']})" if kd_v else "    kd: MISSING")
    for name in ["kg_FLS_Baseline"] + [f"HalfEffectConc_FLSProlif_by{c}" for c in cyts] \
            + [f"Slope_FLSProlif_by{c}" for c in cyts]:
        p = prov.get(name)
        if not p or not p.get("from_literature"):
            missing.append(name)
    print(f"\n  parameters NOT in the real data (top-down fitted, must be inferred): "
          f"{len(missing)}")
    for m in missing:
        print(f"    {m}")

    # cytokine levels (the clamp) + FLS target - from MOESM1
    levels = {c: float(targets[c]["target_model_unit"]) for c in cyts if c in targets}
    fls_target = targets.get("FLS", {}).get("target_model_unit")
    print(f"\n  steady-state targets from MOESM1: FLS={fls_target:g}, "
          f"cytokines={ {c: round(levels[c], 3) for c in levels} }")
    maxes = {c: maxes[c] for c in maxes if c in levels}

    if not (fls_target and kd_v and maxes):
        print("insufficient real data resolved."); return

    # the identifiability demonstration: baseline rate kg is fit to the FLS target; the K's
    # cannot be - try three very different K-sets, refit kg each, all reproduce the SAME target.
    print("\n== can the real data pin the top-down parameters? ==")
    print("  fitting the baseline rate kg to the FLS target, under three different K guesses:")
    print(f"  {'K assumption':22} {'effect':>8} {'fitted kg':>14} {'reproduces FLS?':>16}")
    for label, kfac in [("K = 0.1 x level", 0.1), ("K = 1 x level", 1.0), ("K = 10 x level", 10.0)]:
        ks = {c: max(levels[c] * kfac, 1e-12) for c in maxes}
        eff = effect(maxes, levels, ks)
        kg = fls_target * kd_v / eff                   # 1-D fit: kg absorbs whatever eff is
        fls = kg / kd_v * eff
        print(f"  {label:22} {eff:8.3f} {kg:14.3g} {('yes (%.2g)' % fls):>16}")

    print("\n== conclusion ==")
    print("  The FLS steady-state target pins ONLY the baseline rate kg (a 1-D fit): every K "
          "guess\n  reproduces it, because kg re-absorbs the difference. The shape parameters "
          "(the K's,\n  slopes) are UNIDENTIFIED from the real curated data - and the data that "
          "would identify\n  them (per-cytokine dose-responses) is in the primary-paper figures "
          "(benchmark A), not the\n  supplement. Removing the answer-model oracle exposes exactly "
          "this: the oracle was doing\n  the job of the per-cytokine experiments the real data "
          "does not contain.")


if __name__ == "__main__":
    main()
