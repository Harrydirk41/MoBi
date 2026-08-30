r"""The fair-fit question, in miniature: an under-determined fit can MATCH the training data yet
DIVERGE on held-out - so 'the fit succeeded on training' does not mean the model generalizes.

Using the real data (MOESM2 Max fold-changes + MOESM1 steady-state target and cytokine levels),
we fit the FLS half-effect K's to the disease steady-state FLS target. Because 5 K's face one
aggregate target the fit is under-determined: many K-sets reproduce it (the baseline rate kg
re-absorbs the difference). We take TWO such solutions - both matching training exactly - and
predict FLS at a HELD-OUT operating point (an anti-IL6 therapy that drops IL6 sharply). If the
two training-equivalent fits predict DIFFERENT held-out FLS, the under-determined fit does not
generalize - the honest answer to 'the paper fit it, why can't the agent': the agent CAN fit
(match training), but its particular solution is not pinned, so held-out is not guaranteed.

    python -m examples.run_qsp_fair_fit

Pure - reads only the two JSONs; no MATLAB, no answer-model K values (the K's are what we fit).
"""

from __future__ import annotations

import json
import os

_DATA = os.path.join(os.path.dirname(__file__), "..", "projects", "vantage_ra", "data")


def effect(maxes, levels, ks, cap=10.0):
    s = sum((maxes[c] - 1) * levels[c] / (ks[c] + levels[c]) for c in maxes)
    return 1.0 + min(cap, s)


def main() -> None:
    prov = {p["name"]: p for p in json.load(open(os.path.join(_DATA, "param_provenance.json")))}
    tg = {t["model_species"]: t for t in json.load(open(os.path.join(_DATA,
          "steady_state_targets.json"))) if t.get("model_species")}
    cyts = ["TNFa", "IL1b", "IL17", "TGFb", "IL6"]
    maxes = {c: float(prov[f"FLSProlif_Maxby{c}"]["value_from_reference"]) for c in cyts}
    kd = float(prov["kd_FLS_Baseline"]["value_from_reference"])
    levels = {c: float(tg[c]["target_model_unit"]) for c in cyts if c in tg}
    maxes = {c: maxes[c] for c in maxes if c in levels}
    fls_target = float(tg["FLS"]["target_model_unit"])

    print("== two DIFFERENT K-fits, both matched to the same real training target ==")
    print(f"  training target: FLS steady state = {fls_target:g}")
    # two solutions: 'tight' K (small, cytokines near saturation) and 'loose' K (large)
    solutions = {}
    for label, kfac in [("solution A (small K)", 0.1), ("solution B (large K)", 10.0)]:
        ks = {c: levels[c] * kfac for c in maxes}
        eff = effect(maxes, levels, ks)
        kg = fls_target * kd / eff                     # re-fit kg so BOTH match training exactly
        solutions[label] = (ks, kg)
        print(f"  {label:22} effect={eff:.3f}  kg={kg:.3g}  -> training FLS={kg/kd*eff:.3g}")

    # held-out operating point: an anti-IL6 therapy drops IL6 to 10% (a different milieu, not fit)
    held = dict(levels); held["IL6"] = levels["IL6"] * 0.1
    print(f"\n== held-out operating point (anti-IL6 therapy: IL6 {levels['IL6']:g} -> "
          f"{held['IL6']:g}) ==")
    preds = {}
    for label, (ks, kg) in solutions.items():
        fls = kg / kd * effect(maxes, held, ks)
        preds[label] = fls
        print(f"  {label:22} predicts held-out FLS = {fls:.3g}")

    a, b = list(preds.values())
    disagree = abs(a - b) / max(abs(a), abs(b))
    print(f"\n== result ==")
    print(f"  both fits match TRAINING exactly ({fls_target:g}), but on HELD-OUT they differ by "
          f"{disagree:.0%}")
    print("\n  -> the fit 'succeeding' on training does NOT pin the model: two solutions that are "
          "identical\n     on the data diverge where it matters. This is why 'the paper fit it' "
          "is not the whole story\n     - the agent CAN fit (match training) with the same data, "
          "but its particular under-determined\n     solution is not guaranteed to generalize. "
          "Pinning it needs the per-cytokine dose-responses\n     (benchmark A) or the expert "
          "judgment + held-out validation the paper's team supplied.")


if __name__ == "__main__":
    main()
