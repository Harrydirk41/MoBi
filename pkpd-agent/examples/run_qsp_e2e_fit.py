r"""A) Honest parameter origin: look up ONLY what is literature, FIT the rest yourself.

The end-to-end build looks up each regulator's strength in the curated table. But that table mixes
two kinds of number: values with a real citation (literature - fair to look up) and values with no
citation (top-down FITTED by the original modellers - looking those up would be reading the answer).
This driver splits the IL-6 hub parameters by that provenance and then does the honest thing:

  * LOOK UP the literature-cited Max fold-changes and the measured clearance kcl.
  * FIT the non-literature parameters from data the agent could actually obtain:
      - kg (baseline secretion): pinned uniquely by the disease steady-state target -> identifiable.
      - K  (half-effect concentrations): under-determined by a single steady-state target - many
        K-sets reproduce it. We show two that match TRAINING exactly yet DIVERGE on held-out, so
        the honest verdict is: K needs per-cytokine dose-response (benchmark A) to pin. The model
        is sloppy, so the held-out spread is small - unpinned does not mean unusable.

    python -m examples.run_qsp_e2e_fit

Pure - reads the provenance + target JSONs; reuses the assembly/integration from the end-to-end
driver so the fitted model is actually emitted and run, not just computed on paper.
"""

from __future__ import annotations

import json
import os
import tempfile

from pkpd_agent.engines import model_assembly as MA
from pkpd_agent.engines.sbml_import import sbml_to_network
from examples.run_qsp_end_to_end import assemble, integrate

_DATA = os.path.join(os.path.dirname(__file__), "..", "projects", "vantage_ra", "data")


def main() -> None:
    prov = {p["name"]: p for p in json.load(open(os.path.join(_DATA, "param_provenance.json")))}
    tg = {t["model_species"]: t for t in json.load(open(os.path.join(_DATA,
          "steady_state_targets.json"))) if t.get("model_species")}
    levels = {c: float(tg[c]["target_model_unit"]) for c in tg
              if tg[c].get("target_model_unit") is not None}
    target = float(tg["IL6"]["target_model_unit"])
    kcl = float(prov["kcl_IL6"]["value_from_reference"])

    # ---- provenance split ----
    lit_max, need_fit = {}, []
    for pref in ("IL6SecFLS_Maxby", "IL6SecMacro_Maxby"):
        for n, p in prov.items():
            if n.startswith(pref):
                c = n.split("Maxby")[-1]
                if p.get("from_literature") and p.get("value_from_reference") is not None \
                        and c in levels and c not in lit_max:
                    lit_max[c] = (float(p["value_from_reference"]),
                                  p.get("reference") or p.get("citation"))
    regs = sorted(lit_max)
    print("== parameter origin split for the IL-6 hub ==")
    print("  LOOK UP (literature-cited, fair):")
    for c in regs:
        print(f"    Max_{c:5} = {lit_max[c][0]:<8.4g} [{lit_max[c][1]}]")
    print(f"    kcl_IL6   = {kcl:<8.4g} [measured clearance]")
    print("  MUST FIT (no citation - top-down fitted; looking these up = reading the answer):")
    print("    kg_IL6     (baseline secretion)")
    print(f"    K_<cyt>    (half-effect concentration, one per regulator: {regs})")

    truth_maxes = {c: lit_max[c][0] for c in regs}
    chosen = [{"cytokine": c, "direction": "up" if truth_maxes[c] > 1 else "down"} for c in regs]
    motif = {"proliferation_order": "zeroth", "combination": "product", "cap": None}

    # ---- FIT kg: identifiable from the steady-state target ----
    # assemble() sets K = level (one plausible K-choice) and fits kg to the target.
    spec, _, _ = assemble(chosen, truth_maxes, levels, kcl, target, motif=motif)
    kg = next(p["value"] for p in spec["parameters"] if p["name"] == "kg_IL6")
    ax = os.path.join(tempfile.gettempdir(), "e2e_fit.xml")
    open(ax, "w", encoding="utf-8").write(MA.to_sbml(spec))
    ss = integrate(sbml_to_network(ax), {c: levels[c] for c in levels if c != "IL6"})["IL6"]
    print(f"\n== FIT kg to the steady-state target -> {kg:.4g};  emitted model integrates to "
          f"{ss:.4g} (target {target:g}) ==")
    print("  kg is IDENTIFIABLE: one steady-state equation pins one baseline rate uniquely.")

    # ---- FIT K: under-determined by the same single target ----
    # Analytic (assemble ties K to the operating level, so vary K here directly): with K_c a free
    # half-effect, kg is refit for EACH K-set to hit the target, so all match training exactly;
    # only held-out separates them. effect(x) = prod_c (1 + (mx-1) x/(K_c + x)).
    def eff(kfac, x):
        p = 1.0
        for c in regs:
            K = levels[c] * kfac
            p *= 1.0 + (truth_maxes[c] - 1.0) * x[c] / (K + x[c])
        return p

    print("\n== FIT K from the same steady-state target: UNDER-DETERMINED ==")
    print(f"  {'K choice':<22}{'train IL-6':>11}{'held-out (anti-IL-1)':>22}")
    base = {c: levels[c] for c in regs}
    held = dict(base); held["IL1b"] = levels["IL1b"] * 0.1
    preds = {}
    for label, kfac in [("tight  (K = 0.1x level)", 0.1), ("loose  (K = 10x level)", 10.0)]:
        eb = eff(kfac, base)                               # kg refit so train == target exactly
        tr = target * eff(kfac, base) / eb                 # == target by construction
        ho = target * eff(kfac, held) / eb
        preds[label] = ho
        print(f"  {label:<22}{tr:>11.4g}{ho:>22.4g}")
    a, b = list(preds.values())
    spread = abs(a - b) / max(a, b)
    print(f"\n  both K-fits match TRAINING exactly ({target:g}), but held-out (anti-IL-1) differs "
          f"by {spread:.0%}.")
    print("  -> K is NOT pinned by steady state alone, and here the cost is LARGE: a "
          f"{spread:.0%} swing at")
    print("     anti-IL-1, an intervention on the dominant driver. (For weak or off-target")
    print("     operating points the same K ambiguity costs little - it is operating-point")
    print("     dependent, not uniformly small.) Honest recipe & limit: look up the literature Max")
    print("     values, fit kg (the steady state pins it), but the half-effects genuinely need")
    print("     per-cytokine dose-response (benchmark A) to predict strong-driver therapies - that")
    print("     data, not the answer table, is the real missing ingredient.")


if __name__ == "__main__":
    main()
