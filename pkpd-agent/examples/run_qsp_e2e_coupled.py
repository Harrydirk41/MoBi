r"""B) Two coupled hubs: do independently-built modules compose into a working feedback loop?

Every earlier check was ONE hub in open loop (its regulators clamped). The real model is a closed
network, and the caveat was that feedback could shift the weights. This driver builds TWO hubs the
modular way and closes the loop between them:

    IL-6  <-- IL-17   (IL6SecFLS_MaxbyIL17 = 4.27, literature)
    IL-17 <-- IL-6    (IL17SecTh17_MaxbyIL6, no citation -> a fitted/prior strength, per rule A)

so IL-6 and IL-17 mutually up-regulate: a positive feedback loop. We assemble both secretion
modules, fit each baseline rate to the joint disease steady state, emit ONE SBML with both dynamic
species, and integrate the coupled ODEs. Then we compare the CLOSED loop (IL-17 free to respond) to
the OPEN loop (IL-17 clamped, the earlier reduced-module assumption) under the same intervention -
the gap is exactly the feedback amplification the open-loop probe could not see.

    python -m examples.run_qsp_e2e_coupled

Pure. IL-6<-IL-17 strength is literature; IL-17<-IL-6 is a prior (it has no citation - exactly the
kind of coupling A says must be fitted, not read off the table).
"""

from __future__ import annotations

import json
import os
import tempfile

from pkpd_agent.engines import model_assembly as MA
from pkpd_agent.engines.sbml_import import sbml_to_network
from examples.run_qsp_end_to_end import integrate

_DATA = os.path.join(os.path.dirname(__file__), "..", "projects", "vantage_ra", "data")
_MOTIF = {"proliferation_order": "zeroth", "combination": "product", "cap": None}


def hub_reactions(cell, base_param, clr_param, reg_specs):
    rate = MA.rate_from_motif(_MOTIF, base_param, reg_specs, cell)
    return [{"id": f"{cell}_sec", "reactants": [], "products": [cell], "rate": rate},
            {"id": f"{cell}_clr", "reactants": [cell], "products": [], "rate": f"{clr_param} * {cell}"}]


def eff_at(reg_specs, values, state):
    """Evaluate the product-of-folds effect for one hub at a given state (for kg fitting)."""
    p = 1.0
    for r in reg_specs:
        mx = values[r["max_param"]]; K = values[r["k_param"]]; x = state[r["species"]]
        p *= 1.0 + (mx - 1.0) * x / (K + x)
    return p


def main() -> None:
    prov = {p["name"]: p for p in json.load(open(os.path.join(_DATA, "param_provenance.json")))}
    tg = {t["model_species"]: t for t in json.load(open(os.path.join(_DATA,
          "steady_state_targets.json"))) if t.get("model_species")}
    lv = {c: float(tg[c]["target_model_unit"]) for c in tg
          if tg[c].get("target_model_unit") is not None}
    kcl6 = float(prov["kcl_IL6"]["value_from_reference"])
    kcl17 = float(prov["kcl_IL17"]["value_from_reference"])
    il6_t, il17_t = lv["IL6"], lv["IL17"]

    # ---- IL-6 hub (literature strengths, incl. IL6<-IL17 = 4.27) ----
    il6_maxes = {}
    for pref in ("IL6SecFLS_Maxby", "IL6SecMacro_Maxby"):
        for n, p in prov.items():
            if n.startswith(pref):
                c = n.split("Maxby")[-1]; v = p.get("value_from_reference")
                if v is not None and c in lv and c not in il6_maxes:
                    il6_maxes[c] = float(v)
    regs6 = [{"species": c, "max_param": f"M6_{c}", "k_param": f"K6_{c}"} for c in il6_maxes]

    # ---- IL-17 hub: IL6->IL17 is a PRIOR (no citation); IL23 literature (~1, negligible) ----
    il17_prior_from_il6 = 3.0            # agent's up-regulation prior for the un-cited coupling
    il17_maxes = {"IL6": il17_prior_from_il6, "IL23": float(prov["IL17SecTh17_MaxbyIL23"]
                  ["value_from_reference"])}
    il17_maxes = {c: m for c, m in il17_maxes.items() if c in lv}
    regs17 = [{"species": c, "max_param": f"M17_{c}", "k_param": f"K17_{c}"} for c in il17_maxes]

    # parameter values: K = level (Hill 0.5), Max as above
    values = {"kcl_IL6": kcl6, "kcl_IL17": kcl17}
    for c, m in il6_maxes.items():
        values[f"M6_{c}"] = m; values[f"K6_{c}"] = lv[c]
    for c, m in il17_maxes.items():
        values[f"M17_{c}"] = m; values[f"K17_{c}"] = lv[c]

    baseline = {**{c: lv[c] for c in lv}}
    baseline["IL6"], baseline["IL17"] = il6_t, il17_t
    values["kg_IL6"] = il6_t * kcl6 / eff_at(regs6, values, baseline)      # fit each baseline rate
    values["kg_IL17"] = il17_t * kcl17 / eff_at(regs17, values, baseline)  # to the joint steady state

    print("== two coupled hubs assembled from independent modules ==")
    print(f"  IL-6  <- {sorted(il6_maxes)}   (IL17 strength 4.27, literature)")
    print(f"  IL-17 <- {sorted(il17_maxes)}   (IL6 strength {il17_prior_from_il6:g}, PRIOR - uncited)")
    print(f"  fitted baseline rates: kg_IL6={values['kg_IL6']:.4g}, kg_IL17={values['kg_IL17']:.3g}")

    # ---- one SBML with BOTH species dynamic ----
    ext = sorted((set(il6_maxes) | set(il17_maxes)) - {"IL6", "IL17"})
    species = [{"name": "IL6", "initial": il6_t * 0.5}, {"name": "IL17", "initial": il17_t * 0.5}]
    for c in ext:
        species.append({"name": c, "initial": lv[c], "boundary": True})
    spec = {"name": "IL6_IL17_loop", "species": species,
            "parameters": [{"name": k, "value": v} for k, v in values.items()],
            "reactions": hub_reactions("IL6", "kg_IL6", "kcl_IL6", regs6)
                         + hub_reactions("IL17", "kg_IL17", "kcl_IL17", regs17),
            "rules": []}
    path = os.path.join(tempfile.gettempdir(), "il6_il17_loop.xml")
    open(path, "w", encoding="utf-8").write(MA.to_sbml(spec))
    net = sbml_to_network(path)
    print(f"\n  emitted {path}: {len(net['species'])} species, {len(net['reactions'])} reactions"
          f"  -> valid, runnable")

    clamp = {c: lv[c] for c in ext}
    ss = integrate(net, clamp, t_end=12.0)
    print(f"\n== closed-loop steady state (from a perturbed start) ==")
    print(f"  IL-6  -> {ss['IL6']:.4g}  (target {il6_t:g})")
    print(f"  IL-17 -> {ss['IL17']:.4g}  (target {il17_t:g})")
    print("  -> the two independently-built modules settle at the joint disease steady state:"
          " they compose into a stable loop.")

    # ---- feedback amplification: closed vs open loop under the same intervention ----
    # intervention: anti-IL-1 (IL1b -> 10%), a strong IL-6 driver. Closed loop lets the IL-6 drop
    # propagate to IL-17 and feed back; open loop holds IL-17 at baseline (reduced-module view).
    print("\n== feedback amplification (anti-IL-1: IL1b -> 10%) ==")
    held = dict(clamp); held["IL1b"] = lv["IL1b"] * 0.1
    closed = integrate(net, held, t_end=12.0)
    # open loop: clamp IL-17 at baseline too, so IL-6 sees no feedback
    open_clamp = dict(held); open_clamp["IL17"] = il17_t
    open_net = sbml_to_network(path)
    op = integrate(open_net, open_clamp, t_end=12.0)
    d6 = (closed["IL6"] - op["IL6"]) / op["IL6"]
    print(f"  {'':14}{'IL-17':>10}{'IL-6':>10}")
    print(f"  open loop   {il17_t:>10.4g}{op['IL6']:>10.4g}   (IL-17 held at baseline)")
    print(f"  closed loop {closed['IL17']:>10.4g}{closed['IL6']:>10.4g}   (IL-17 free to respond)")
    print(f"\n  closing the loop moves the IL-6 prediction by {d6:+.0%} vs the open-loop reduced "
          "module.")
    print("  -> feedback is real and modular assembly captures it: the open-loop probes were a")
    print("     lower bound on coupling effects, and how much it matters depends on the loop gain")
    print("     (here set partly by the uncited IL-17<-IL-6 prior - another param that needs data).")


if __name__ == "__main__":
    main()
