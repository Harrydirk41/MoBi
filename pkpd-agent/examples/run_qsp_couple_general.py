r"""B) Coupled feedback loop - general, any mutually-regulating node pair, no hardcoding.

Discovers a FEEDBACK PAIR from the model itself: two buildable nodes X, Y where the model says X
regulates Y and Y regulates X (mutual regulation = a loop). Builds both secretion modules
independently, closes the loop, fits each baseline rate to the joint steady state, emits ONE SBML
with both species dynamic, integrates the coupled ODEs, and compares the CLOSED loop (both free) to
the OPEN loop (one clamped) under a shared upstream knockdown - the feedback amplification a
one-node probe cannot see.

    python -m examples.run_qsp_couple_general --model ra            # first discovered pair
    python -m examples.run_qsp_couple_general --model ra --pair IL6,IL17
    python -m examples.run_qsp_couple_general --model ra --list     # list all feedback pairs

Pure. Cross-strengths are looked up by the model's Maxby convention where cited, a prior otherwise.
"""

from __future__ import annotations

import argparse
import os
import tempfile

from pkpd_agent.engines import model_assembly as MA
from pkpd_agent.engines.sbml_import import sbml_to_network
from examples.run_qsp_end_to_end import integrate
from examples.run_qsp_build_general import load_model, discover_nodes, lookup_max

_MOTIF = {"proliferation_order": "zeroth", "combination": "product", "cap": None}


def feedback_pairs(nodes):
    pairs = []
    for x in nodes:
        for y in nodes[x]["regulators"]:
            if y in nodes and x in nodes[y]["regulators"] and x < y:
                pairs.append((x, y))
    return pairs


def hub_reactions(node, kg, kcl, reg_specs):
    rate = MA.rate_from_motif(_MOTIF, kg, reg_specs, node)
    return [{"id": f"{node}_sec", "reactants": [], "products": [node], "rate": rate},
            {"id": f"{node}_clr", "reactants": [node], "products": [], "rate": f"{kcl} * {node}"}]


def eff_prod(reg_specs, values, state):
    p = 1.0
    for r in reg_specs:
        p *= 1.0 + (values[r["max_param"]] - 1.0) * state[r["species"]] / \
            (values[r["k_param"]] + state[r["species"]])
    return p


def build_hub(node, prov, levels, nodes):
    """One node's secretion module (its model regulators; cross-strength prior if uncited)."""
    specs, vals = [], {}
    for c in nodes[node]["regulators"]:
        mx = lookup_max(prov, node, c)
        if mx is None:
            mx = 3.0                                        # uncited coupling -> up-prior
        specs.append({"species": c, "max_param": f"M_{node}_{c}", "k_param": f"K_{node}_{c}"})
        vals[f"M_{node}_{c}"] = mx; vals[f"K_{node}_{c}"] = levels[c]
    kcl = float(prov[nodes[node]["clearance"]]["value_from_reference"])
    return specs, vals, kcl


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="ra")
    ap.add_argument("--pair", default=None, help="X,Y; default = first discovered feedback pair")
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()

    prov, levels, cytokines = load_model(args.model)
    nodes = discover_nodes(prov, levels, cytokines)
    pairs = feedback_pairs(nodes)
    if args.list or not pairs:
        print(f"== feedback pairs discovered in '{args.model}' (mutual regulation) ==")
        for x, y in pairs:
            print(f"  {x} <-> {y}")
        if not pairs:
            print("  (none)")
        return

    if args.pair:
        a, b = [s.strip() for s in args.pair.split(",")]
    else:
        a, b = pairs[0]
    if (a, b) not in pairs and (b, a) not in pairs:
        print(f"{a}<->{b} is not a feedback pair; discovered: {pairs}"); return
    print(f"== coupled feedback loop: {a} <-> {b} (auto-discovered) ==")

    sa, va, kcla = build_hub(a, prov, levels, nodes)
    sb, vb, kclb = build_hub(b, prov, levels, nodes)
    values = {**va, **vb, f"kcl_{a}": kcla, f"kcl_{b}": kclb}
    baseline = {c: levels[c] for c in levels}
    values[f"kg_{a}"] = levels[a] * kcla / eff_prod(sa, values, baseline)
    values[f"kg_{b}"] = levels[b] * kclb / eff_prod(sb, values, baseline)
    print(f"  {a} <- {[s['species'] for s in sa]};  {b} <- {[s['species'] for s in sb]}")
    print(f"  fitted baseline rates kg_{a}={values[f'kg_{a}']:.3g}, kg_{b}={values[f'kg_{b}']:.3g}")

    ext = sorted({s["species"] for s in sa + sb} - {a, b})
    species = [{"name": a, "initial": levels[a] * 0.5}, {"name": b, "initial": levels[b] * 0.5}]
    for c in ext:
        species.append({"name": c, "initial": levels[c], "boundary": True})
    spec = {"name": f"{a}_{b}_loop", "species": species,
            "parameters": [{"name": k, "value": v} for k, v in values.items()],
            "reactions": hub_reactions(a, f"kg_{a}", f"kcl_{a}", sa)
                         + hub_reactions(b, f"kg_{b}", f"kcl_{b}", sb), "rules": []}
    path = os.path.join(tempfile.gettempdir(), f"loop_{a}_{b}.xml")
    open(path, "w", encoding="utf-8").write(MA.to_sbml(spec))
    net = sbml_to_network(path)

    clamp = {c: levels[c] for c in ext}
    ss = integrate(net, clamp, t_end=12.0)
    print(f"\n  closed-loop steady state: {a}={ss[a]:.4g} (target {levels[a]:g}), "
          f"{b}={ss[b]:.4g} (target {levels[b]:g})")

    # feedback amplification: knock the strongest EXTERNAL driver of a; compare closed vs open (b clamped)
    driver = max(ext, key=lambda c: abs((lookup_max(prov, a, c) or 1) - 1)) if ext else None
    if driver is None:
        print("  (no external driver to perturb)"); return
    held = dict(clamp); held[driver] = levels[driver] * 0.1
    closed = integrate(net, held, t_end=12.0)
    open_clamp = dict(held); open_clamp[b] = levels[b]
    op = integrate(sbml_to_network(path), open_clamp, t_end=12.0)
    d = (closed[a] - op[a]) / op[a] if op[a] else 0.0
    print(f"\n== feedback amplification (anti-{driver}) ==")
    print(f"  open loop ({b} clamped):  {a}={op[a]:.4g}")
    print(f"  closed loop ({b} free):   {a}={closed[a]:.4g},  {b}={closed[b]:.4g}")
    print(f"  -> closing the loop moves {a} by {d:+.0%} vs the open-loop reduced module: the "
          "feedback\n     the one-node probes bounded from below, on an auto-discovered pair, same "
          "code.")


if __name__ == "__main__":
    main()
