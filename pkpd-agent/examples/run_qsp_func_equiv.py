r"""Functional equivalence vs topological similarity, on a load-bearing HUB node (IL-6 secretion).

The topology benchmark scores a reconstructed model by how much its EDGE SET matches the paper's
graph (precision/recall). But 'different from the paper' is not the same as 'wrong': the model is
sloppy and non-unique, so a differently-wired model can still reproduce the same input->output
behaviour. This experiment measures the thing that actually matters - FUNCTIONAL equivalence - and
shows the two metrics dissociate.

Key idea (pure, no LLM, no fabrication): the functional cost of a topological difference does NOT
depend on guessing what an LLM would pick. It is a property of the model - a missed (or spurious)
regulator only changes behaviour where that regulator is ACTIVE at the operating point you care
about. So for each true IL-6 regulator we build the model that OMITS exactly that edge (an agent
that missed it), refit the baseline, then perturb that pathway at a held-out operating point and
measure how far IL-6 moves versus the full paper model. Symmetrically we price a SPURIOUS edge
(over-inclusion). The result is a per-edge table: which structural differences are functionally
free and which are catastrophic.

    python -m examples.run_qsp_func_equiv                 # pure, computes the cost table
    ANTHROPIC_API_KEY=... python -m examples.run_qsp_func_equiv --agent   # also runs the real
                                                                          # agent and prices ITS
                                                                          # actual chosen structure

Reads only the two provenance/target JSONs. Both the 'agent' and 'paper' reduced modules use the
IDENTICAL functional form (product of Hill fold-changes), so the ONLY thing that differs is the
regulator SET - isolating the cost of the topological difference itself.
"""

from __future__ import annotations

import argparse
import json
import os

_DATA = os.path.join(os.path.dirname(__file__), "..", "projects", "vantage_ra", "data")


def fold(mx: float, L: float, K: float) -> float:
    """One regulator's Hill fold-change on a secretion rate. mx>1 up-regulates, mx<1 down."""
    return 1.0 + (mx - 1.0) * L / (K + L)


def il6(regs, maxes, levels, ks, kg, kcl) -> float:
    """Reduced IL-6 steady state: baseline secretion kg times the product of fold-changes,
    over clearance. Same form for every model - only the regulator set `regs` changes."""
    prod = 1.0
    for c in regs:
        prod *= fold(maxes[c], levels[c], ks[c])
    return kg * prod / kcl


def fit_kg(regs, maxes, levels, ks, kcl, target) -> float:
    """Pin baseline secretion so this model reproduces the disease steady-state IL-6 exactly.
    Every model (full, edge-dropped, edge-added) is refit here, so all agree at baseline by
    construction - divergence can only appear at a held-out operating point."""
    prod = 1.0
    for c in regs:
        prod *= fold(maxes[c], levels[c], ks[c])
    return target * kcl / prod


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--agent", action="store_true",
                    help="also call the real LLM to choose IL-6 regulators and price its structure")
    ap.add_argument("--drop", type=float, default=0.1,
                    help="held-out therapy drops the targeted cytokine to this fraction")
    args = ap.parse_args()

    prov = {p["name"]: p for p in json.load(open(os.path.join(_DATA, "param_provenance.json")))}
    tg = {t["model_species"]: t for t in json.load(open(os.path.join(_DATA,
          "steady_state_targets.json"))) if t.get("model_species")}

    # ---- the paper's TRUE IL-6 regulators (union of the FLS and macrophage secretion terms) ----
    truth_maxes = {}
    for pref in ("IL6SecFLS_Maxby", "IL6SecMacro_Maxby"):
        for n, p in prov.items():
            if n.startswith(pref):
                c = n.split("Maxby")[-1]
                v = p.get("value_from_reference")
                if v is not None and c not in truth_maxes:       # first real value wins
                    truth_maxes[c] = float(v)
    levels = {c: float(tg[c]["target_model_unit"]) for c in tg
              if tg[c].get("target_model_unit") is not None}
    truth = [c for c in truth_maxes if c in levels]              # need a measured level to include
    ks = {c: levels[c] for c in levels}                          # K = level (Hill 0.5): one fit
    kcl = float(prov["kcl_IL6"]["value_from_reference"])
    target = float(tg["IL6"]["target_model_unit"])

    print("== load-bearing hub: IL-6 secretion ==")
    print(f"  paper's true regulators (max fold-change):")
    for c in sorted(truth, key=lambda c: -abs(truth_maxes[c] - 1)):
        d = "up" if truth_maxes[c] > 1 else "down"
        kind = "canonical" if c in ("TNFa", "IL1b", "IL17") else "phenomenological/weak"
        print(f"    {c:6} x{truth_maxes[c]:<6.3g} {d:4} ({kind})")

    kg_full = fit_kg(truth, truth_maxes, levels, ks, kcl, target)
    base_full = il6(truth, truth_maxes, levels, ks, kg_full, kcl)
    print(f"\n  full paper module reproduces baseline IL-6 = {base_full:.3g} (target {target:g})")

    # ---------- price each MISSED edge: an agent that omitted exactly one true regulator ----------
    print(f"\n== functional cost of a MISSED edge (held-out: anti-<cytokine> drops it to "
          f"{args.drop:g}x) ==")
    print(f"  {'missed':7} {'topo':>5}  {'paper IL-6':>11} {'agent IL-6':>11} {'func diff':>10}")
    rows = []
    for miss in sorted(truth, key=lambda c: -abs(truth_maxes[c] - 1)):
        sub = [c for c in truth if c != miss]
        kg_sub = fit_kg(sub, truth_maxes, levels, ks, kcl, target)   # refit to same baseline
        held = dict(levels); held[miss] = levels[miss] * args.drop   # therapy hits the missed path
        paper = il6(truth, truth_maxes, held, ks, kg_full, kcl)      # full model responds
        agent = il6(sub, truth_maxes, held, ks, kg_sub, kcl)         # agent's model can't respond
        diff = abs(agent - paper) / max(agent, paper)
        rows.append((miss, diff))
        print(f"  {miss:7} {'-1edge':>5}  {paper:>11.3g} {agent:>11.3g} {diff:>9.0%}")
    cheap = [c for c, d in rows if d < 0.05]
    dear = [c for c, d in rows if d >= 0.20]
    print(f"\n  -> missing {dear} is functionally CATASTROPHIC; missing {cheap} is nearly FREE.")
    print("     Same 'one missing edge' topological error; the functional cost ranges "
          f"{min(d for _,d in rows):.0%}..{max(d for _,d in rows):.0%} depending purely on whether")
    print("     the missed regulator is active & strong at the operating point you test.")

    # ---------- price a SPURIOUS edge: an agent that added a distractor not in the model ----------
    print(f"\n== functional cost of a SPURIOUS edge (over-inclusion), held-out perturbs the "
          f"distractor ==")
    print(f"  {'added':7} {'prior':>6}  {'paper IL-6':>11} {'agent IL-6':>11} {'func diff':>10}")
    for d, mx in [("VEGF", 1.5), ("GMCSF", 1.5), ("IL12", 1.5)]:
        if d not in levels:
            continue
        ext = truth + [d]
        em, ek = dict(truth_maxes), dict(ks)
        em[d] = mx; ek[d] = levels[d]
        kg_ext = fit_kg(ext, em, levels, ek, kcl, target)
        held = dict(levels); held[d] = levels[d] * args.drop        # therapy hits the fake pathway
        paper = il6(truth, truth_maxes, held, ks, kg_full, kcl)     # true model: no such edge, flat
        agent = il6(ext, em, held, ek, kg_ext, kcl)                 # agent's model responds falsely
        diff = abs(agent - paper) / max(agent, paper)
        print(f"  {d:7} x{mx:<5.3g}  {paper:>11.3g} {agent:>11.3g} {diff:>9.0%}")
    print("  -> a spurious edge is silent until you operate on it; then it fabricates a response.")

    print("\n== verdict ==")
    print("  Topological similarity and functional equivalence are DIFFERENT metrics. An agent")
    print("  model that scores low on edge-overlap is NOT necessarily wrong: its differences are")
    print("  free wherever they touch weak/inactive regulators, and only bite where they touch a")
    print("  strong regulator that moves at the clinical operating point of interest. So 'the")
    print("  agent's own-understanding model may well be right' is TRUE - qualified by exactly")
    print("  WHERE it differs, not by how many edges differ.")

    if args.agent:
        cyts = sorted(c for c, t in tg.items()
                      if t.get("kind") == "cytokine" and c in levels)
        _price_real_agent(cyts, truth, truth_maxes, levels, ks, kcl, target, args.drop)


def _price_real_agent(cyts, truth, truth_maxes, levels, ks, kcl, target, drop) -> None:
    """When a key is present: let the real LLM choose IL-6 regulators from biology, then price
    ITS actual structure (missed + spurious edges) on the same held-out sweep. No fabrication -
    this path only runs when the model actually answers."""
    from pkpd_agent.config import AgentConfig
    from pkpd_agent.engines import llm_tasks as LT, model_assembly as MA
    cfg = AgentConfig(mock=False)
    if not cfg.anthropic_key_present():
        print("\n[--agent] ANTHROPIC_API_KEY not set; skipping the real-agent pricing."); return
    print(f"\n== real agent chooses IL-6 secretion regulators (from {len(cyts)} cytokines) ==")
    regs = MA.propose_regulators("IL6", cyts, "secretion", LT.default_call(cfg))
    chosen = [r["cytokine"] for r in regs]
    dirn = {r["cytokine"]: (r.get("direction") or "up") for r in regs}
    missed = [c for c in truth if c not in chosen]
    spurious = [c for c in chosen if c not in truth]
    prec = len(set(chosen) & set(truth)) / len(chosen) if chosen else 0
    rec = len(set(chosen) & set(truth)) / len(truth) if truth else 0
    print(f"  chose: {chosen}")
    print(f"  vs truth {sorted(truth)}: recall {rec:.2f}, precision {prec:.2f}; "
          f"missed {missed}, spurious {spurious}")
    # Build the agent's ACTUAL reduced module: a real value where it can look one up (a true
    # regulator, so the table pins it and even self-corrects a wrong direction guess); for a
    # spurious edge there is no table entry, so it falls back to its OWN stated direction as a
    # prior (up -> 1.5, down -> 0.6).
    akeys = [c for c in chosen if c in levels]
    am = {c: (truth_maxes[c] if c in truth_maxes
              else (1.5 if dirn.get(c) == "up" else 0.6)) for c in akeys}
    ak = {c: levels[c] for c in akeys}
    kg_a = fit_kg(akeys, am, levels, ak, kcl, target)
    kg_full = fit_kg(truth, truth_maxes, levels, ks, kcl, target)
    # Worst-case held-out error over BOTH failure surfaces: perturbing each true pathway (exposes
    # a MISSED edge) AND each spurious pathway (exposes an OVER-INCLUDED edge that fabricates a
    # response the true model does not have). Skipping the spurious sweep under-reports the error.
    worst_name, worst = "", 0.0
    for pert in sorted(set(truth) | set(spurious)):
        if pert not in levels:
            continue
        held = dict(levels); held[pert] = levels[pert] * drop
        paper = il6(truth, truth_maxes, held, ks, kg_full, kcl)
        agent = il6(akeys, am, held, ak, kg_a, kcl)
        d = abs(agent - paper) / max(agent, paper)
        if d > worst:
            worst_name, worst = pert, d
    where = "spurious edge" if worst_name in spurious else "missed/mis-set edge"
    print(f"  worst-case held-out IL-6 error of the agent's ACTUAL chosen structure: "
          f"{worst:.0%} (perturbing {worst_name or 'n/a'}, a {where})")


if __name__ == "__main__":
    main()
