r"""Stage-1 benchmark: an LLM agent reconstructs the RA disease network from scratch.

Given only the cast of cells and cytokines, the agent proposes the regulatory wiring;
it is scored against the real Vantage RA model's edges. This is the one genuinely
Stage-1, genuinely reasoning-heavy task - no simulation, no run-and-compare, and the
edge space is far too large to brute-force.

The answer key:
  * BEST: run the MATLAB dump once to get the COMPLETE wiring, then point --network at it
        >> sb_load('...\Vantage RA QSP Model v1.0.sbproj'); sb_network_json('network.json')
    python -m examples.run_llm_ra_network --network network.json
  * BOOTSTRAP (no MATLAB): parse the partial key drawn in the SimBiology diagram
    python -m examples.run_llm_ra_network --sbproj "..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj"

No SimBiology engine is needed to RUN the benchmark - only ANTHROPIC_API_KEY.
"""

from __future__ import annotations

import argparse
import json

from pkpd_agent.config import AgentConfig
from pkpd_agent.engines import ra_network as N
from pkpd_agent.llm import LLMPolicy
from pkpd_agent.loop import DecisionLoop
from pkpd_agent.state import Decision, Finish, ModelingSession, Observation
from pkpd_agent.tools.registry import ToolRegistry
from pkpd_agent.tools.ra_network_loop_tools import register_ra_network_loop_tools


_CONVENTIONS = (
    "\n\nIMPORTANT - how THIS QSP model is structured (its modeling grammar, distinct "
    "from textbook cause-and-effect):\n"
    "1. The network is BIPARTITE cell<->cytokine. Cells SECRETE cytokines (cell->"
    "cytokine edges); cytokines DRIVE cell proliferation / influx / apoptosis (cytokine"
    "->cell edges). AVOID cytokine->cytokine shortcuts - where textbook biology says "
    "'TNFa induces IL-6', the model routes it through a cell (TNFa->cell, cell->IL6).\n"
    "2. Every secreting cell also carries a NEGATIVE self-feedback edge on its own "
    "cytokine (e.g. Macro -| its own IL-6/TNFa, Th1 -| its own IFN-g): saturating / "
    "self-limiting terms the model needs to stay bounded. Include these.\n"
    "3. Chemokines (MIP3, RANTES, MCP1) and adhesion recruit BROADLY - fan them out "
    "across many leukocyte compartments, not one.\n"
    "4. TGF-b and IL-10 are context-dependent, not globally suppressive: some targets "
    "positive, some negative - reason per target rather than assuming all-inhibitory."
)


def _system_prompt(conventions: bool = False) -> str:
    base = (
        "You are an immunologist reconstructing the regulatory network of a rheumatoid-"
        "arthritis QSP model. You are given the cast of cells and cytokines and must "
        "propose the directed, signed edges - which cell/cytokine up- or down-regulates "
        "the secretion, proliferation, or influx of which other node. Reason from "
        "established RA biology: macrophage and FLS as the main TNF-a / IL-6 / IL-1b "
        "sources; the Th17 / IL-17 axis; Th1 / IFN-g; regulatory (Treg, anti-"
        "inflammatory) suppression as negative edges; cytokine-driven cell influx and "
        "proliferation. Build the draft in batches with network_propose, use the "
        "structural feedback to find gaps, then network_finalize exactly once. Balance "
        "recall (find the real edges) against precision (do not propose every pair)."
    )
    return base + (_CONVENTIONS if conventions else "")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--network", help="network.json from sb_network_json.m (full key)")
    ap.add_argument("--sbproj", help="sbproj to parse the diagram key from (bootstrap)")
    ap.add_argument("--max-steps", type=int, default=20)
    ap.add_argument("--conventions", action="store_true",
                    help="prime the agent with the model's bipartite + negative-feedback "
                         "conventions (tests whether its gap was grammar, not biology)")
    ap.add_argument("--model", default=None)
    ap.add_argument("--effort", default=None)
    args = ap.parse_args()

    if args.network:
        truth = N.parse_truth(args.network)
        src = args.network
    elif args.sbproj:
        truth = N.parse_truth_from_diagram(args.sbproj)
        src = f"{args.sbproj} (diagram subset)"
    else:
        ap.error("give --network (full key) or --sbproj (diagram bootstrap)")
    print(f"answer key: {len(truth)} regulatory edges from {src}\n")

    cfg = AgentConfig(mock=False, max_steps=args.max_steps)
    if args.model:
        cfg.model = args.model
    if args.effort:
        cfg.effort = args.effort
    if not cfg.anthropic_key_present():
        print("ANTHROPIC_API_KEY not set.")
        return

    registry = ToolRegistry()
    register_ra_network_loop_tools(registry, cfg, {"truth": truth})

    goal = ("Reconstruct the RA disease network: propose the directed signed regulatory "
            "edges among the cells and cytokines, then finalize to be scored.")
    print(f"prompt: {'convention-primed' if args.conventions else 'biology-only (baseline)'}\n")
    policy = LLMPolicy(cfg, registry, _system_prompt(args.conventions))
    loop = DecisionLoop(config=cfg, registry=registry, policy=policy)

    def show(ev):
        if isinstance(ev, Decision):
            if ev.text:
                print(f"\n[reason] {ev.text[:1200]}")
            for c in ev.calls:
                n = len((c.arguments or {}).get("edges", [])) if c.name == "network_propose" else ""
                print(f"  -> {c.name} {('('+str(n)+' edges)') if n != '' else ''}")
        elif isinstance(ev, Observation):
            print(f"  <- {ev.tool}: {ev.content.get('message', '')}")
        elif isinstance(ev, Finish):
            print(f"\n=== AGENT SUMMARY ===\n{ev.text}")

    session = loop.run(goal, ModelingSession(goal=goal), on_event=show)

    final = session.get("net_final")
    print("\n== SCORE vs the real model ==")
    if not final:
        print("  (agent did not finalize)")
        return
    topo, sa = final["topology"], final["sign_aware"]
    print(f"  edges proposed : {final['n_edges']}   (truth {topo['n_truth']})")
    print(f"  TOPOLOGY  : P {topo['precision']}  R {topo['recall']}  F1 {topo['f1']}"
          f"   ({topo['hit']} found, {topo['missed']} missed, {topo['extra']} extra)")
    print(f"  SIGN-AWARE: P {sa['precision']}  R {sa['recall']}  F1 {sa['f1']}")
    print("\n  missed real edges (recall gaps):")
    for e in topo["missed_edges"][:25]:
        print(f"    {e}")
    print("\n  extra edges not in the model (may be defensible biology):")
    for e in topo["extra_edges"][:25]:
        print(f"    {e}")


if __name__ == "__main__":
    main()
