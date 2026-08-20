r"""Stage-1 Layer 5: an LLM estimates the RA model's parameter values from physiology.

Given each parameter's name, units, and cell context, the agent predicts a value; it is
scored by order-of-magnitude error against the model (ESM2), split into the easy
dimensionless fold-effects and the harder dimensional parameters (rates, secretion,
concentrations), and compared to a naive unit-geomean baseline. Beating that baseline on
the dimensional parameters is the real signal that the LLM contributes physiological
knowledge, not just units bookkeeping.

No SimBiology engine needed - the answer key is the committed ESM2 extract.

    python -m examples.run_llm_ra_params            # uses the bundled key
    python -m examples.run_llm_ra_params --max-steps 12
"""

from __future__ import annotations

import argparse

from pkpd_agent.config import AgentConfig
from pkpd_agent.engines import ra_params as P
from pkpd_agent.llm import LLMPolicy
from pkpd_agent.loop import DecisionLoop
from pkpd_agent.state import Decision, Finish, ModelingSession, Observation
from pkpd_agent.tools.registry import ToolRegistry
from pkpd_agent.tools.ra_params_loop_tools import register_ra_params_loop_tools


def _system_prompt() -> str:
    return (
        "You are a quantitative systems pharmacologist setting physiologically grounded "
        "priors for a rheumatoid-arthritis QSP model's parameters. For each parameter you "
        "are given its name, units, and cell context; predict a numeric value. Reason "
        "from physiology: turnover rates (1/day) from cell/cytokine half-lives, EC50-type "
        "concentrations from typical cytokine levels, and dimensionless Max-fold effects "
        "as modest multipliers around 1. You are scored on order-of-magnitude accuracy; "
        "the dimensional parameters (rates, secretion, concentrations) are where real "
        "physiological knowledge matters. Estimate every parameter with param_estimate "
        "(batches are fine), then param_finalize once."
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--key", default=None, help="override the ESM2 params JSON")
    ap.add_argument("--max-steps", type=int, default=14)
    ap.add_argument("--model", default=None)
    ap.add_argument("--effort", default=None)
    args = ap.parse_args()

    truth = P.load_truth(args.key) if args.key else P.load_truth()
    print(f"answer key: {len(truth)} model parameters (ESM2)\n")

    cfg = AgentConfig(mock=False, max_steps=args.max_steps)
    if args.model:
        cfg.model = args.model
    if args.effort:
        cfg.effort = args.effort
    if not cfg.anthropic_key_present():
        print("ANTHROPIC_API_KEY not set.")
        return

    registry = ToolRegistry()
    register_ra_params_loop_tools(registry, cfg, {"truth": truth})

    goal = ("Estimate every RA model parameter's value from its name, units, and cell "
            "context, then finalize to be scored order-of-magnitude vs the real model.")
    policy = LLMPolicy(cfg, registry, _system_prompt())
    loop = DecisionLoop(config=cfg, registry=registry, policy=policy)

    def show(ev):
        if isinstance(ev, Decision):
            if ev.text:
                print(f"\n[reason] {ev.text[:1000]}")
            for c in ev.calls:
                n = len((c.arguments or {}).get("predictions", [])) \
                    if c.name == "param_estimate" else ""
                print(f"  -> {c.name} {('('+str(n)+' values)') if n != '' else ''}")
        elif isinstance(ev, Observation):
            print(f"  <- {ev.tool}: {ev.content.get('message', '')}")
        elif isinstance(ev, Finish):
            print(f"\n=== AGENT SUMMARY ===\n{ev.text}")

    session = loop.run(goal, ModelingSession(goal=goal), on_event=show)

    final = session.get("param_final")
    print("\n== SCORE vs the real model (order of magnitude) ==")
    if not final:
        print("  (agent did not finalize)")
        return

    def line(tag, d):
        if d.get("n"):
            print(f"  {tag:16s} n={d['n']:3d}  median log10 err {d['median_log10_err']}"
                  f"   within 3x {d['within_3x']}  10x {d['within_10x']}")
    line("OVERALL", final["overall"])
    line("dimensionless", final["dimensionless"])
    line("DIMENSIONAL", final["dimensional"])
    line("naive baseline", final["naive_unit_geomean_baseline"])
    print(f"\n  beats naive baseline (overall median): {final['beats_baseline']}")
    print("\n  worst order-of-magnitude misses:")
    for w in final["worst_misses"][:10]:
        print(f"    {w['name']:32s} true {w['true']:<10g} pred {w['pred']:<10g} "
              f"({w['log10_err']} dex off)")


if __name__ == "__main__":
    main()
