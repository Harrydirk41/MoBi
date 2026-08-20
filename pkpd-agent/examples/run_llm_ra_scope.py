r"""Stage-2: an LLM chooses the RA model's SCOPE (which cells and mediators to include).

Given only the disease and the modeling goal, the agent proposes the cast; scored against
the Vantage RA model's real 26 nodes. Precision is the real signal - it measures whether
the LLM shares the model's parsimony or just dumps the RA textbook.

No SimBiology engine needed - the answer key is the model's node list.

    python -m examples.run_llm_ra_scope
    python -m examples.run_llm_ra_scope --repeat 5
"""

from __future__ import annotations

import argparse

from pkpd_agent.config import AgentConfig
from pkpd_agent.llm import LLMPolicy
from pkpd_agent.loop import DecisionLoop
from pkpd_agent.state import Decision, Finish, ModelingSession, Observation
from pkpd_agent.tools.registry import ToolRegistry
from pkpd_agent.tools.ra_scope_loop_tools import register_ra_scope_loop_tools


_CONVENTIONS = (
    "\n\nIMPORTANT - scope to the ENDPOINTS, not to RA pathology in general:\n"
    "1. ACR/DAS28-CRP read out tender/swollen joints + CRP/ESR (systemic inflammation) - "
    "NOT radiographic joint damage. So LEAVE OUT the bone/cartilage EROSION axis "
    "(chondrocytes, osteoclasts, RANKL, MMPs): the endpoints never see it.\n"
    "2. DO INCLUDE the cell-recruitment/trafficking layer - endothelium, adhesion "
    "molecules (CAM), and chemokines (MCP1, MIP3, RANTES): influx into the synovium is "
    "modeled mechanistically and drives the cell densities the endpoints depend on.\n"
    "3. CRP and IL-8 are not separate nodes (CRP is read off IL-6; IL-8/neutrophils are "
    "omitted). Autoantibodies and BAFF ARE in (B-cell-targeted therapy arm)."
)


def _system_prompt(conventions: bool = False) -> str:
    base = (
        "You are a QSP modeler scoping a rheumatoid-arthritis model. Choose the cast: the "
        "cell types and soluble mediators to include so the model can simulate late-phase "
        "trials and reproduce ACR / DAS28-CRP endpoints. A good QSP model is PARSIMONIOUS "
        "- include the cells and mediators that drive the modeled synovial biology and the "
        "endpoints (the main effector and regulatory leukocytes, the structural cells of "
        "the joint, and the cytokines/chemokines/growth factors central to RA and its "
        "therapies), and deliberately leave out the long tail of RA-associated molecules "
        "that would not change the trial-level behavior. Propose the cast with "
        "scope_propose, then scope_finalize once. Precision is scored - do not dump the "
        "whole immunology textbook."
    )
    return base + (_CONVENTIONS if conventions else "")


def _report_variance(tag: str, xs: list) -> None:
    n = len(xs)
    mean = sum(xs) / n
    sd = (sum((x - mean) ** 2 for x in xs) / (n - 1)) ** 0.5 if n > 1 else 0.0
    print(f"\n  {tag:16s} over {n} runs: mean {mean:.3f}  sd {sd:.3f}  "
          f"min {min(xs):.3f}  max {max(xs):.3f}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--max-steps", type=int, default=12)
    ap.add_argument("--repeat", type=int, default=1)
    ap.add_argument("--conventions", action="store_true",
                    help="prime the endpoint-focus + trafficking-layer conventions the "
                         "agent diagnosed missing (tests judgment vs recall)")
    ap.add_argument("--model", default=None)
    ap.add_argument("--effort", default=None)
    args = ap.parse_args()

    cfg = AgentConfig(mock=False, max_steps=args.max_steps)
    if args.model:
        cfg.model = args.model
    if args.effort:
        cfg.effort = args.effort
    if not cfg.anthropic_key_present():
        print("ANTHROPIC_API_KEY not set.")
        return

    registry = ToolRegistry()
    register_ra_scope_loop_tools(registry, cfg, {})
    goal = ("Choose the RA model's cast - the cells and mediators to include - then "
            "finalize to be scored against the real model.")
    verbose = args.repeat == 1

    def show(ev):
        if not verbose:
            return
        if isinstance(ev, Decision):
            if ev.text:
                print(f"\n[reason] {ev.text[:1200]}")
            for c in ev.calls:
                n = len((c.arguments or {}).get("nodes", [])) if c.name == "scope_propose" else ""
                print(f"  -> {c.name} {('('+str(n)+' names)') if n != '' else ''}")
        elif isinstance(ev, Observation):
            print(f"  <- {ev.tool}: {ev.content.get('message', '')}")
        elif isinstance(ev, Finish):
            print(f"\n=== AGENT SUMMARY ===\n{ev.text}")

    runs = []
    for i in range(args.repeat):
        policy = LLMPolicy(cfg, registry, _system_prompt(args.conventions))
        loop = DecisionLoop(config=cfg, registry=registry, policy=policy)
        session = loop.run(goal, ModelingSession(goal=goal), on_event=show)
        final = session.get("scope_final")
        if not final:
            print(f"  run {i + 1}: agent did not finalize")
            continue
        runs.append(final)
        if not verbose:
            print(f"  run {i + 1}/{args.repeat}: F1 {final['f1']} "
                  f"(P {final['precision']} R {final['recall']}), "
                  f"{final['hit']} hit / {final['extra']} extra")

    if not runs:
        print("\nno finalized runs.")
        return

    if verbose:
        f = runs[0]
        print("\n== SCOPE vs the real model (26 nodes) ==")
        print(f"  F1 {f['f1']}   precision {f['precision']}   recall {f['recall']}")
        print(f"  {f['hit']} of {f['hit'] + f['missed']} real nodes included; "
              f"{f['extra']} extra")
        print(f"\n  missed (real nodes left out): {f['missed_nodes']}")
        print(f"\n  extra (not in the model): {f['extra_nodes']}")
        print(f"  ...of which are real RA mediators the model excludes: "
              f"{f['extra_known_mediators']}")
    else:
        _report_variance("SCOPE F1", [r["f1"] for r in runs])
        _report_variance("precision", [r["precision"] for r in runs])
        _report_variance("recall", [r["recall"] for r in runs])


if __name__ == "__main__":
    main()
