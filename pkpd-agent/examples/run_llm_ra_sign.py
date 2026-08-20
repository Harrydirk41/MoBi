r"""Isolated sign prediction: given the model's true edges, can the LLM sign them?

The weakest layer (signs) tested on its own - the agent is handed the real (unsigned)
edges and only decides activate vs inhibit; scored vs the model against the majority-class
baseline. Uses the same answer key as the topology task.

    python -m examples.run_llm_ra_sign --network network.json
    python -m examples.run_llm_ra_sign --network network.json --repeat 5
    python -m examples.run_llm_ra_sign --sbproj "..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj"
"""

from __future__ import annotations

import argparse

from pkpd_agent.config import AgentConfig
from pkpd_agent.engines import ra_network as N
from pkpd_agent.llm import LLMPolicy
from pkpd_agent.loop import DecisionLoop
from pkpd_agent.state import Decision, Finish, ModelingSession, Observation
from pkpd_agent.tools.registry import ToolRegistry
from pkpd_agent.tools.ra_sign_loop_tools import register_ra_sign_loop_tools


def _system_prompt() -> str:
    return (
        "You are an immunologist annotating a rheumatoid-arthritis QSP model. You are given "
        "its real regulatory edges and must decide, for each, whether the source ACTIVATES "
        "(+1) or INHIBITS (-1) the target. Most edges are pro-inflammatory activations; the "
        "inhibitory ones are the anti-inflammatory mediators (TGF-b, IL-10, Treg) and the "
        "model's negative self-feedback loops where a cell down-regulates its own cytokine "
        "to stay bounded. Sign every edge with sign_predict, then sign_finalize once."
    )


def _report_variance(tag: str, xs: list) -> None:
    n = len(xs)
    mean = sum(xs) / n
    sd = (sum((x - mean) ** 2 for x in xs) / (n - 1)) ** 0.5 if n > 1 else 0.0
    print(f"\n  {tag:16s} over {n} runs: mean {mean:.3f}  sd {sd:.3f}  "
          f"min {min(xs):.3f}  max {max(xs):.3f}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--network", help="network.json (full key)")
    ap.add_argument("--sbproj", help="sbproj (diagram bootstrap key)")
    ap.add_argument("--max-steps", type=int, default=12)
    ap.add_argument("--repeat", type=int, default=1)
    ap.add_argument("--model", default=None)
    ap.add_argument("--effort", default=None)
    args = ap.parse_args()

    if args.network:
        truth = N.parse_truth(args.network)
    elif args.sbproj:
        truth = N.parse_truth_from_diagram(args.sbproj)
    else:
        ap.error("give --network or --sbproj")
    print(f"answer key: {len(truth)} edges to sign\n")

    cfg = AgentConfig(mock=False, max_steps=args.max_steps)
    if args.model:
        cfg.model = args.model
    if args.effort:
        cfg.effort = args.effort
    if not cfg.anthropic_key_present():
        print("ANTHROPIC_API_KEY not set.")
        return

    registry = ToolRegistry()
    register_ra_sign_loop_tools(registry, cfg, {"truth": truth})
    goal = "Sign every edge of the RA network, then finalize to be scored vs the model."
    verbose = args.repeat == 1

    def show(ev):
        if not verbose:
            return
        if isinstance(ev, Decision):
            if ev.text:
                print(f"\n[reason] {ev.text[:1000]}")
            for c in ev.calls:
                n = len((c.arguments or {}).get("edges", [])) if c.name == "sign_predict" else ""
                print(f"  -> {c.name} {('('+str(n)+' signs)') if n != '' else ''}")
        elif isinstance(ev, Observation):
            print(f"  <- {ev.tool}: {ev.content.get('message', '')}")
        elif isinstance(ev, Finish):
            print(f"\n=== AGENT SUMMARY ===\n{ev.text}")

    runs = []
    for i in range(args.repeat):
        policy = LLMPolicy(cfg, registry, _system_prompt())
        loop = DecisionLoop(config=cfg, registry=registry, policy=policy)
        session = loop.run(goal, ModelingSession(goal=goal), on_event=show)
        final = session.get("sign_final")
        if not final:
            print(f"  run {i + 1}: agent did not finalize")
            continue
        runs.append(final)
        if not verbose:
            print(f"  run {i + 1}/{args.repeat}: accuracy {final['accuracy']} "
                  f"({final['correct']}/{final['n']}), majority baseline "
                  f"{final['majority_baseline']}, beats {final['beats_majority']}")

    if not runs:
        print("\nno finalized runs.")
        return

    if verbose:
        f = runs[0]
        print("\n== SIGN accuracy vs the model ==")
        print(f"  accuracy {f['accuracy']} ({f['correct']}/{f['n']})")
        print(f"  majority-class baseline {f['majority_baseline']} "
              f"({f['frac_positive']} of edges activate)   beats majority: "
              f"{f['beats_majority']}")
    else:
        _report_variance("sign accuracy", [r["accuracy"] for r in runs])
        print(f"    ^ majority-class baseline = {runs[0]['majority_baseline']}  "
              f"<- the bar to beat")


if __name__ == "__main__":
    main()
