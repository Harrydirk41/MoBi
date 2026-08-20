r"""Does the LLM know which parameters matter? Rank sensitivity vs the paper's GSA.

From a pool of real model parameters (the GSA top-20 hidden among distractors), the agent
picks/ranks the most influential on DAS28-CRP; scored against the paper's global
sensitivity analysis (Fig 9), vs a blind-pick random baseline. No engine / answer key file
needed - the GSA top-20 is bundled.

    python -m examples.run_llm_ra_sensitivity
    python -m examples.run_llm_ra_sensitivity --repeat 5
"""

from __future__ import annotations

import argparse

from pkpd_agent.config import AgentConfig
from pkpd_agent.llm import LLMPolicy
from pkpd_agent.loop import DecisionLoop
from pkpd_agent.state import Decision, Finish, ModelingSession, Observation
from pkpd_agent.tools.registry import ToolRegistry
from pkpd_agent.tools.ra_sensitivity_loop_tools import register_ra_sensitivity_loop_tools


def _system_prompt() -> str:
    return (
        "You are a QSP modeler doing a global sensitivity analysis in your head. Given a "
        "pool of model parameters, identify which most drive variance in the DAS28-CRP "
        "disease-severity readout. Reason about leverage: parameters that set the size of "
        "the abundant, disease-driving cell populations (baseline growth/influx rates) and "
        "the fractional drivers of the dominant cytokine axes move the whole system; "
        "narrow downstream effect strengths and drug-binding constants move it less. Rank "
        "the ~20 most sensitive (most first) with sens_rank, then sens_finalize once."
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
    ap.add_argument("--max-steps", type=int, default=10)
    ap.add_argument("--repeat", type=int, default=1)
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
    register_ra_sensitivity_loop_tools(registry, cfg, {})
    goal = ("Rank the parameters most sensitive for DAS28-CRP from the pool, then finalize "
            "to be scored against the model's global sensitivity analysis.")
    verbose = args.repeat == 1

    def show(ev):
        if not verbose:
            return
        if isinstance(ev, Decision):
            if ev.text:
                print(f"\n[reason] {ev.text[:1200]}")
            for c in ev.calls:
                n = len((c.arguments or {}).get("ranked", [])) if c.name == "sens_rank" else ""
                print(f"  -> {c.name} {('('+str(n)+' ranked)') if n != '' else ''}")
        elif isinstance(ev, Observation):
            print(f"  <- {ev.tool}: {ev.content.get('message', '')}")
        elif isinstance(ev, Finish):
            print(f"\n=== AGENT SUMMARY ===\n{ev.text}")

    runs = []
    for i in range(args.repeat):
        policy = LLMPolicy(cfg, registry, _system_prompt())
        loop = DecisionLoop(config=cfg, registry=registry, policy=policy)
        session = loop.run(goal, ModelingSession(goal=goal), on_event=show)
        final = session.get("sens_final")
        if not final:
            print(f"  run {i + 1}: agent did not finalize")
            continue
        runs.append(final)
        if not verbose:
            print(f"  run {i + 1}/{args.repeat}: recall {final['recall']} "
                  f"({final['hit']}/20), precision {final['precision']}, "
                  f"spearman {final['spearman_on_hits']}, beats random {final['beats_random']}")

    if not runs:
        print("\nno finalized runs.")
        return

    if verbose:
        f = runs[0]
        print("\n== SENSITIVITY vs the paper's GSA (top-20) ==")
        print(f"  recall {f['recall']} ({f['hit']}/20)   precision {f['precision']}   "
              f"F1 {f['f1']}")
        print(f"  rank correlation on hits: {f['spearman_on_hits']}")
        print(f"  random blind-pick baseline recall: {f['random_baseline_recall']}   "
              f"beats random: {f['beats_random']}")
        print(f"\n  missed GSA top-20: {f['missed_top20']}")
        print(f"\n  picked distractors (not in GSA top-20): {f['picked_distractors']}")
    else:
        _report_variance("recall", [r["recall"] for r in runs])
        print(f"    ^ random blind-pick baseline recall ≈ "
              f"{runs[0]['random_baseline_recall']}  <- the bar to beat")
        _report_variance("precision", [r["precision"] for r in runs])
        sp = [r["spearman_on_hits"] for r in runs if r["spearman_on_hits"] is not None]
        if sp:
            _report_variance("rank corr (hits)", sp)


if __name__ == "__main__":
    main()
