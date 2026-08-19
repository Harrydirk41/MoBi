r"""LLM-driven RA Stage-3: the agent GENERATES a virtual population.

Where run_llm_ra_fit calibrates a parameter, this drives the agent over virtual-
population generation: it chooses which disease-driver parameters to vary and over
what bounds, samples a candidate cohort, simulates each patient to the untreated
disease baseline, and matches the baseline DAS28-CRP distribution to the active-RA
clinical target. Too-wide bounds spill patients out of the active band; too-narrow
kills diversity. Finding the bounds that reproduce a realistic, diverse population
is the task.

    set ANTHROPIC_API_KEY=...
    python -m examples.run_llm_ra_vpop_gen ^
        --sbproj "..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj" ^
        --n 60 --max-steps 12

Each vpop_sample simulates the whole candidate cohort, so keep --n modest.
"""

from __future__ import annotations

import argparse
import json
import os

from pkpd_agent.config import AgentConfig
from pkpd_agent.engines.simbiology import SimBiologyEngine
from pkpd_agent.engines import osp_ra_trial
from pkpd_agent.llm import LLMPolicy
from pkpd_agent.loop import DecisionLoop
from pkpd_agent.state import Decision, Finish, ModelingSession, Observation
from pkpd_agent.tools.registry import ToolRegistry
from pkpd_agent.tools.ra_vpop_loop_tools import register_ra_vpop_loop_tools


def _system_prompt(target: dict) -> str:
    return (
        "You are a QSP modeler GENERATING a virtual population for a rheumatoid-"
        "arthritis model with the SimBiology engine. The model structure is fixed. "
        "You choose which disease-driver parameters to vary and over what bounds; "
        "each candidate patient draws those parameters, is simulated to its "
        "untreated disease baseline, and its baseline DAS28-CRP is recorded. Your "
        "goal is a population whose baseline DAS28-CRP matches the clinical target "
        f"(mean {target['mean']}, sd {target['sd']}) and mostly falls in the active-"
        f"disease band {target['band']}.\n\n"
        "Work like a modeler building a Vpop:\n"
        "1. Call vpop_inspect: the disease-driver parameters (meaning, nominal, "
        "observed span), and the target distribution.\n"
        "2. Reason about what drives severity: the pro-inflammatory amplification "
        "factors (TNF, IL-6, IL-17, ...) and cell-baseline growth rates push DAS28 "
        "up; these factors span orders of magnitude, so sample them on a LOG scale.\n"
        "3. Call vpop_sample with {bounds:{param:[lo,hi,scale]}}. Read the yield "
        "(fraction inside the active band), the accepted mean/sd, and the distance "
        "to target. You need enough drivers over wide enough bounds for phenotypic "
        "DIVERSITY (a realistic sd, not a spike), but not so wide that most patients "
        "fall outside the active band (low yield).\n"
        "4. Iterate the bounds: if the mean is too low, widen/raise the inflammatory "
        "factors; if the spread is too small, widen the bounds; if yield is poor, "
        "pull the extremes in. Balance yield against diversity.\n"
        "5. Call vpop_finalize with your committed design, then finish with: which "
        "parameters you varied and why, the bounds, and the population you produced."
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sbproj", required=True)
    ap.add_argument("--n", type=int, default=60, help="candidates per sample (default 60)")
    ap.add_argument("--baseline-day", type=float, default=200.0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--max-steps", type=int, default=12)
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

    target = osp_ra_trial.VPOP_TARGET

    sb = SimBiologyEngine()
    try:
        print("== starting MATLAB/SimBiology engine ==", flush=True)
        sb.start()
        print(f"== loading {os.path.basename(args.sbproj)} ==", flush=True)
        sb.load_project(args.sbproj)

        registry = ToolRegistry()
        register_ra_vpop_loop_tools(registry, cfg, {
            "sb": sb, "n_samples": args.n, "baseline_day": args.baseline_day,
            "seed": args.seed, "target": target})

        goal = ("Generate a virtual population whose baseline DAS28-CRP matches the "
                f"clinical target {target}. Start with vpop_inspect, then set "
                "sampling bounds and iterate.")
        policy = LLMPolicy(cfg, registry, _system_prompt(target))
        loop = DecisionLoop(config=cfg, registry=registry, policy=policy)

        print(f"\n== LLM RA vpop-generation loop (n={args.n}/sample, "
              f"target {target}) ==\n", flush=True)

        def show(ev):
            if isinstance(ev, Decision):
                if ev.text:
                    print(f"\n[reason] {ev.text[:2000]}")
                for c in ev.calls:
                    print(f"  -> {c.name} "
                          f"{json.dumps(c.arguments, ensure_ascii=False)[:300]}")
                    if c.name == "vpop_sample":
                        print("     ...simulating the candidate cohort...", flush=True)
            elif isinstance(ev, Observation):
                print(f"  <- {ev.tool}: {ev.content.get('message', '')}", flush=True)
            elif isinstance(ev, Finish):
                print(f"\n=== AGENT SUMMARY ===\n{ev.text}")

        session = loop.run(goal, ModelingSession(goal=goal), on_event=show)

        print("\n== SCORING [vpop-generation] ==")
        final = session.get("vpop_final")
        if not final:
            print("[warn] agent did not vpop_finalize; best design was",
                  session.get("vpop_best_bounds"))
        else:
            s = final["score"]
            print(f"target distribution: mean {target['mean']} sd {target['sd']} "
                  f"band {target['band']}")
            print(f"committed bounds: {json.dumps(final['bounds'], ensure_ascii=False)}")
            print(f"population: {s.get('n_accepted')}/{s.get('n')} in band "
                  f"({s.get('yield_pct')}%), DAS28 mean {s.get('accepted_mean')} "
                  f"sd {s.get('accepted_sd')}")
            print(f"distribution distance to target: {s.get('distribution_distance')}")
    finally:
        sb.stop()


if __name__ == "__main__":
    main()
