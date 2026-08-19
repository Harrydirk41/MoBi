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


def _system_prompt(target: dict, numeric: bool) -> str:
    head = (
        "You are a QSP modeler GENERATING a virtual population for a rheumatoid-"
        "arthritis model with the SimBiology engine. The model structure is fixed. "
        "You choose which disease-driver parameters to vary and over what bounds; "
        "each candidate patient draws those parameters, is simulated to its untreated "
        "disease baseline, and its baseline DAS28-CRP is recorded. The goal is a "
        f"population whose baseline DAS28-CRP matches the clinical target (mean "
        f"{target['mean']}, sd {target['sd']}) inside the active band {target['band']}.\n\n"
        "Start with vpop_inspect (the disease drivers - pro-inflammatory "
        "amplification factors and cell-baseline growth rates - and the target). "
        "These factors span orders of magnitude, so sample on a LOG scale, and vary "
        "several for genuine phenotypic diversity.\n\n")
    if numeric:
        body = (
            "BUILD THE POPULATION NUMERICALLY. Matching a distribution is a numerical "
            "selection problem, not something to hand-tune bounds toward:\n"
            "1. Choose the drivers and set WIDE log bounds that SPAN the target range "
            "(cover both mild and severe, so the pool brackets the target mean).\n"
            "2. Call vpop_select - a prevalence-weighting routine (the standard QSP "
            "method, what the paper's GA implements) reweights the pool to the target "
            "moments and returns the reweighted mean/sd and the EFFECTIVE SAMPLE SIZE.\n"
            "3. Judge the result: a healthy effective sample size means a diverse, "
            "on-target population; a LOW one means the pool barely covers the target "
            "(the reweighting is leaning on a few patients) - widen the bounds and "
            "re-select. Use vpop_sample only to probe where raw bounds land.\n"
            "4. Call vpop_finalize with your committed design, then finish with the "
            "drivers you varied, the bounds, the reweighted population, and what the "
            "effective sample size says about its diversity.")
    else:
        body = (
            "Build it by tuning the sampling bounds (no numerical selector available):\n"
            "1. Call vpop_sample with {bounds:{param:[lo,hi,scale]}}; read the yield "
            "(fraction in band), the accepted mean/sd, and the distance to target.\n"
            "2. Iterate: mean too low -> raise the inflammatory factors; spread too "
            "small -> widen; yield poor -> pull the extremes in. Balance yield vs "
            "diversity.\n"
            "3. Call vpop_finalize with your committed design, then finish with which "
            "parameters you varied and why, the bounds, and the population produced.")
    return head + body


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sbproj", required=True)
    ap.add_argument("--n", type=int, default=60, help="candidates per sample (default 60)")
    ap.add_argument("--baseline-day", type=float, default=200.0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--selector", choices=["numeric", "agent", "both"], default="both",
                    help="numeric = agent designs, prevalence-weighting selects "
                         "(vpop_select); agent = hand-tune bounds (vpop_sample only); "
                         "both = agent has both tools")
    ap.add_argument("--n-pool", type=int, default=80,
                    help="numeric selection: candidate pool size (default 80)")
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

        enable_select = args.selector in ("numeric", "both")
        registry = ToolRegistry()
        register_ra_vpop_loop_tools(registry, cfg, {
            "sb": sb, "n_samples": args.n, "baseline_day": args.baseline_day,
            "seed": args.seed, "target": target, "enable_select": enable_select,
            "n_pool": args.n_pool})

        goal = ("Generate a virtual population whose baseline DAS28-CRP matches the "
                f"clinical target {target}. Start with vpop_inspect, then build it.")
        policy = LLMPolicy(cfg, registry, _system_prompt(target, numeric=enable_select))
        loop = DecisionLoop(config=cfg, registry=registry, policy=policy)

        print(f"\n== LLM RA vpop-generation loop [selector={args.selector}] "
              f"(target {target}) ==\n", flush=True)

        def show(ev):
            if isinstance(ev, Decision):
                if ev.text:
                    print(f"\n[reason] {ev.text[:2000]}")
                for c in ev.calls:
                    print(f"  -> {c.name} "
                          f"{json.dumps(c.arguments, ensure_ascii=False)[:300]}")
                    if c.name in ("vpop_sample", "vpop_select"):
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
            via = final.get("via", "raw")
            print(f"target distribution: mean {target['mean']} sd {target['sd']} "
                  f"band {target['band']}")
            print(f"committed bounds ({via}): "
                  f"{json.dumps(final['bounds'], ensure_ascii=False)}")
            if via == "numeric":
                print(f"reweighted population: pool {s.get('n_pool')}, "
                      f"{s.get('n_inband')} in band; DAS28 mean {s.get('weighted_mean')} "
                      f"sd {s.get('weighted_sd')}; effective sample size "
                      f"{s.get('effective_sample_size')} ({int(s.get('ess_fraction',0)*100)}%)")
            else:
                print(f"population: {s.get('n_accepted')}/{s.get('n')} in band "
                      f"({s.get('yield_pct')}%), DAS28 mean {s.get('accepted_mean')} "
                      f"sd {s.get('accepted_sd')}")
            print(f"distribution distance to target: {s.get('distribution_distance')}")
    finally:
        sb.stop()


if __name__ == "__main__":
    main()
