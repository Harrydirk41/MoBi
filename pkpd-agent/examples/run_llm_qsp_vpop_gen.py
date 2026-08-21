r"""LLM-driven QSP virtual-population generation.

Model-agnostic (--model loads projects/<name>/tasks.json). The agent chooses which
disease-driver parameters to vary and over what bounds so the untreated baseline
severity distribution matches a clinical target. With the selector enabled it samples
a WIDE pool and a prevalence-weighting routine reweights it to the target moments.

    set ANTHROPIC_API_KEY=...
    python -m examples.run_llm_qsp_vpop_gen --model ra ^
        --sbproj "..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj" --max-steps 10
"""

from __future__ import annotations

import argparse
import json
import os

from pkpd_agent.config import AgentConfig
from pkpd_agent.engines.simbiology import SimBiologyEngine
from pkpd_agent.engines import qsp_config
from pkpd_agent.llm import LLMPolicy
from pkpd_agent.loop import DecisionLoop
from pkpd_agent.state import Decision, Finish, ModelingSession, Observation
from pkpd_agent.tools.registry import ToolRegistry
from pkpd_agent.tools.qsp_vpop_loop_tools import register_qsp_vpop_loop_tools


def _system_prompt(cfg, target, numeric: bool) -> str:
    sev = cfg.severity_readout
    method = ("Prefer vpop_select: sample WIDE bounds that span the target and let the "
              "prevalence-weighting routine reweight the pool - judge the result by the "
              "effective sample size. " if numeric else
              "Use vpop_sample to probe bounds and read the raw distribution. ")
    return (
        f"You are a QSP modeler building a virtual population for the {cfg.name}. "
        f"Choose which disease-driver parameters to vary and over what bounds so the "
        f"untreated baseline {sev} distribution matches the clinical target "
        f"{target}. Too-wide bounds push patients out of the active band; too-narrow "
        f"kills diversity. {method}Finalize the design you recommend.")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="vantage_ra")
    ap.add_argument("--sbproj", required=True)
    ap.add_argument("--n-samples", type=int, default=60)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--max-steps", type=int, default=10)
    ap.add_argument("--llm-model", dest="llm_model", default=None)
    ap.add_argument("--effort", default=None)
    ap.add_argument("--selector", choices=["on", "off"], default="on")
    args = ap.parse_args()

    tcfg = qsp_config.get(args.model)
    cfg = AgentConfig(mock=False, max_steps=args.max_steps)
    if args.llm_model:
        cfg.model = args.llm_model
    if args.effort:
        cfg.effort = args.effort
    if not cfg.anthropic_key_present():
        print("ANTHROPIC_API_KEY not set.")
        return

    target = tcfg.vpop_target
    enable_select = args.selector == "on"

    sb = SimBiologyEngine()
    try:
        print("== starting MATLAB/SimBiology engine ==", flush=True)
        sb.start()
        print(f"== loading {os.path.basename(args.sbproj)} ==", flush=True)
        sb.load_project(args.sbproj)

        registry = ToolRegistry()
        register_qsp_vpop_loop_tools(registry, cfg, {
            "cfg": tcfg, "sb": sb, "n_samples": args.n_samples, "seed": args.seed,
            "enable_select": enable_select})

        goal = (f"Build a virtual population whose baseline {tcfg.severity_readout} "
                f"matches the clinical target {target}. Start with vpop_inspect, then "
                "build it and finalize.")
        policy = LLMPolicy(cfg, registry, _system_prompt(tcfg, target, enable_select))
        loop = DecisionLoop(config=cfg, registry=registry, policy=policy)

        print(f"\n== LLM vpop-generation loop [selector={args.selector}] "
              f"(target {target}) ==\n", flush=True)

        def show(ev):
            if isinstance(ev, Decision):
                if ev.text:
                    print(f"\n[reason] {ev.text[:2000]}")
                for c in ev.calls:
                    print(f"  -> {c.name} "
                          f"{json.dumps(c.arguments, ensure_ascii=False)[:300]}")
                    if c.name in ("vpop_sample", "vpop_select"):
                        print("     ...sampling + simulating...", flush=True)
            elif isinstance(ev, Observation):
                print(f"  <- {ev.tool}: {ev.content.get('message', '')}", flush=True)
            elif isinstance(ev, Finish):
                print(f"\n=== AGENT SUMMARY ===\n{ev.text}")

        session = loop.run(goal, ModelingSession(goal=goal), on_event=show)

        final = session.get("vpop_final")
        print("\n== SCORING ==")
        if not final:
            print("agent did not finalize a design; best bounds:",
                  session.get("vpop_best_bounds"))
        else:
            s = final["score"]
            via = final.get("via", "raw")
            print(f"target distribution: mean {target['mean']} sd {target['sd']} "
                  f"band {target['band']}")
            print(f"committed bounds ({via}): {final['bounds']}")
            print(f"distribution distance to target: {s.get('distribution_distance')}")
    finally:
        sb.stop()


if __name__ == "__main__":
    main()
