r"""LLM-driven QSP calibration: fit a PD parameter so the model reproduces a trial.

Model-agnostic (--model loads projects/<name>/tasks.json). The protocol is FIXED (a
given drug arm); the unknown is a model PD parameter. The agent sets up a numerical
1-D fit and interprets it (identifiability, comparison to the literature reference).

    set ANTHROPIC_API_KEY=...
    python -m examples.run_llm_qsp_fit --model ra ^
        --sbproj "..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj" ^
        --vpop "..\RA-QSP-Model\Vpop1.xlsx" --limit 50 --max-steps 10
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
from pkpd_agent.tools.qsp_fit_loop_tools import register_qsp_fit_loop_tools


def _system_prompt(cfg) -> str:
    return (
        f"You are a QSP modeler calibrating the {cfg.name}. A drug arm is FIXED; you "
        "estimate a model PD parameter so the simulated response matches an OBSERVED "
        "clinical target. You have fit_inspect (the arm, the parameter, the target), "
        "fit_optimize (set up a bounded numerical fit - the optimizer minimizes), "
        "fit_try (probe one value), and fit_finalize (commit). Prefer the optimizer "
        "over hand-searching. Your job is choosing what to fit, its bounds and scale, "
        "and interpreting the result: is the parameter identifiable from the profile, "
        "are residuals structural, how does the fitted value compare to the "
        "literature? Finalize the parameter set you recommend.")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="ra")
    ap.add_argument("--sbproj", required=True)
    ap.add_argument("--vpop", required=True)
    ap.add_argument("--arm", default=None, help="fixed drug arm (default: project's)")
    ap.add_argument("--limit", type=int, default=50)
    ap.add_argument("--max-steps", type=int, default=10)
    ap.add_argument("--llm-model", dest="llm_model", default=None)
    ap.add_argument("--effort", default=None)
    ap.add_argument("--no-optimizer", action="store_true")
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

    arm = args.arm or tcfg.fit_default_arm
    target = tcfg.trial_target(**tcfg.fit_target) if tcfg.fit_target else {}

    sb = SimBiologyEngine()
    try:
        print("== starting MATLAB/SimBiology engine ==", flush=True)
        sb.start()
        print(f"== loading {os.path.basename(args.sbproj)} ==", flush=True)
        sb.load_project(args.sbproj)

        registry = ToolRegistry()
        register_qsp_fit_loop_tools(registry, cfg, {
            "cfg": tcfg, "sb": sb, "vpop": args.vpop, "arm": arm, "target": target,
            "limit": args.limit, "enable_optimize": not args.no_optimizer})

        goal = (f"Calibrate the model's PD parameter(s) on the fixed arm '{arm}' to "
                f"match the observed target {target}. Start with fit_inspect, then fit "
                "and finalize.")
        policy = LLMPolicy(cfg, registry, _system_prompt(tcfg))
        loop = DecisionLoop(config=cfg, registry=registry, policy=policy)

        print(f"\n== LLM calibration loop (target {target}) ==\n", flush=True)

        def show(ev):
            if isinstance(ev, Decision):
                if ev.text:
                    print(f"\n[reason] {ev.text[:2000]}")
                for c in ev.calls:
                    print(f"  -> {c.name} "
                          f"{json.dumps(c.arguments, ensure_ascii=False)[:300]}")
                    if c.name in ("fit_try", "fit_optimize"):
                        print("     ...simulating...", flush=True)
            elif isinstance(ev, Observation):
                print(f"  <- {ev.tool}: {ev.content.get('message', '')}", flush=True)
            elif isinstance(ev, Finish):
                print(f"\n=== AGENT SUMMARY ===\n{ev.text}")

        session = loop.run(goal, ModelingSession(goal=goal), on_event=show)

        final = session.get("fit_final")
        print("\n== SCORING ==")
        if not final:
            print("agent did not finalize a fit.")
        else:
            sc = final["score"]
            print(f"target (observed): {target}")
            print(f"committed fit: {final['overrides']}")
            print(f"predicted: {final['predicted']}")
            print(f"fit quality: MAE {sc.get('acr_mae_pp')} pp vs observed")
            for name, p in (sc.get("parameters") or {}).items():
                print(f"parameter recovery: {name} fitted {p.get('fitted')} vs "
                      f"reference {p.get('reference')} "
                      f"(log10-fold {p.get('log10_fold_from_ref')})")
    finally:
        sb.stop()


if __name__ == "__main__":
    main()
