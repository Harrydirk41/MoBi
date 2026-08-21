r"""One general agent, all the downstream QSP modeling capabilities in one session.

Registers every task loop's tools into ONE registry and hands a single agent a
meta-goal: exercise the model across trial design, calibration, virtual-population
generation, drug design, and the held-out validation. The agent orchestrates -
picking tools across all five families in one continuous session.

Model-agnostic: everything specific to the model (drugs, parameters, targets, the
readout column roles) is loaded from projects/<name>/tasks.json via --model. Point
it at a different QSP model by adding a projects/<name>/ folder.

Honest scope: this is "one agent, many capabilities", NOT a from-scratch build
pipeline. The tasks are independent probes of the finished model (each starts from
the shipped, calibrated model); the drug-design task mutates the model but restores
it, so the tasks don't contaminate each other.

    set ANTHROPIC_API_KEY=...
    python -m examples.run_llm_qsp_full --model ra ^
        --sbproj "..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj" ^
        --vpop   "..\RA-QSP-Model\Vpop1.xlsx" --limit 40 --max-steps 30

This is LONG - five tasks x several population sims each. Keep --limit modest.
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
from pkpd_agent.tools.qsp_trial_loop_tools import register_qsp_trial_loop_tools
from pkpd_agent.tools.qsp_fit_loop_tools import register_qsp_fit_loop_tools
from pkpd_agent.tools.qsp_vpop_loop_tools import register_qsp_vpop_loop_tools
from pkpd_agent.tools.qsp_design_loop_tools import register_qsp_design_loop_tools
from pkpd_agent.tools.qsp_validate_loop_tools import register_qsp_validate_loop_tools


def _system_prompt(cfg) -> str:
    return (
        f"You are a QSP modeler with a single session on the validated {cfg.name} "
        "model (SimBiology). You have FIVE families of tools and should exercise each "
        "in turn, committing an answer for each before moving on:\n\n"
        "1. TRIAL DESIGN (trial_inspect / trial_run / trial_finalize): predict the "
        "response of first-line inadequate responders escalated to a second-line "
        "therapy; choose the drug and sequencing.\n"
        "2. CALIBRATION (fit_inspect / fit_optimize / fit_try / fit_finalize): fit a "
        "PD parameter so the model reproduces a real trial; set up the numerical fit "
        "and interpret it (do not hand-search).\n"
        "3. VPOP GENERATION (vpop_inspect / vpop_select / vpop_sample / "
        "vpop_finalize): design a sampling over disease drivers with WIDE bounds and "
        "let the numerical selector reweight to the target distribution.\n"
        "4. DRUG DESIGN (design_inspect / design_try / design_finalize): design a new "
        "pathway inhibitor - screen pathways, pick the best target and efficacy.\n"
        "5. VALIDATION (validate_inspect / validate_run / validate_finalize): "
        "reproduce the held-out validation - build the refractory (multi-line "
        "inadequate-responder) population and compare the test therapy's response to "
        "a real trial.\n\n"
        "Work efficiently - a couple of tool calls per task is enough. For the "
        "numerical tasks, use the optimizer/selector rather than hand-searching. After "
        "finalizing all five, finish with a concise cross-task summary: what you found "
        "in each, and one honest caveat about the model that showed up across tasks.")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="ra", help="project name (projects/<name>/tasks.json)")
    ap.add_argument("--sbproj", required=True)
    ap.add_argument("--vpop", required=True)
    ap.add_argument("--limit", type=int, default=40)
    ap.add_argument("--max-steps", type=int, default=30)
    ap.add_argument("--llm-model", dest="llm_model", default=None)
    ap.add_argument("--effort", default=None)
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

    sb = SimBiologyEngine()
    try:
        print("== starting MATLAB/SimBiology engine ==", flush=True)
        sb.start()
        print(f"== loading {os.path.basename(args.sbproj)} ==", flush=True)
        sb.load_project(args.sbproj)

        registry = ToolRegistry()
        register_qsp_trial_loop_tools(registry, cfg, {
            "cfg": tcfg, "sb": sb, "vpop": args.vpop, "limit": args.limit,
            "stop_time": 700.0})
        register_qsp_fit_loop_tools(registry, cfg, {
            "cfg": tcfg, "sb": sb, "vpop": args.vpop, "arm": tcfg.fit_default_arm,
            "target": tcfg.trial_target(**tcfg.fit_target) if tcfg.fit_target else {},
            "limit": args.limit, "stop_time": 700.0, "enable_optimize": True})
        register_qsp_vpop_loop_tools(registry, cfg, {
            "cfg": tcfg, "sb": sb, "n_samples": max(args.limit, 60), "seed": 1,
            "enable_select": True, "n_pool": 80})
        register_qsp_design_loop_tools(registry, cfg, {
            "cfg": tcfg, "sb": sb, "sbproj": args.sbproj, "vpop": args.vpop,
            "limit": args.limit})
        register_qsp_validate_loop_tools(registry, cfg, {
            "cfg": tcfg, "sb": sb, "vpop": args.vpop, "limit": max(args.limit, 60)})

        print(f"   registered {len(registry)} tools across 5 task families", flush=True)

        goal = (f"Exercise the {tcfg.name} across all five capabilities - trial "
                "design, calibration, Vpop generation, drug design, and the held-out "
                "validation - committing an answer for each, then summarize.")
        policy = LLMPolicy(cfg, registry, _system_prompt(tcfg))
        loop = DecisionLoop(config=cfg, registry=registry, policy=policy)

        print(f"\n== ONE-AGENT session (limit {args.limit}, "
              f"max {args.max_steps} steps) ==\n", flush=True)

        slow = {"trial_run", "fit_try", "fit_optimize", "vpop_sample",
                "vpop_select", "design_try", "validate_run"}

        def show(ev):
            if isinstance(ev, Decision):
                if ev.text:
                    print(f"\n[reason] {ev.text[:1500]}")
                for c in ev.calls:
                    print(f"  -> {c.name} "
                          f"{json.dumps(c.arguments, ensure_ascii=False)[:220]}")
                    if c.name in slow:
                        print("     ...simulating...", flush=True)
            elif isinstance(ev, Observation):
                print(f"  <- {ev.tool}: {ev.content.get('message', '')}", flush=True)
            elif isinstance(ev, Finish):
                print(f"\n=== AGENT SUMMARY ===\n{ev.text}")

        session = loop.run(goal, ModelingSession(goal=goal), on_event=show)

        print("\n== COMMITTED ANSWERS ==")
        labels = [("trial", "trial_final"), ("calibration", "fit_final"),
                  ("vpop", "vpop_final"), ("drug design", "design_final"),
                  ("validation", "val_final")]
        for name, key in labels:
            v = session.get(key)
            print(f"  {name:12}: {'committed' if v else 'not finalized'}")
    finally:
        sb.stop()


if __name__ == "__main__":
    main()
