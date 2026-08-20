r"""One general agent, all the RA modeling capabilities in a single session.

Registers every task loop's tools into ONE registry and hands a single agent a
meta-goal: exercise the model across trial design, calibration, virtual-population
generation, drug design, and the held-out validation. The agent orchestrates -
picking tools across all five families in one continuous session - rather than five
separate scripts.

Honest scope: this is "one agent, many capabilities", NOT a from-scratch build
pipeline. The tasks are independent probes of the finished model (each starts from
the shipped, calibrated model); the drug-design task mutates the model but restores
it, so the tasks don't contaminate each other. It shows an agent can carry all the
decisions around a validated QSP model in one loop - not that it builds the model.

    set ANTHROPIC_API_KEY=...
    python -m examples.run_llm_ra_full ^
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
from pkpd_agent.engines import osp_ra_trial
from pkpd_agent.engines import ra_clinical_reference as ra_clin
from pkpd_agent.llm import LLMPolicy
from pkpd_agent.loop import DecisionLoop
from pkpd_agent.state import Decision, Finish, ModelingSession, Observation
from pkpd_agent.tools.registry import ToolRegistry
from pkpd_agent.tools.ra_trial_loop_tools import register_ra_trial_loop_tools
from pkpd_agent.tools.ra_fit_loop_tools import register_ra_fit_loop_tools
from pkpd_agent.tools.ra_vpop_loop_tools import register_ra_vpop_loop_tools
from pkpd_agent.tools.ra_design_loop_tools import register_ra_design_loop_tools
from pkpd_agent.tools.ra_validate_loop_tools import register_ra_validate_loop_tools

_FLAGSHIP = "MTX_15mg_Q1W_SC_t200;TCZ8mgkg_Q4W_IV_t200@285"


def _system_prompt() -> str:
    return (
        "You are a QSP modeler with a single session on a validated rheumatoid-"
        "arthritis model (SimBiology). You have FIVE families of tools and should "
        "exercise each in turn, committing an answer for each before moving on:\n\n"
        "1. TRIAL DESIGN (ra_inspect / ra_run_trial / ra_finalize): predict the "
        "response of MTX-inadequate responders escalated to a second-line therapy; "
        "choose the drug and sequencing.\n"
        "2. CALIBRATION (fit_inspect / fit_optimize / fit_try / fit_finalize): fit "
        "the PD parameter KD_TCZ so the model reproduces a real trial; set up the "
        "numerical fit and interpret it (do not hand-search).\n"
        "3. VPOP GENERATION (vpop_inspect / vpop_select / vpop_sample / "
        "vpop_finalize): design a sampling over disease drivers with WIDE bounds and "
        "let the numerical selector reweight to the target distribution.\n"
        "4. DRUG DESIGN (design_inspect / design_try / design_finalize): design a "
        "new anti-cytokine biologic - screen pathways, pick the best target and a "
        "sensible efficacy.\n"
        "5. VALIDATION (validate_inspect / validate_run / validate_finalize): "
        "reproduce the held-out validation - build the refractory (multi-line "
        "inadequate-responder) population and compare TCZ's response to a real "
        "trial.\n\n"
        "Work efficiently - a couple of tool calls per task is enough; you do not "
        "need to exhaustively optimize each one. For the numerical tasks, use the "
        "optimizer/selector rather than hand-searching. After finalizing all five, "
        "finish with a concise cross-task summary: what you found in each, and one "
        "honest caveat about the model that showed up across tasks."
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sbproj", required=True)
    ap.add_argument("--vpop", required=True)
    ap.add_argument("--limit", type=int, default=40)
    ap.add_argument("--max-steps", type=int, default=30)
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

    sb = SimBiologyEngine()
    try:
        print("== starting MATLAB/SimBiology engine ==", flush=True)
        sb.start()
        print(f"== loading {os.path.basename(args.sbproj)} ==", flush=True)
        sb.load_project(args.sbproj)

        registry = ToolRegistry()
        # 1. trial design
        register_ra_trial_loop_tools(registry, cfg, {
            "sb": sb, "vpop": args.vpop, "limit": args.limit, "stop_time": 700.0,
            "objective": "Predict the second-line response in MTX-inadequate responders.",
            "calibrated_arms": [{
                "arm": "MTX first-line monotherapy",
                "protocol": {"first_line": ["MTX_15mg_Q1W_SC_t200"]},
                "known_rates_day284_full_pop_n300": {"ACR20": 33.7, "ACR50": 17.7,
                                                     "ACR70": 2.3, "remission": 18.0}}]})
        # 2. calibration
        register_ra_fit_loop_tools(registry, cfg, {
            "sb": sb, "vpop": args.vpop, "arm": _FLAGSHIP,
            "target": ra_clin.trial_target("TCZ", 24, "raw"),
            "fit_params": ["KD_TCZ"], "limit": args.limit, "stop_time": 700.0,
            "enable_optimize": True})
        # 3. vpop generation
        register_ra_vpop_loop_tools(registry, cfg, {
            "sb": sb, "n_samples": max(args.limit, 60), "baseline_day": 200.0,
            "seed": 1, "target": osp_ra_trial.VPOP_TARGET,
            "enable_select": True, "n_pool": 80})
        # 4. drug design
        register_ra_design_loop_tools(registry, cfg, {
            "sb": sb, "sbproj": args.sbproj, "vpop": args.vpop,
            "background": "MTX_15mg_Q1W_SC_t200", "start_day": 200.0,
            "readout_day": 284.0, "limit": args.limit})
        # 5. validation
        register_ra_validate_loop_tools(registry, cfg, {
            "sb": sb, "vpop": args.vpop, "limit": max(args.limit, 60),
            "comparator": ra_clin.REFRACTORY_TCZ})

        print(f"   registered {len(registry)} tools across 5 task families",
              flush=True)

        goal = ("Exercise the RA model across all five capabilities - trial design, "
                "calibration, Vpop generation, drug design, and the held-out "
                "validation - committing an answer for each, then summarize.")
        policy = LLMPolicy(cfg, registry, _system_prompt())
        loop = DecisionLoop(config=cfg, registry=registry, policy=policy)

        print(f"\n== ONE-AGENT RA session (limit {args.limit}, "
              f"max {args.max_steps} steps) ==\n", flush=True)

        slow = {"ra_run_trial", "fit_try", "fit_optimize", "vpop_sample",
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
        labels = [("trial", "ra_final"), ("calibration", "fit_final"),
                  ("vpop", "vpop_final"), ("drug design", "design_final"),
                  ("validation", "val_final")]
        for name, key in labels:
            v = session.get(key)
            print(f"  {name:12}: {'committed' if v else 'not finalized'}")
    finally:
        sb.stop()


if __name__ == "__main__":
    main()
