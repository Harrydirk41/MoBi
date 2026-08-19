r"""LLM-driven RA Stage-2: the agent DESIGNS a new anti-cytokine biologic.

The most open-ended task: the agent designs a drug from scratch - it picks which
cytokine PATHWAY to inhibit and how strongly - and the harness edits the model
structurally (adds the drug as a suppression of that pathway's driver), simulates,
and reports the ACR. The agent screens pathways and tunes the drug, discovering
which target makes the best RA drug (it should recover that IL-6 / TNF dominate,
matching real biology, and that IL-17 underperforms as it does clinically).

    set ANTHROPIC_API_KEY=...
    python -m examples.run_llm_ra_design ^
        --sbproj "..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj" ^
        --vpop   "..\RA-QSP-Model\Vpop1.xlsx" --limit 40 --max-steps 12

Each design_try reloads the model + simulates the cohort, so keep --limit modest.
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
from pkpd_agent.tools.ra_design_loop_tools import register_ra_design_loop_tools


def _system_prompt() -> str:
    return (
        "You are a translational pharmacologist DESIGNING a new biologic for "
        "rheumatoid arthritis on a QSP model with the SimBiology engine. You design "
        "the drug: choose which cytokine PATHWAY it inhibits and how strongly "
        "(efficacy). The model is edited to add your drug and the clinical response "
        "is simulated on a methotrexate background.\n\n"
        "Work like a drug hunter:\n"
        "1. Call design_inspect: the targetable pathways (each with mechanism and "
        "real-world analogue) and the ACR readout.\n"
        "2. Establish the MTX-alone baseline (design_try with efficacy 0) so you can "
        "measure each drug's ADDED benefit.\n"
        "3. SCREEN the pathways: run each candidate target at a meaningful efficacy "
        "(e.g. 0.8) and compare the ACR. Reason from mechanism about which cytokine "
        "should dominate a DAS28-CRP readout, but let the model rank them.\n"
        "4. For the best pathway, TUNE efficacy - find where the benefit saturates, "
        "and report the smallest efficacy that reaches the plateau (a lower "
        "efficacious dose is a better, safer drug).\n"
        "5. SANITY-CHECK against biology: does the ranking match what is known in RA "
        "(IL-6 and TNF effective, IL-17 weak)? Flag any pathway the model rates very "
        "differently from clinical reality.\n"
        "6. Call design_finalize with your chosen {target, efficacy}, then finish "
        "with: the drug you designed, its mechanism and real-world analogue, the "
        "response, and how the model's pathway ranking compares to clinical reality."
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sbproj", required=True)
    ap.add_argument("--vpop", required=True)
    ap.add_argument("--background", default="MTX_15mg_Q1W_SC_t200",
                    help="background therapy dose (blank for monotherapy)")
    ap.add_argument("--start-day", type=float, default=200.0)
    ap.add_argument("--readout-day", type=float, default=284.0)
    ap.add_argument("--limit", type=int, default=40)
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

    sb = SimBiologyEngine()
    try:
        print("== starting MATLAB/SimBiology engine ==", flush=True)
        sb.start()
        print(f"== loading {os.path.basename(args.sbproj)} ==", flush=True)
        sb.load_project(args.sbproj)

        registry = ToolRegistry()
        register_ra_design_loop_tools(registry, cfg, {
            "sb": sb, "sbproj": args.sbproj, "vpop": args.vpop,
            "background": args.background, "start_day": args.start_day,
            "readout_day": args.readout_day, "limit": args.limit})

        goal = ("Design a new anti-cytokine biologic for RA that maximizes the ACR "
                "response on an MTX background. Start with design_inspect, baseline "
                "MTX alone, screen the pathways, tune, and finalize.")
        policy = LLMPolicy(cfg, registry, _system_prompt())
        loop = DecisionLoop(config=cfg, registry=registry, policy=policy)

        print(f"\n== LLM RA drug-design loop ({args.vpop}, limit {args.limit}) ==\n",
              flush=True)

        def show(ev):
            if isinstance(ev, Decision):
                if ev.text:
                    print(f"\n[reason] {ev.text[:2000]}")
                for c in ev.calls:
                    print(f"  -> {c.name} "
                          f"{json.dumps(c.arguments, ensure_ascii=False)[:300]}")
                    if c.name == "design_try":
                        print("     ...editing model + simulating...", flush=True)
            elif isinstance(ev, Observation):
                print(f"  <- {ev.tool}: {ev.content.get('message', '')}", flush=True)
            elif isinstance(ev, Finish):
                print(f"\n=== AGENT SUMMARY ===\n{ev.text}")

        session = loop.run(goal, ModelingSession(goal=goal), on_event=show)

        print("\n== RESULT [drug-design] ==")
        final = session.get("design_final")
        if not final:
            print("[warn] agent did not design_finalize; best was",
                  session.get("design_best"))
        else:
            d = osp_ra_trial.DESIGN_TARGETS.get(final["target"], {})
            fl = final["response"]
            print(f"designed drug: anti-{d.get('pathway', final['target'])} "
                  f"(target {final['target']}, efficacy {final['efficacy']})")
            print(f"real-world analogue: {d.get('analogue', 'n/a')}")
            print(f"response: ACR20 {fl.get('ACR20')}% ACR50 {fl.get('ACR50')}% "
                  f"ACR70 {fl.get('ACR70')}% remission {fl.get('remission')}%")
    finally:
        sb.stop()


if __name__ == "__main__":
    main()
