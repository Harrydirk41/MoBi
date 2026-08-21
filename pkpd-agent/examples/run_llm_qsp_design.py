r"""LLM-driven QSP drug design: the agent designs a new pathway inhibitor.

Model-agnostic (--model loads projects/<name>/tasks.json). The agent chooses which
disease-driver pathway to inhibit and how hard; the harness edits the model to add
the drug, simulates, and reports the response. It screens pathways to discover the
best target and a sensible efficacy.

    set ANTHROPIC_API_KEY=...
    python -m examples.run_llm_qsp_design --model ra ^
        --sbproj "..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj" ^
        --vpop "..\RA-QSP-Model\Vpop1.xlsx" --limit 40 --max-steps 12
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
from pkpd_agent.tools.qsp_design_loop_tools import register_qsp_design_loop_tools


def _system_prompt(cfg) -> str:
    return (
        f"You are a QSP modeler designing a new drug on the {cfg.name}. You choose a "
        "disease-driver pathway to inhibit and how strongly (efficacy); the harness "
        "edits the model to add the drug and simulates the response. You have "
        "design_inspect (the targetable pathways with mechanism and real-world "
        "analogue), design_try (build + simulate a drug), and design_finalize "
        "(commit). Run efficacy 0 once for the background-alone baseline, then screen "
        "the pathways at a meaningful efficacy to find the best target, then tune. "
        "Reason about which pathway should dominate from the biology, and check "
        "whether the model agrees. Finalize the drug you recommend.")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="ra")
    ap.add_argument("--sbproj", required=True)
    ap.add_argument("--vpop", required=True)
    ap.add_argument("--limit", type=int, default=40)
    ap.add_argument("--max-steps", type=int, default=12)
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
        register_qsp_design_loop_tools(registry, cfg, {
            "cfg": tcfg, "sb": sb, "sbproj": args.sbproj, "vpop": args.vpop,
            "limit": args.limit})

        goal = ("Design the pathway inhibitor that gives the best response on the "
                "given background. Start with design_inspect, screen the pathways, "
                "then tune and finalize.")
        policy = LLMPolicy(cfg, registry, _system_prompt(tcfg))
        loop = DecisionLoop(config=cfg, registry=registry, policy=policy)

        print(f"\n== LLM drug-design loop ({args.vpop}, limit {args.limit}) ==\n",
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

        final = session.get("design_final")
        print("\n== RESULT ==")
        if not final:
            print("agent did not finalize a design; best:", session.get("design_best"))
        else:
            d = tcfg.design_targets.get(final["target"], {})
            print(f"designed drug: anti-{d.get('pathway', final['target'])} "
                  f"(target {final['target']}, efficacy {final['efficacy']})")
            print(f"real-world analogue: {d.get('analogue', 'n/a')}")
            print(f"response: {final['response']}")
    finally:
        sb.stop()


if __name__ == "__main__":
    main()
