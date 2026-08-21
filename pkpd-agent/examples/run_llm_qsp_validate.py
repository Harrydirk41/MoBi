r"""LLM-driven QSP held-out validation: reproduce the paper's refractory prediction.

Model-agnostic (--model loads projects/<name>/tasks.json). The agent constructs the
refractory subpopulation (run each prior therapy, intersect the inadequate
responders), gives the test therapy, and compares its response in that subgroup to a
real trial.

    set ANTHROPIC_API_KEY=...
    python -m examples.run_llm_qsp_validate --model ra ^
        --sbproj "..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj" ^
        --vpop "..\RA-QSP-Model\Vpop1.xlsx" --limit 60 --max-steps 10
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
from pkpd_agent.tools.qsp_validate_loop_tools import register_qsp_validate_loop_tools


def _system_prompt(cfg, comparator) -> str:
    return (
        f"You are a QSP modeler reproducing the held-out validation of the {cfg.name}. "
        "The refractory population is not a shipped flag - you construct it by running "
        "the prior therapies and intersecting their inadequate responders, then give "
        "the test therapy and read its response in that subgroup. You have "
        "validate_inspect (goal, available arms, IR convention), validate_run (run the "
        "selection and score against the real comparator), and validate_finalize "
        f"(commit). The real comparator is {comparator.get('trial')}. Design the "
        "selection - which prior therapies define refractoriness, and the IR criteria "
        "- and finalize the design you recommend.")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="vantage_ra")
    ap.add_argument("--sbproj", required=True)
    ap.add_argument("--vpop", required=True)
    ap.add_argument("--limit", type=int, default=60)
    ap.add_argument("--max-steps", type=int, default=10)
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

    comparator = tcfg.refractory_target

    sb = SimBiologyEngine()
    try:
        print("== starting MATLAB/SimBiology engine ==", flush=True)
        sb.start()
        print(f"== loading {os.path.basename(args.sbproj)} ==", flush=True)
        sb.load_project(args.sbproj)

        registry = ToolRegistry()
        register_qsp_validate_loop_tools(registry, cfg, {
            "cfg": tcfg, "sb": sb, "vpop": args.vpop, "limit": args.limit,
            "comparator": comparator})

        goal = ("Reproduce the held-out validation: build the refractory population, "
                "give the test therapy, and compare to the real comparator. Start with "
                "validate_inspect, then run and finalize.")
        policy = LLMPolicy(cfg, registry, _system_prompt(tcfg, comparator))
        loop = DecisionLoop(config=cfg, registry=registry, policy=policy)

        print(f"\n== LLM held-out validation loop (limit {args.limit}) ==\n", flush=True)

        def show(ev):
            if isinstance(ev, Decision):
                if ev.text:
                    print(f"\n[reason] {ev.text[:2000]}")
                for c in ev.calls:
                    print(f"  -> {c.name} "
                          f"{json.dumps(c.arguments, ensure_ascii=False)[:300]}")
                    if c.name == "validate_run":
                        print("     ...running arms + classifying...", flush=True)
            elif isinstance(ev, Observation):
                print(f"  <- {ev.tool}: {ev.content.get('message', '')}", flush=True)
            elif isinstance(ev, Finish):
                print(f"\n=== AGENT SUMMARY ===\n{ev.text}")

        session = loop.run(goal, ModelingSession(goal=goal), on_event=show)

        final = session.get("val_final")
        print("\n== SCORING ==")
        if not final:
            print("agent did not finalize a validation.")
        else:
            print(f"real comparator: {comparator.get('trial')}")
            print(f"refractory population: n={final['n_refractory']}")
            print(f"predicted test response: {final['predicted']}")
            print(f"MAE vs real trial: {final['mae']} pp")
    finally:
        sb.stop()


if __name__ == "__main__":
    main()
