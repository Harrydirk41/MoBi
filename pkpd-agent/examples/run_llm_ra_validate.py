r"""LLM-driven RA Stage-6: reproduce the paper's HELD-OUT validation.

The paper validated its model by predicting tocilizumab's response in a REFRACTORY
subpopulation (patients who inadequately respond to prior MTX and ADA) and comparing
to a real trial. That subpopulation is not a shipped model flag - it must be built
by running each prior therapy and intersecting the non-responders. This drives the
agent to DESIGN that two-stage inadequate-responder selection, then read TCZ's
response in it and compare to the real refractory-population trial.

    set ANTHROPIC_API_KEY=...
    python -m examples.run_llm_ra_validate ^
        --sbproj "..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj" ^
        --vpop   "..\RA-QSP-Model\Vpop1.xlsx" --limit 60 --max-steps 8

Each validate_run does THREE population simulations (MTX, ADA, TCZ), so it is the
slowest task - keep --limit modest and --max-steps low.
"""

from __future__ import annotations

import argparse
import json
import os

from pkpd_agent.config import AgentConfig
from pkpd_agent.engines.simbiology import SimBiologyEngine
from pkpd_agent.engines import ra_clinical_reference as ra_clin
from pkpd_agent.llm import LLMPolicy
from pkpd_agent.loop import DecisionLoop
from pkpd_agent.state import Decision, Finish, ModelingSession, Observation
from pkpd_agent.tools.registry import ToolRegistry
from pkpd_agent.tools.ra_validate_loop_tools import register_ra_validate_loop_tools


def _system_prompt(comparator: dict) -> str:
    return (
        "You are a QSP modeler running a model's HELD-OUT VALIDATION on the "
        "SimBiology engine. The validation predicts tocilizumab's response in a "
        "REFRACTORY population - patients who inadequately respond to prior therapies "
        "- and compares it to a real clinical trial. The refractory population is not "
        "given; you must construct it.\n\n"
        "1. Call validate_inspect: the goal, the available prior-therapy and TCZ "
        "arms, and the inadequate-responder convention.\n"
        "2. DESIGN the selection: refractoriness means failing prior therapy. To be a "
        "strong validation, define it as failing MORE than one line - run MTX and ADA "
        "as prior therapies and take patients who inadequately respond to BOTH (the "
        "intersection), using the clinical criterion ACR<50 with active disease "
        "DAS28-CRP>3.2.\n"
        "3. Call validate_run with the prior_arms and the TCZ arm; read how many "
        "patients each therapy fails, how many fail all (the refractory n), and TCZ's "
        "ACR20/50/70 in that subgroup vs the real trial.\n"
        "4. JUDGE the match: the paper reported a good match for ACR50/70 and "
        "remission with ACR20 harder to match. Compare your prediction to the real "
        f"trial ({comparator.get('trial')}: ACR20 {comparator.get('ACR20')}, "
        f"ACR50 {comparator.get('ACR50')}, ACR70 {comparator.get('ACR70')}). A small "
        "refractory n means a noisy estimate - note it.\n"
        "5. Call validate_finalize, then finish with: how you defined the refractory "
        "population, its size, the predicted vs real response, and whether the model "
        "validates."
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sbproj", required=True)
    ap.add_argument("--vpop", required=True)
    ap.add_argument("--limit", type=int, default=60)
    ap.add_argument("--max-steps", type=int, default=8)
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

    comparator = ra_clin.REFRACTORY_TCZ

    sb = SimBiologyEngine()
    try:
        print("== starting MATLAB/SimBiology engine ==", flush=True)
        sb.start()
        print(f"== loading {os.path.basename(args.sbproj)} ==", flush=True)
        sb.load_project(args.sbproj)

        registry = ToolRegistry()
        register_ra_validate_loop_tools(registry, cfg, {
            "sb": sb, "vpop": args.vpop, "limit": args.limit,
            "comparator": comparator})

        goal = ("Reproduce the model's held-out validation: predict TCZ's response in "
                "the refractory (multi-line inadequate responder) population and "
                "compare to the real trial. Start with validate_inspect.")
        policy = LLMPolicy(cfg, registry, _system_prompt(comparator))
        loop = DecisionLoop(config=cfg, registry=registry, policy=policy)

        print(f"\n== LLM RA held-out validation loop (limit {args.limit}) ==\n",
              flush=True)

        def show(ev):
            if isinstance(ev, Decision):
                if ev.text:
                    print(f"\n[reason] {ev.text[:2000]}")
                for c in ev.calls:
                    print(f"  -> {c.name} "
                          f"{json.dumps(c.arguments, ensure_ascii=False)[:300]}")
                    if c.name == "validate_run":
                        print("     ...simulating MTX, ADA, and TCZ arms...", flush=True)
            elif isinstance(ev, Observation):
                print(f"  <- {ev.tool}: {ev.content.get('message', '')}", flush=True)
            elif isinstance(ev, Finish):
                print(f"\n=== AGENT SUMMARY ===\n{ev.text}")

        session = loop.run(goal, ModelingSession(goal=goal), on_event=show)

        print("\n== SCORING [held-out validation] ==")
        final = session.get("val_final")
        if not final:
            print("[warn] agent did not validate_finalize.")
        else:
            print(f"real comparator: {comparator.get('trial')}")
            print(f"  ACR20 {comparator.get('ACR20')} ACR50 {comparator.get('ACR50')} "
                  f"ACR70 {comparator.get('ACR70')}")
            print(f"refractory population: n={final['n_refractory']} "
                  f"(failed {', '.join(final['prior_arms'])})")
            print(f"predicted TCZ response: {final['predicted']}")
            print(f"MAE vs real trial: {final['mae']} pp")
    finally:
        sb.stop()


if __name__ == "__main__":
    main()
