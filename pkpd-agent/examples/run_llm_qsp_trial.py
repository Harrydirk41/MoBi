r"""LLM-driven QSP virtual-trial: the agent designs the protocol, SimBiology runs it.

Model-agnostic (--model loads projects/<name>/tasks.json). The agent reads the
disease/readout, the drug formulary, and the calibrated reference arms; reasons about
which therapy to escalate first-line inadequate responders to; runs the protocol
across the virtual population; and predicts the second-line response. The disease
model and virtual population are fixed; the unknown is the PROTOCOL.

Two objectives:
  * --objective predict  (default): reproduce the held-out second-line rates.
  * --objective min-dose: find the LOWEST second-line dose clearing an endpoint bar.

    set ANTHROPIC_API_KEY=...
    python -m examples.run_llm_qsp_trial --model ra ^
        --sbproj "..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj" ^
        --vpop   "..\RA-QSP-Model\Vpop1.xlsx" --limit 50 --max-steps 10
"""

from __future__ import annotations

import argparse
import json
import os
import re

from pkpd_agent.config import AgentConfig
from pkpd_agent.engines.simbiology import SimBiologyEngine
from pkpd_agent.engines import qsp_config, qsp_tasks
from pkpd_agent.llm import LLMPolicy
from pkpd_agent.loop import DecisionLoop
from pkpd_agent.state import Decision, Finish, ModelingSession, Observation
from pkpd_agent.tools.registry import ToolRegistry
from pkpd_agent.tools.qsp_trial_loop_tools import register_qsp_trial_loop_tools


def _system_prompt(cfg, objective: str, thr: float) -> str:
    common = (
        f"You are a clinical pharmacologist running an in-silico trial on the "
        f"{cfg.name} ({cfg.disease}) with the SimBiology engine. The disease model "
        "and the virtual population are GIVEN and fixed; you design the treatment "
        "PROTOCOL. You have trial_inspect (disease/readout, timeline, drug formulary, "
        "calibrated reference), trial_run (run a protocol and read the response), and "
        "trial_finalize (commit the protocol you recommend; that is what gets "
        "scored).\n\nWork like a scientist, not a search: validate your harness "
        "against the calibrated reference first, then reason from MECHANISM about the "
        "escalation. Be suspicious of clinically implausible results or artifacts of "
        "how you set the trial up - diagnose rather than report them. A run with an "
        "empty subgroup arm is a mis-specified protocol, not a result; fix it. "
        "Finalize the protocol you would actually give a patient.\n\n")
    goals = {
        "predict": ("OBJECTIVE: predict the response of first-line inadequate "
                    "responders escalated to a second-line therapy. Choose the "
                    "therapy and protocol, run it, and commit your predicted "
                    "second-line rates."),
        "min-dose": (f"OBJECTIVE: find the LOWEST-DOSE second-line regimen still "
                     f"achieving the primary endpoint >= {thr:g}% in the inadequate-"
                     "responder subgroup. Use dose_scale to titrate. Maxing the dose "
                     "wastes drug; too little misses. Find the smallest passing "
                     "dose_scale, then commit that protocol."),
    }
    return common + goals[objective]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="vantage_ra")
    ap.add_argument("--sbproj", required=True)
    ap.add_argument("--vpop", required=True)
    ap.add_argument("--limit", type=int, default=50)
    ap.add_argument("--stop-time", type=float, default=700.0)
    ap.add_argument("--max-steps", type=int, default=8)
    ap.add_argument("--llm-model", dest="llm_model", default=None)
    ap.add_argument("--effort", default=None)
    ap.add_argument("--objective", choices=["predict", "min-dose"], default="predict")
    ap.add_argument("--min-dose-target", type=float, default=35.0,
                    help="min-dose objective: the primary-endpoint %% bar to clear")
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

    # scoring target for the predict objective: the model's clinical reference
    target = tcfg.trial_target(**tcfg.fit_target) if tcfg.fit_target else {}
    primary = next(iter(tcfg.run_columns.get("first_line", {"ACR20": 1})), "ACR20")

    sb = SimBiologyEngine()
    try:
        print("== starting MATLAB/SimBiology engine ==", flush=True)
        sb.start()
        print(f"== loading {os.path.basename(args.sbproj)} ==", flush=True)
        sb.load_project(args.sbproj)

        registry = ToolRegistry()
        register_qsp_trial_loop_tools(registry, cfg, {
            "cfg": tcfg, "sb": sb, "vpop": args.vpop, "limit": args.limit,
            "stop_time": args.stop_time})

        goal = (tcfg.trial_objective + " Start by calling trial_inspect and validating "
                "the harness on the calibrated arm, then design and run the second-line "
                "protocol.")
        policy = LLMPolicy(cfg, registry,
                           _system_prompt(tcfg, args.objective, args.min_dose_target))
        loop = DecisionLoop(config=cfg, registry=registry, policy=policy)

        print(f"\n== LLM virtual-trial loop [{args.objective}] ({args.vpop}, "
              f"limit {args.limit}/run) ==\n", flush=True)

        def show(ev):
            if isinstance(ev, Decision):
                if ev.text:
                    print(f"\n[reason] {ev.text[:2000]}")
                for c in ev.calls:
                    print(f"  -> {c.name} "
                          f"{json.dumps(c.arguments, ensure_ascii=False)[:400]}")
                    if c.name == "trial_run":
                        print("     ...simulating the population...", flush=True)
            elif isinstance(ev, Observation):
                print(f"  <- {ev.tool}: {ev.content.get('message', '')}", flush=True)
            elif isinstance(ev, Finish):
                print(f"\n=== AGENT SUMMARY ===\n{ev.text}")

        session = loop.run(goal, ModelingSession(goal=goal), on_event=show)

        final = session.get("trial_final")
        committed = final is not None
        if final is None:
            hist = session.get("trial_history") or []
            final = next((h for h in reversed(hist)
                          if (h.get("second_line") or {}).get("n_subgroup", 0) > 0), None)
        print(f"\n== SCORING [{args.objective}] ==")
        if not committed:
            print("[warn] agent did not call trial_finalize - scoring its last "
                  "non-empty run.")
        if not final:
            print("no scorable protocol run.")
        elif args.objective == "predict":
            pred = final["second_line"]
            print(f"agent protocol: {final['protocol']}")
            print(f"agent prediction (subgroup n={pred.get('n_subgroup')}): {pred}")
            if target:
                score = qsp_tasks.score_flagship(pred, target)
                print(f"\nscore vs clinical reference {target}: MAE "
                      f"{score.get('mae_pp')} pp over {score.get('n_endpoints')} endpoints")
        else:
            pred = final["second_line"]
            m = re.search(r"\*([0-9.]+)", final["protocol"])
            scale = float(m.group(1)) if m else 1.0
            score = qsp_tasks.score_min_dose(pred, scale, args.min_dose_target, primary)
            print(f"agent protocol: {final['protocol']}  (dose_scale={scale:g})")
            print(f"achieved: {primary} {pred.get(primary)}% -> {score['verdict']}")
    finally:
        sb.stop()


if __name__ == "__main__":
    main()
