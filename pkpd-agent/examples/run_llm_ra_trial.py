r"""LLM-driven RA virtual-trial: the agent designs the protocol, SimBiology runs it.

The SimBiology counterpart of run_llm_ddi. Where the DDI loop tunes an interaction,
this drives the agent over an in-silico TRIAL of the Vantage RA QSP model: it reads
the disease/readout, the drug formulary, and the calibrated reference arms; reasons
about which therapy to give patients who fail first-line MTX and how to give it
(a mechanistically distinct biologic, switched in after the first-line readout);
runs the protocol across the virtual population; and predicts the second-line
response. The disease model and virtual population are fixed; the unknown is the
PROTOCOL.

The held-out truth (TCZ-in-MTX-inadequate-responders ACR20/50/70) is NOT shown to
the agent - it must pick the therapy from mechanism, and this script scores the
final prediction against the held-out numbers after the loop.

    set ANTHROPIC_API_KEY=...
    python -m examples.run_llm_ra_trial ^
        --sbproj "..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj" ^
        --vpop   "..\RA-QSP-Model\Vpop1.xlsx" ^
        --limit 50 --max-steps 8

Note: each ra_run_trial simulates the whole (subsampled) population, so keep
--limit modest (50) during the loop; confirm the winner at --limit 300 with
examples.run_ra_vpop.
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
from pkpd_agent.tools.ra_trial_loop_tools import register_ra_trial_loop_tools


def _system_prompt() -> str:
    return (
        "You are a clinical pharmacologist running an in-silico trial on a "
        "rheumatoid-arthritis QSP model with the SimBiology engine. The disease "
        "model and the virtual population are GIVEN and fixed; you design the "
        "treatment PROTOCOL.\n\n"
        "The task: predict the response of patients who INADEQUATELY RESPOND to "
        "first-line methotrexate (MTX) when they are escalated to a second-line "
        "therapy. Each round:\n"
        "1. Call ra_inspect: the DAS28/ACR readout, the trial timeline (baseline "
        "day 199, first-line readout day 284, second-line readout day 600), the "
        "drug formulary (MTX plus the biologics, each with its MECHANISM), and the "
        "calibrated reference arms with their known rates.\n"
        "2. VALIDATE your harness: run first-line MTX (ra_run_trial with only "
        "first_line) and confirm the ACR20/remission match the reference arm. This "
        "proves the population and readout are wired correctly.\n"
        "3. REASON about the second line. MTX-inadequate responders need a therapy "
        "that hits a DIFFERENT node of the inflammatory network than MTX - a "
        "biologic. Pick the one whose mechanism best matches an escalation for "
        "MTX-IR RA, and give it as a SWITCH: set second_line and a switch_day just "
        "after the day-284 readout (e.g. 285), NOT concurrently from day 200 (which "
        "would conflate the arms and give a dead second-line, ACR50/70 near 0).\n"
        "4. Call ra_run_trial with first_line MTX + your second_line + switch_day, "
        "and read the second-line ACR20/50/70 among the MTX non-responders.\n"
        "5. SANITY-CHECK: are the second-line rates clinically plausible for a "
        "biologic in MTX-IR RA (ACR20 roughly one-half, ACR50 roughly one-quarter, "
        "ACR70 lower)? If the second-line arm is empty or dead, fix the protocol "
        "(give MTX first-line, switch the biologic in later) and re-run.\n"
        "ENGINE/EMPTY-ARM RESULTS ARE FAILURES: if a run warns the second-line arm "
        "is empty, the protocol is wrong - fix it and re-run; do not finish on it.\n"
        "6. COMMIT: call ra_finalize with the protocol you actually recommend. This "
        "is the run that gets scored - if you ran controls or dose-response checks "
        "afterwards, ra_finalize makes sure your ANSWER is scored, not the last "
        "thing you happened to run. Finalize the therapy+dose you would give a "
        "patient, then finish.\n"
        "Finish with: the therapy you chose and WHY (mechanism), the protocol "
        "(doses + switch day), and your predicted second-line ACR20/50/70."
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sbproj", required=True)
    ap.add_argument("--vpop", required=True)
    ap.add_argument("--limit", type=int, default=50,
                    help="patients per trial run during the loop (default 50)")
    ap.add_argument("--stop-time", type=float, default=700.0)
    ap.add_argument("--max-steps", type=int, default=8)
    ap.add_argument("--model", default=None)
    ap.add_argument("--effort", default=None)
    # held-out target (the paper's validation; default = the validated model output)
    ap.add_argument("--target-acr20", type=float, default=44.9)
    ap.add_argument("--target-acr50", type=float, default=23.5)
    ap.add_argument("--target-acr70", type=float, default=14.0)
    args = ap.parse_args()

    cfg = AgentConfig(mock=False, max_steps=args.max_steps)
    if args.model:
        cfg.model = args.model
    if args.effort:
        cfg.effort = args.effort
    if not cfg.anthropic_key_present():
        print("ANTHROPIC_API_KEY not set.")
        return

    target = {"ACR20": args.target_acr20, "ACR50": args.target_acr50,
              "ACR70": args.target_acr70}

    # the calibrated first-line reference the agent validates its harness against.
    # These are FULL-POPULATION (n=300) rates; a subsampled run (--limit) will
    # differ by sampling noise (at n=50 the standard error on a ~40% rate is ~7pp
    # and the resolution is 2pp), so treat a match as "same ballpark", not exact.
    calibrated_arms = [{
        "arm": "MTX first-line monotherapy",
        "protocol": {"first_line": ["MTX_15mg_Q1W_SC_t200"]},
        "known_rates_day284_full_pop_n300": {"ACR20": 33.7, "ACR50": 17.7,
                                             "ACR70": 2.3, "remission": 18.0},
        "note": "run this first; expect these +/- sampling noise at your --limit",
    }]

    sb = SimBiologyEngine()
    try:
        print("== starting MATLAB/SimBiology engine ==", flush=True)
        sb.start()
        print(f"== loading {os.path.basename(args.sbproj)} ==", flush=True)
        sb.load_project(args.sbproj)

        registry = ToolRegistry()
        register_ra_trial_loop_tools(registry, cfg, {
            "sb": sb, "vpop": args.vpop, "calibrated_arms": calibrated_arms,
            "objective": "Predict the ACR20/50/70 response of MTX-inadequate "
                         "responders escalated to a second-line therapy.",
            "limit": args.limit, "stop_time": args.stop_time})

        goal = ("Predict the second-line response in MTX-inadequate responders. "
                "Start by calling ra_inspect, validate the harness on first-line "
                "MTX, then choose and run the second-line protocol.")
        policy = LLMPolicy(cfg, registry, _system_prompt())
        loop = DecisionLoop(config=cfg, registry=registry, policy=policy)

        print(f"\n== LLM RA virtual-trial loop ({args.vpop}, "
              f"limit {args.limit}/run) ==\n", flush=True)

        def show(ev):
            if isinstance(ev, Decision):
                if ev.text:
                    print(f"\n[reason] {ev.text[:2000]}")
                for c in ev.calls:
                    print(f"  -> {c.name} "
                          f"{json.dumps(c.arguments, ensure_ascii=False)[:400]}")
                    if c.name == "ra_run_trial":
                        print("     ...simulating the population...", flush=True)
            elif isinstance(ev, Observation):
                print(f"  <- {ev.tool}: {ev.content.get('message', '')}", flush=True)
            elif isinstance(ev, Finish):
                print(f"\n=== AGENT SUMMARY ===\n{ev.text}")

        session = loop.run(goal, ModelingSession(goal=goal), on_event=show)

        # score the agent's COMMITTED protocol (ra_finalize); fall back to the last
        # non-empty run only if it never finalized (and say so).
        final = session.get("ra_final")
        committed = final is not None
        if final is None:
            hist = session.get("ra_history") or []
            final = next((h for h in reversed(hist)
                          if (h.get("second_line") or {}).get("n_MTX_IR", 0) > 0), None)
        print("\n== SCORING vs held-out truth ==")
        print(f"held-out target (second-line, MTX-IR): {target}")
        if not committed:
            print("[warn] agent did not call ra_finalize - scoring its last "
                  "non-empty run, which may be a control rather than its answer.")
        if not final:
            print("no scorable protocol run (the agent never produced a non-empty "
                  "second-line arm).")
        else:
            pred = final["second_line"]
            score = osp_ra_trial.score_flagship(pred, target)
            print(f"agent protocol: {final['protocol']}")
            print(f"agent prediction (MTX-IR n={pred.get('n_MTX_IR')}): "
                  f"ACR20 {pred.get('ACR20')}% ACR50 {pred.get('ACR50')}% "
                  f"ACR70 {pred.get('ACR70')}%")
            print(f"score: MAE {score.get('mae_pp')} percentage points over "
                  f"{score.get('n_endpoints')} endpoints")
            for k, v in (score.get("per_endpoint") or {}).items():
                print(f"   {k}: predicted {v['predicted']} vs target {v['target']} "
                      f"(|err| {v['abs_error_pp']} pp)")
    finally:
        sb.stop()


if __name__ == "__main__":
    main()
