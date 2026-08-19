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

Two objectives:
  * --objective predict  (default): reproduce the held-out TCZ-in-MTX-IR rates.
  * --objective min-dose: find the LOWEST second-line dose that clears an ACR20
    bar (--min-dose-acr20). This one is not brute-forceable by "pick the max
    response" - the agent must titrate dose_scale to the efficient point.

The tool descriptions no longer hint the drug, the switch timing or the dose -
the agent derives them from mechanism and experiment.

    set ANTHROPIC_API_KEY=...
    python -m examples.run_llm_ra_trial ^
        --sbproj "..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj" ^
        --vpop   "..\RA-QSP-Model\Vpop1.xlsx" ^
        --limit 50 --max-steps 10 --objective min-dose

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
from pkpd_agent.engines import ra_clinical_reference as ra_clin
from pkpd_agent.llm import LLMPolicy
from pkpd_agent.loop import DecisionLoop
from pkpd_agent.state import Decision, Finish, ModelingSession, Observation
from pkpd_agent.tools.registry import ToolRegistry
from pkpd_agent.tools.ra_trial_loop_tools import register_ra_trial_loop_tools


_COMMON = (
    "You are a clinical pharmacologist running an in-silico trial on a "
    "rheumatoid-arthritis QSP model with the SimBiology engine. The disease model "
    "and the virtual population are GIVEN and fixed; you design the treatment "
    "PROTOCOL. You have three tools: ra_inspect (the disease/readout, trial "
    "timeline, drug formulary with mechanisms, and the calibrated MTX reference), "
    "ra_run_trial (run a protocol - first_line / second_line dose names, optional "
    "switch_day and dose_scale - and read the model's response rates), and "
    "ra_finalize (commit the protocol you recommend; that is what gets scored).\n\n"
    "Work like a scientist, not a search: the numbers are cheap but the protocol "
    "and its interpretation are yours to get right. Nothing tells you the correct "
    "drug, timing or dose - derive them. Start by VALIDATING your harness against "
    "the calibrated MTX reference so you trust the pipeline. Then reason from "
    "MECHANISM about the escalation. Watch how a protocol interacts with the "
    "model's day-284 MTX-non-responder classification and its day-600 readout, and "
    "be suspicious of any result that is clinically implausible or an artifact of "
    "how you set the trial up - diagnose it rather than reporting it. A run whose "
    "MTX-IR arm is empty is a mis-specified protocol, not a result; fix it. "
    "Finalize the protocol you would actually give a patient before finishing.\n\n")

_GOAL = {
    "predict": (
        "OBJECTIVE: predict the ACR20/50/70 response of patients who inadequately "
        "respond to first-line methotrexate when they are escalated to a "
        "second-line therapy. Choose the therapy and protocol, run it, and commit "
        "your predicted second-line rates. Judge for yourself whether your numbers "
        "are clinically plausible for the therapy you chose."),
    "min-dose": (
        "OBJECTIVE: for the methotrexate-inadequate responders, find the "
        "LOWEST-DOSE second-line regimen that still achieves ACR20 >= {thr:g}% in "
        "that subgroup. Use dose_scale to titrate the second-line amount (1.0 = the "
        "labeled dose). Simply maxing the dose meets the target but wastes drug and "
        "risks toxicity; too little misses it. Find the smallest dose_scale that "
        "clears the bar, then commit that protocol. Report the drug, the dose_scale, "
        "and the response you achieved."),
}


def _system_prompt(objective: str, thr: float) -> str:
    return _COMMON + _GOAL[objective].format(thr=thr)


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
    ap.add_argument("--objective", choices=["predict", "min-dose"], default="predict",
                    help="predict = reproduce the held-out TCZ-in-MTX-IR rates; "
                         "min-dose = find the lowest second-line dose clearing an "
                         "ACR20 target (the harder, non-brute-forceable task)")
    ap.add_argument("--min-dose-acr20", type=float, default=35.0,
                    help="min-dose objective: the ACR20%% bar to clear (default 35)")
    # scoring target for the predict objective
    ap.add_argument("--target-source", choices=["clinical", "model"], default="clinical",
                    help="clinical = score against the REAL TCZ trial (ROSE, ESM1); "
                         "model = score against the model's own reference output")
    ap.add_argument("--target-week", type=int, default=24, choices=[12, 24],
                    help="clinical target: which ROSE readout week (default 24)")
    # model-output target (used when --target-source model)
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

    model_target = {"ACR20": args.target_acr20, "ACR50": args.target_acr50,
                    "ACR70": args.target_acr70}
    clin_raw = ra_clin.trial_target("TCZ", args.target_week, "raw")
    clin_pcorr = ra_clin.trial_target("TCZ", args.target_week, "pcorr")
    target = clin_raw if args.target_source == "clinical" else model_target

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

        obj_text = (
            "Predict the ACR20/50/70 response of MTX-inadequate responders "
            "escalated to a second-line therapy."
            if args.objective == "predict" else
            f"Find the lowest-dose second-line regimen reaching ACR20 >= "
            f"{args.min_dose_acr20:g}% in the MTX-inadequate responders.")
        registry = ToolRegistry()
        register_ra_trial_loop_tools(registry, cfg, {
            "sb": sb, "vpop": args.vpop, "calibrated_arms": calibrated_arms,
            "objective": obj_text,
            "limit": args.limit, "stop_time": args.stop_time})

        goal = (obj_text + " Start by calling ra_inspect and validating the harness "
                "on first-line MTX, then design and run the second-line protocol.")
        policy = LLMPolicy(cfg, registry,
                           _system_prompt(args.objective, args.min_dose_acr20))
        loop = DecisionLoop(config=cfg, registry=registry, policy=policy)

        print(f"\n== LLM RA virtual-trial loop [{args.objective}] ({args.vpop}, "
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
        print(f"\n== SCORING [{args.objective}] ==")
        if not committed:
            print("[warn] agent did not call ra_finalize - scoring its last "
                  "non-empty run, which may be a control rather than its answer.")
        if not final:
            print("no scorable protocol run (the agent never produced a non-empty "
                  "second-line arm).")
        elif args.objective == "predict":
            pred = final["second_line"]
            print(f"agent protocol: {final['protocol']}")
            print(f"agent prediction (MTX-IR n={pred.get('n_MTX_IR')}): "
                  f"ACR20 {pred.get('ACR20')}% ACR50 {pred.get('ACR50')}% "
                  f"ACR70 {pred.get('ACR70')}%")
            src = ("REAL TCZ trial " + ra_clin.describe("TCZ") +
                   f", week {args.target_week}, raw"
                   if args.target_source == "clinical" else "the model's own output")
            score = osp_ra_trial.score_flagship(pred, target)
            print(f"\nPRIMARY score vs {src}:")
            print(f"   target {target}  ->  MAE {score.get('mae_pp')} pp "
                  f"over {score.get('n_endpoints')} endpoints")
            for k, v in (score.get("per_endpoint") or {}).items():
                print(f"     {k}: pred {v['predicted']} vs {v['target']} "
                      f"(|err| {v['abs_error_pp']} pp)")
            # transparency: show all three reference frames side by side
            print("\nreference frames (agent prediction vs each):")
            frames = [("model output", model_target),
                      (f"ROSE wk{args.target_week} raw", clin_raw),
                      (f"ROSE wk{args.target_week} placebo-corrected", clin_pcorr)]
            for name, ref in frames:
                if not ref:
                    continue
                mae = osp_ra_trial.score_flagship(pred, ref).get("mae_pp")
                cells = " ".join(f"{k}={ref.get(k)}" for k in ("ACR20", "ACR50", "ACR70"))
                print(f"   {name:38} {cells:34} MAE {mae} pp")
        else:  # min-dose
            pred = final["second_line"]
            # recover the committed dose_scale from the protocol string (NAME*scale@..)
            import re as _re
            m = _re.search(r"\*([0-9.]+)", final["protocol"])
            scale = float(m.group(1)) if m else 1.0
            score = osp_ra_trial.score_min_dose(pred, scale, args.min_dose_acr20)
            print(f"target: ACR20 >= {args.min_dose_acr20:g}% at the lowest dose")
            print(f"agent protocol: {final['protocol']}  (dose_scale={scale:g})")
            print(f"achieved: ACR20 {pred.get('ACR20')}% in MTX-IR "
                  f"n={pred.get('n_MTX_IR')}  -> {score['verdict']}")
            if score["target_met"]:
                print(f"response-per-dose (efficiency, higher=leaner): "
                      f"{score['response_per_dose']}  "
                      f"(a full-dose x1.0 answer scores {pred.get('ACR20')}; "
                      f"a leaner passing dose scores higher)")
    finally:
        sb.stop()


if __name__ == "__main__":
    main()
