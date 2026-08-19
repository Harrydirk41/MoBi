r"""LLM-driven RA Stage-4 CALIBRATION: the agent fits a PD parameter to trial data.

Where run_llm_ra_trial designs a protocol, this drives the agent over the paper's
own hard step - parameter estimation. The drug arm is fixed (MTX then TCZ in the
MTX-inadequate responders); the tocilizumab potency parameter KD_TCZ is treated as
UNKNOWN, and the agent estimates it so the model reproduces the OBSERVED ROSE-trial
ACR. Same inverse-problem loop as the OSP PK fitting, on the RA QSP model.

    set ANTHROPIC_API_KEY=...
    python -m examples.run_llm_ra_fit ^
        --sbproj "..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj" ^
        --vpop   "..\RA-QSP-Model\Vpop1.xlsx" --limit 50 --max-steps 12

Each fit_try simulates the whole (subsampled) population, so keep --limit modest.
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
from pkpd_agent.tools.ra_fit_loop_tools import register_ra_fit_loop_tools


def _system_prompt(param: str) -> str:
    return (
        "You are a QSP modeler CALIBRATING a rheumatoid-arthritis model with the "
        "SimBiology engine. The model structure, the virtual population and the "
        "treatment protocol are all FIXED. One PD parameter is UNKNOWN and you must "
        f"estimate it - here, {param}, tocilizumab's binding affinity - so that the "
        "model's simulated response reproduces an OBSERVED clinical trial result.\n\n"
        "This is a parameter-identification problem, so work like one:\n"
        "1. Call fit_inspect: the fixed arm, the parameter (its unit, meaning, and "
        "plausible range), and the observed target ACR you must hit.\n"
        "2. Reason about direction from the biology: for a dissociation constant, a "
        "SMALLER value means tighter binding and a STRONGER effect. Bracket the "
        "target first with a low and a high value, confirm the response moves the "
        "way you expect, then converge.\n"
        "3. Call fit_try with {overrides:{" + param + ": value}} and read the ACR "
        "error (MAE) vs the target. Because the parameter is log-scale, step by "
        "FACTORS (10x, 3x, ...) not by small increments; then bisect in log space.\n"
        "4. Iterate until the ACR MAE is small and stable. Watch for a plateau - if "
        "the response saturates, a wide range of values fits equally and you should "
        "report the identifiable edge, not an arbitrary point in the flat region.\n"
        "5. Call fit_finalize with your committed value, then finish with: the "
        "estimated value, how tightly it is identified, and how it compares to the "
        "literature reference."
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sbproj", required=True)
    ap.add_argument("--vpop", required=True)
    ap.add_argument("--param", default="KD_TCZ",
                    help="PD parameter to calibrate (default KD_TCZ)")
    ap.add_argument("--arm", default="MTX_15mg_Q1W_SC_t200;TCZ8mgkg_Q4W_IV_t200@285",
                    help="fixed drug arm (the protocol is not the unknown here)")
    ap.add_argument("--target-week", type=int, default=24, choices=[12, 24])
    ap.add_argument("--limit", type=int, default=50)
    ap.add_argument("--stop-time", type=float, default=700.0)
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

    target = ra_clin.trial_target("TCZ", args.target_week, "raw")

    sb = SimBiologyEngine()
    try:
        print("== starting MATLAB/SimBiology engine ==", flush=True)
        sb.start()
        print(f"== loading {os.path.basename(args.sbproj)} ==", flush=True)
        sb.load_project(args.sbproj)

        registry = ToolRegistry()
        register_ra_fit_loop_tools(registry, cfg, {
            "sb": sb, "vpop": args.vpop, "arm": args.arm, "target": target,
            "fit_params": [args.param], "limit": args.limit,
            "stop_time": args.stop_time})

        goal = (f"Estimate {args.param} so the model reproduces the observed ACR "
                f"{target}. Start with fit_inspect, bracket the target, then converge "
                "and finalize.")
        policy = LLMPolicy(cfg, registry, _system_prompt(args.param))
        loop = DecisionLoop(config=cfg, registry=registry, policy=policy)

        print(f"\n== LLM RA calibration loop [{args.param}] "
              f"(target ROSE wk{args.target_week} {target}) ==\n", flush=True)

        def show(ev):
            if isinstance(ev, Decision):
                if ev.text:
                    print(f"\n[reason] {ev.text[:2000]}")
                for c in ev.calls:
                    print(f"  -> {c.name} "
                          f"{json.dumps(c.arguments, ensure_ascii=False)[:300]}")
                    if c.name == "fit_try":
                        print("     ...simulating the population...", flush=True)
            elif isinstance(ev, Observation):
                print(f"  <- {ev.tool}: {ev.content.get('message', '')}", flush=True)
            elif isinstance(ev, Finish):
                print(f"\n=== AGENT SUMMARY ===\n{ev.text}")

        session = loop.run(goal, ModelingSession(goal=goal), on_event=show)

        print("\n== SCORING [calibration] ==")
        final = session.get("fit_final")
        if not final:
            best = session.get("fit_best_overrides")
            print("[warn] agent did not fit_finalize; best-MAE run was", best)
            final = None
        if final:
            sc = final["score"]
            print(f"target (observed ROSE wk{args.target_week}): {target}")
            print(f"committed fit: {final['overrides']}")
            print(f"predicted ACR: {final['predicted']}")
            print(f"fit quality: ACR MAE {sc.get('acr_mae_pp')} pp vs observed")
            for name, p in (sc.get("parameters") or {}).items():
                lf = p.get("log10_fold_from_ref")
                print(f"parameter recovery: {name} fitted {p.get('fitted'):.3g} vs "
                      f"literature {p.get('reference'):.3g}  "
                      f"({'%+.2f' % lf if lf is not None else 'n/a'} log10-fold)")
    finally:
        sb.stop()


if __name__ == "__main__":
    main()
