r"""Stage-2b: Claude builds a PBPK model on the real OSP engine, closed loop.

Claude inspects a blanked PBPK model, edits parameters and structure, runs the
model headless in PK-Sim, reads the fit (GMFE + per-route bias), and revises -
until it hits the GMFE target or stops improving. Every run is a real PK-Sim
simulation; the score is measured against the observed clinical data.

    set ANTHROPIC_API_KEY=...
    set PKPD_PKSIM_CLI=C:\Program Files\Open Systems Pharmacology\PK-Sim 12.3\PKSim.CLI.exe

    python -m examples.run_llm_task ^
        --snapshot ..\OSP-PBPK-Model-Library\Alfentanil\benchmark\Alfentanil-Model.blanked.json ^
        --input    ..\OSP-PBPK-Model-Library\Alfentanil\json_input\Alfentanil-Model.input.json ^
        --target 1.6 --max-steps 10

Each osp_try_model call runs all 12 simulations (a few minutes), so keep
--max-steps modest.
"""

from __future__ import annotations

import argparse
import json
import os

from pkpd_agent.config import AgentConfig
from pkpd_agent.engines.osp_cli import OSPCli
from pkpd_agent.llm import LLMPolicy
from pkpd_agent.loop import DecisionLoop
from pkpd_agent.state import Decision, Finish, ModelingSession, Observation
from pkpd_agent.tools.registry import ToolRegistry
from pkpd_agent.tools.osp_loop_tools import register_osp_loop_tools


def _system_prompt(target: float) -> str:
    return (
        "You are a PBPK modeler using the OSP PK-Sim engine. You are given a "
        "whole-body physiologically based model that is set up but NOT yet "
        "fitted, and clinical plasma concentration-time data. Your job is to "
        "edit the model so simulated plasma concentrations match the observed "
        "data.\n\n"
        "Workflow:\n"
        "1. Call osp_inspect once to see the current parameters, distribution/"
        "permeability methods, processes, the literature priors, and the "
        "observed data (routes/studies/doses).\n"
        "2. Reason mechanistically, then call osp_try_model with an edit spec. "
        "Use the per-route BIAS it returns to decide direction:\n"
        "   - bias > 1 means the model OVER-predicts concentration -> exposure "
        "too high -> increase clearance (or lower bioavailability/permeability "
        "for oral).\n"
        "   - bias < 1 means UNDER-predicts -> decrease clearance, or for oral "
        "raise intestinal permeability.\n"
        "   - IV arms isolate distribution + clearance; oral arms add "
        "absorption. Fix IV (clearance, distribution) before oral (absorption).\n"
        "3. Change only a few parameters at a time so you can attribute the "
        "effect. Keep every value physically plausible (fraction unbound in "
        "(0,1]; clearance below hepatic blood flow; positive permeabilities).\n"
        f"4. Iterate until the overall GMFE is <= {target} or it stops "
        "improving, then finish with a short summary of your final model "
        "(the parameter values and any structural choices) and the GMFE.\n\n"
        "Be decisive and quantitative: a fold-change in clearance moves exposure "
        "roughly the inverse fold. Do not call osp_try_model with no edits."
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--snapshot", required=True, help="blanked snapshot to build from")
    ap.add_argument("--input", required=True, help="clean-input JSON (observed + objective)")
    ap.add_argument("--pksim", default=None, help="PKSim.CLI.exe (else PKPD_PKSIM_CLI)")
    ap.add_argument("--target", type=float, default=1.6, help="GMFE target")
    ap.add_argument("--max-steps", type=int, default=10)
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

    cli = OSPCli(pksim_cli_path=args.pksim or cfg.pksim_cli_path,
                 timeout_s=cfg.pksim_timeout_s)
    if not cli.pksim_cli_path or not os.path.exists(cli.pksim_cli_path):
        print(f"PKSim.CLI not found at {cli.pksim_cli_path!r}; set --pksim or "
              "PKPD_PKSIM_CLI.")
        return

    with open(args.input, encoding="utf-8") as fh:
        inp = json.load(fh)
    observed = inp["given_data"]["clinical_observed_data"]

    registry = ToolRegistry()
    register_osp_loop_tools(registry, cfg, {
        "cli": cli, "snapshot_path": args.snapshot,
        "observed": observed, "input": inp})

    goal = (f"{inp.get('objective','Fit the PBPK model to the observed data.')}\n\n"
            f"Target: overall GMFE <= {args.target}. Start by calling osp_inspect.")
    policy = LLMPolicy(cfg, registry, _system_prompt(args.target))
    loop = DecisionLoop(config=cfg, registry=registry, policy=policy)

    print(f"== LLM PBPK loop on {os.path.basename(args.snapshot)} "
          f"(target GMFE {args.target}, max {args.max_steps} steps) ==\n")
    session = ModelingSession(goal=goal)
    session = loop.run(goal, session)

    # replay the decision/observation trace
    for ev in session.transcript:
        if isinstance(ev, Decision):
            if ev.text:
                print(f"\n[reason] {ev.text[:600]}")
            for c in ev.calls:
                ed = c.arguments.get("edits")
                print(f"  -> {c.name}" + (f" {json.dumps(ed, ensure_ascii=False)}"
                                          if ed else ""))
        elif isinstance(ev, Observation):
            msg = ev.content.get("message", "")
            print(f"  <- {ev.tool}: {msg}")
            br = ev.content.get("by_route")
            if br:
                for r, m in br.items():
                    print(f"       {r}: GMFE {m.get('gmfe')} bias {m.get('bias')} "
                          f"within2fold {m.get('within_2fold_pct')}%")
            for f in ev.findings:
                print(f"     [{f.level.upper()}] {f.message}")
        elif isinstance(ev, Finish):
            print(f"\n=== SUMMARY ===\n{ev.text}")

    print(f"\nbest GMFE reached: {session.get('osp_best_gmfe')}")
    print(f"best edits: {json.dumps(session.get('osp_best_edits'), ensure_ascii=False)}")


if __name__ == "__main__":
    main()
