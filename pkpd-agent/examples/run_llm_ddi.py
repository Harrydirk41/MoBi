r"""LLM-driven DDI prediction loop: the agent tunes the interaction, PK-Sim runs it.

Where run_llm_build fits a single compound's disposition, this drives the agent
over a drug-drug interaction: it reads the interaction structure (perpetrator,
victim, mechanism), decides/estimates the interaction parameters (Ki, kinact /
K_kinact_half, EC50 / Emax), runs every control/treatment arm headless, and
checks the predicted vs observed interaction ratio (AUCR). The victim's own model
is fixed; the unknown is the interaction.

    set ANTHROPIC_API_KEY=...
    set PKPD_PKSIM_CLI=C:\Program Files\Open Systems Pharmacology\PK-Sim 12.3\PKSim.CLI.exe

    python -m examples.run_llm_ddi ^
        --snapshot ..\OSP-PBPK-Model-Library\Erythromycin\benchmark\Erythromycin-Model.ddi_blanked.json ^
        --input    ..\OSP-PBPK-Model-Library\Erythromycin\json_input\Erythromycin-Model.ddi_input.json ^
        --victim Midazolam --target 1.5 --max-steps 6
"""

from __future__ import annotations

import argparse
import json
import os

from pkpd_agent.config import AgentConfig
from pkpd_agent.engines.osp_cli import OSPCli
from pkpd_agent.engines import osp_ddi
from pkpd_agent.llm import LLMPolicy
from pkpd_agent.loop import DecisionLoop
from pkpd_agent.state import Decision, Finish, ModelingSession, Observation
from pkpd_agent.tools.registry import ToolRegistry
from pkpd_agent.tools.osp_ddi_loop_tools import register_osp_ddi_loop_tools


def _system_prompt(target: float) -> str:
    return (
        "You are a PBPK modeler predicting a drug-drug interaction (DDI) with the "
        "OSP PK-Sim engine. The victim and perpetrator disposition models are "
        "GIVEN and fixed; you determine only the INTERACTION.\n\n"
        "Each round:\n"
        "1. Call ddi_inspect: the perpetrator(s), the victim, each interaction "
        "mechanism (competitive inhibition -> Ki; mechanism-based inhibition -> "
        "kinact & K_kinact_half; induction -> EC50 & Emax) with its current "
        "parameters, the control/treatment pairs, the observed interaction ratios "
        "(AUCR), and the identifiability guidance.\n"
        "2. DECIDE the interaction parameters. Use in-vitro / literature values "
        "where you have them (a reported Ki, kinact). Act on the identifiability "
        "guidance: a two-parameter mechanism (kinact/K_kinact_half, EC50/Emax) is "
        "underdetermined by a SINGLE ratio - fix the better-known parameter (e.g. "
        "K_kinact_half at its in-vitro value) and estimate the other, rather than "
        "floating both.\n"
        "3. Call ddi_try_model with interaction_parameters to run the arms and get "
        "the predicted vs observed AUCR.\n"
        "4. CHECK: is the AUCR GMFE acceptable? Which arms are off, and in which "
        "direction (predicted interaction too weak -> increase inhibition potency, "
        "i.e. LOWER Ki / RAISE kinact)? Iterate.\n"
        "ENGINE ERRORS ARE HARD FAILURES: if a run errors, the model did not run - "
        "fix the offending edit and re-run; do not finish on a failed run.\n"
        f"5. Iterate until the AUCR GMFE <= {target}, then finish with a summary: "
        "the mechanism, the interaction parameters you set/estimated, and the GMFE."
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--snapshot", required=True)
    ap.add_argument("--input", required=True)
    ap.add_argument("--victim", default=None)
    ap.add_argument("--pksim", default=None)
    ap.add_argument("--target", type=float, default=1.5)
    ap.add_argument("--max-steps", type=int, default=6)
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
        print(f"PKSim.CLI not found at {cli.pksim_cli_path!r}.")
        return

    with open(args.input, encoding="utf-8") as fh:
        inp = json.load(fh)
    with open(args.snapshot, encoding="utf-8") as fh:
        snap = json.load(fh)
    ddi = osp_ddi.analyze_ddi(snap)
    if not ddi:
        print("Not a DDI snapshot.")
        return
    victim = args.victim or (ddi["victims"][0] if ddi.get("victims") else None)
    observed_ratios = ((inp.get("given_data") or {}).get("observed_interaction_ratios")
                       or [])

    registry = ToolRegistry()
    register_osp_ddi_loop_tools(registry, cfg, {
        "cli": cli, "snapshot_path": args.snapshot, "ddi": ddi,
        "victim": victim, "observed_ratios": observed_ratios, "input": inp})

    goal = (f"{inp.get('objective', 'Predict the DDI.')}\n\n"
            f"Target: AUCR GMFE <= {args.target}. Start by calling ddi_inspect, "
            "then set the interaction parameters and call ddi_try_model.")
    policy = LLMPolicy(cfg, registry, _system_prompt(args.target))
    loop = DecisionLoop(config=cfg, registry=registry, policy=policy)

    print(f"== LLM DDI loop on {os.path.basename(args.snapshot)} "
          f"(perpetrator -> {victim}, target AUCR GMFE {args.target}) ==\n")

    def show(ev):
        if isinstance(ev, Decision):
            if ev.text:
                print(f"\n[reason] {ev.text[:2000]}")
            for c in ev.calls:
                print(f"  -> {c.name} {json.dumps(c.arguments, ensure_ascii=False)[:400]}")
                if c.name == "ddi_try_model":
                    print("     ...running control + treatment arms...", flush=True)
        elif isinstance(ev, Observation):
            c = ev.content
            print(f"  <- {ev.tool}: {c.get('message', '')}", flush=True)
            for a in (c.get("per_arm") or [])[:4]:
                print(f"       {a.get('treatment')}: pred {a.get('predicted_aucr')} "
                      f"obs {a.get('observed_aucr')} fold {a.get('fold_error')}")
        elif isinstance(ev, Finish):
            print(f"\n=== SUMMARY ===\n{ev.text}")

    session = loop.run(goal, ModelingSession(goal=goal), on_event=show)
    print(f"\nbest AUCR GMFE reached: {session.get('ddi_best_gmfe')}")
    print(f"best interaction: "
          f"{json.dumps(session.get('ddi_best_edits'), ensure_ascii=False)}")


if __name__ == "__main__":
    main()
