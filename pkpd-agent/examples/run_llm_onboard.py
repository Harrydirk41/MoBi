r"""Conversational onboarding: a modeler describes the model, an AGENT builds the project.

The GUI-free entry point. The modeler writes (or dictates) a few plain sentences about
their model; an agent drives the whole build itself - inspects the model, builds the
config from the description, reads its own validation report, fixes what it can infer,
and saves projects/<name>/{spec,tasks}.json, reporting in plain English what still needs
the modeler. No JSON editing, no schema, no hand-run pipeline steps.

    set ANTHROPIC_API_KEY=...
    python -m examples.run_llm_onboard --network network.json --name my_model ^
        --describe describe.txt

`describe.txt` is free text, e.g.: "Psoriasis QSP. Severity readout PASI, active band
6-20. Drug secukinumab, dose SEC_300mg_Q4W. Match UNCOVER-2 PASI75 at week 12 = 77%.
Disease drivers are the F_* amplification factors; calibrate KD_SEC (reference 1e-10 M)."

Needs ANTHROPIC_API_KEY. This is config-building only - no MATLAB. After it saves, run
the task agents (run_llm_qsp_full --model <name> ...).
"""

from __future__ import annotations

import argparse
import json
import os

from pkpd_agent.config import AgentConfig
from pkpd_agent.engines import llm_tasks as LT
from pkpd_agent.llm import LLMPolicy
from pkpd_agent.loop import DecisionLoop
from pkpd_agent.state import Decision, Finish, ModelingSession, Observation
from pkpd_agent.tools.registry import ToolRegistry
from pkpd_agent.tools.onboarding_loop_tools import register_onboarding_tools


def _system_prompt(name: str) -> str:
    return (
        "You are onboarding a QSP model for a modeler who does not edit JSON. You have "
        "four tools: onboard_inspect (see the model), onboard_build (build the config "
        "from the modeler's description and validate it), onboard_set (fix one flagged "
        "field by dotted path), and onboard_save (write the project once there are no "
        "ERRORS).\n\nWork like this: inspect the model, build the config, then READ the "
        "validation report. Fix every ERROR you can infer from the model or the "
        "description with onboard_set (e.g. a vpop_driver that is not a real parameter "
        "is usually a near-miss of a real one - correct it). Leave WARNINGS about "
        "external clinical data the description did not give - those are the modeler's "
        "to fill. Save when clean, then tell the modeler in plain language what you "
        f"built for '{name}' and exactly which numbers they still need to provide.")


def _read_describe(arg: str) -> str:
    if arg and os.path.isfile(arg):
        with open(arg, encoding="utf-8") as fh:
            return fh.read()
    return arg or ""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--network", required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--describe", default="")
    ap.add_argument("--max-steps", type=int, default=12)
    ap.add_argument("--llm-model", dest="llm_model", default=None)
    ap.add_argument("--out-dir", default=None)
    args = ap.parse_args()

    with open(args.network, encoding="utf-8") as fh:
        network = json.load(fh)
    description = _read_describe(args.describe)

    cfg = AgentConfig(mock=False, max_steps=args.max_steps)
    if args.llm_model:
        cfg.model = args.llm_model
    if not cfg.anthropic_key_present():
        print("ANTHROPIC_API_KEY not set.")
        return

    registry = ToolRegistry()
    register_onboarding_tools(registry, cfg, {
        "network": network, "description": description, "name": args.name,
        "call": LT.default_call(cfg), "out_dir": args.out_dir})

    goal = (f"Onboard the model as project '{args.name}': inspect it, build the config "
            "from this description, fix any validation errors, and save. "
            f"DESCRIPTION:\n{description or '(none given - build a skeleton to fill)'}")
    policy = LLMPolicy(cfg, registry, _system_prompt(args.name))
    loop = DecisionLoop(config=cfg, registry=registry, policy=policy)

    print(f"\n== ONBOARDING agent for '{args.name}' ==\n", flush=True)

    def show(ev):
        if isinstance(ev, Decision):
            if ev.text:
                print(f"\n[reason] {ev.text[:1500]}")
            for c in ev.calls:
                print(f"  -> {c.name} {json.dumps(c.arguments, ensure_ascii=False)[:200]}")
        elif isinstance(ev, Observation):
            print(f"  <- {ev.tool}: {ev.content.get('message', '')}", flush=True)
        elif isinstance(ev, Finish):
            print(f"\n=== AGENT SUMMARY ===\n{ev.text}")

    session = loop.run(goal, ModelingSession(goal=goal), on_event=show)
    saved = session.get("onboard_saved")
    print(f"\n{'saved to ' + saved if saved else 'not saved (see report above)'}")


if __name__ == "__main__":
    main()
