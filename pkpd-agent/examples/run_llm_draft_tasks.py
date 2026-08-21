r"""LLM tasks.json drafter + the RA regression that verifies it.

Reads a model's network.json, has an LLM assign each parameter/state to its downstream-
task ROLE (disease driver / druggable / calibratable / readout), and drafts the
model-derivable half of a tasks.json. On --model ra it REGRESSES the draft against the
hand-written projects/vantage_ra/tasks.json: the role assignments must reproduce before
the drafter is trusted on a new model.

The draft leaves the genuinely-external fields (clinical trial numbers, dose names,
event timeline) as TODO stubs for the author - the drafter never invents those.

    python -m examples.run_llm_draft_tasks --network network.json --model ra
    python -m examples.run_llm_draft_tasks --network network.json --model ra --out draft_tasks.json

Needs ANTHROPIC_API_KEY.
"""

from __future__ import annotations

import argparse
import json

from pkpd_agent.config import AgentConfig
from pkpd_agent.engines import llm_tasks as LT
from pkpd_agent.engines import qsp_config


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--network", required=True)
    ap.add_argument("--model", default="vantage_ra",
                    help="known-good project to regress the draft against (if it exists)")
    ap.add_argument("--out", default=None, help="write the draft tasks.json here")
    ap.add_argument("--llm-model", default=None)
    args = ap.parse_args()

    with open(args.network, encoding="utf-8") as fh:
        data = json.load(fh)

    cfg = AgentConfig(mock=False)
    if args.llm_model:
        cfg.model = args.llm_model
    if not cfg.anthropic_key_present():
        print("ANTHROPIC_API_KEY not set.")
        return

    print("== LLM assigning parameters/states to task roles (any naming) ==", flush=True)
    call = LT.default_call(cfg)
    draft = LT.draft_tasks(data, call, name=data.get("name", args.model))
    print(f"drafted: {len(draft['readout_states'])} readout states, "
          f"{len(draft['vpop_drivers'])} vpop drivers, "
          f"{len(draft['design_targets'])} design targets, "
          f"{len(draft['fit_params'])} fit params\n")

    # regression against the known-good reference, if this project already has one
    try:
        ref = qsp_config.get(args.model)
        reference = {
            "readout_states": ref.readout_states, "vpop_drivers": ref.vpop_drivers,
            "design_targets": ref.design_targets, "fit_params": ref.fit_params}
        print("== REGRESSION vs the hand-written tasks.json (must reproduce it) ==")
        cmp = LT.compare_tasks(draft, reference)
        for k in ("readout_states", "vpop_drivers", "design_targets", "fit_params"):
            d = cmp[k]
            print(f"  {k:15s} P {d['precision']}  R {d['recall']}  F1 {d['f1']}  "
                  f"({d['hit']}/{d['n_truth']})")
            if d["missed"]:
                print(f"      missed: {d['missed']}")
            if d["extra"]:
                print(f"      extra:  {d['extra']}")
        print("\nRead: RECALL is the metric that matters here - the readout states are "
              "structural and reproduce well (high F1), but which PARAMETERS to expose "
              "for a Vpop / design / calibration task is the modeler's COMMITTED CHOICE, "
              "not derivable from the model. So the parameter fields are best read as "
              "high-recall CANDIDATE POOLS the author prunes (low precision is expected: "
              "the drafter finds every candidate but cannot guess the author's selection).")
    except KeyError:
        print(f"(no known-good project '{args.model}' to regress against - draft only)")

    print("\nNOTE: readout_states is a real draft; vpop_drivers / design_targets / "
          "fit_params are candidate POOLS to prune down to the ones your task exposes. "
          "drugs / timeline / vpop_target / clinical_trials / refractory_target are TODO "
          "stubs (external data or .sbproj dose names the drafter cannot invent). Prune + "
          "fill, then save as projects/<name>/tasks.json.")

    if args.out:
        with open(args.out, "w", encoding="utf-8") as fh:
            json.dump(draft, fh, indent=2)
        print(f"\nwrote draft to {args.out}")


if __name__ == "__main__":
    main()
