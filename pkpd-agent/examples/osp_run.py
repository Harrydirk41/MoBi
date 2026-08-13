r"""Stage-1 proof: run the REAL PBPK model headless and grade it.

Pipeline (no LLM, no GUI):
    snapshot JSON  --PKSim.CLI snap+export-->  simulated profiles
    map simulations -> observed datasets by (study, route, dose)
    build a submission  ->  grade_submission.py  ->  GMFE + plausibility

This proves the whole engine + grader on real PK-Sim output before we wire the
LLM loop on top. Run it on the reference snapshot first: the GMFE should be
good (this IS the published model), which validates the plumbing end to end.

    set PKPD_PKSIM_CLI=C:\Program Files\Open Systems Pharmacology\PK-Sim 12.3\PKSim.CLI.exe
    python -m examples.osp_run ^
        --snapshot ..\OSP-PBPK-Model-Library\Alfentanil\json\Alfentanil-Model.json ^
        --input    ..\OSP-PBPK-Model-Library\Alfentanil\json_input\Alfentanil-Model.input.json ^
        --key      ..\OSP-PBPK-Model-Library\Alfentanil\answer_key\Alfentanil-Model.answer_key.json

    # tune parameters and re-run (the agent's action, done by hand):
    #   --overrides "{\"Lipophilicity\": 2.2, \"Intrinsic clearance\": 0.4}"
"""

from __future__ import annotations

import argparse
import json
import os
import sys

from pkpd_agent.config import AgentConfig
from pkpd_agent.engines.osp_cli import OSPCli
from pkpd_agent.engines.osp_score import map_predictions

_HERE = os.path.dirname(os.path.abspath(__file__))
_ALF = os.path.normpath(os.path.join(_HERE, "..", "..",
                                     "OSP-PBPK-Model-Library", "Alfentanil"))
sys.path.insert(0, _ALF)
import grade_submission as G          # noqa: E402  (path injected above)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--snapshot", required=True, help="OSP snapshot JSON to run")
    ap.add_argument("--input", required=True, help="clean-input JSON (observed + rubric)")
    ap.add_argument("--key", default=None, help="answer key (auxiliary closeness)")
    ap.add_argument("--overrides", default=None,
                    help='JSON dict of compound param overrides, '
                         'e.g. {"Lipophilicity": 2.0}')
    ap.add_argument("--edits", default=None,
                    help='Full edit spec: a JSON file path or inline JSON with '
                         'keys parameters / calculation_methods / processes')
    ap.add_argument("--pksim", default=None, help="path to PKSim.CLI.exe "
                    "(else config / PKPD_PKSIM_CLI)")
    ap.add_argument("--reason", action="store_true", help="run the Claude judge too")
    ap.add_argument("--submission", default="osp_submission.json")
    ap.add_argument("--keep-workdir", action="store_true")
    args = ap.parse_args()

    cfg = AgentConfig(mock=False)
    cli = OSPCli(pksim_cli_path=args.pksim or cfg.pksim_cli_path,
                 timeout_s=cfg.pksim_timeout_s, keep_workdir=args.keep_workdir)

    overrides = json.loads(args.overrides) if args.overrides else None
    edits = None
    if args.edits:
        if os.path.exists(args.edits):
            with open(args.edits, encoding="utf-8") as fh:
                edits = json.load(fh)
        else:
            edits = json.loads(args.edits)
    print(f"running PK-Sim on {os.path.basename(args.snapshot)} ...")
    res = cli.build_and_run(args.snapshot, edits=edits, param_overrides=overrides)
    if not res["ok"]:
        print("ENGINE FAILED:", res["message"])
        for lg in res.get("logs", []):
            print(f"  [{lg['cmd']}] rc={lg['returncode']} {lg['stderr'][:300]}")
        sys.exit(1)
    print(f"  {res['message']}")
    applied = res.get("edits_applied") or {}
    if applied.get("parameters"):
        print(f"  parameters set: {applied['parameters']}")
    if applied.get("calculation_methods"):
        print(f"  methods set: {applied['calculation_methods']}")
    if applied.get("processes"):
        print(f"  processes: {applied['processes']}")
    if applied.get("not_found"):
        print(f"  NOT FOUND (ignored): {applied['not_found']}")

    # map simulations -> observed datasets by (study, route, dose)
    with open(args.input, encoding="utf-8") as fh:
        observed = json.load(fh)["given_data"]["clinical_observed_data"]
    predicted_profiles, unmatched_obs = map_predictions(res["profiles"], observed)
    print(f"  matched {len(predicted_profiles)}/{len(observed)} "
          "observed datasets to simulations")
    if unmatched_obs:
        print(f"  UNMATCHED observed ({len(unmatched_obs)}): "
              f"{unmatched_obs[:3]}{' ...' if len(unmatched_obs) > 3 else ''}")

    submission = {
        "submission": {
            "structural_model": {
                "engine": "OSP PK-Sim (PKSim.CLI, headless from snapshot)",
                "source_snapshot": os.path.basename(args.snapshot),
                "notes": "whole-body PBPK; simulations run as defined in the snapshot",
            },
            "parameters": [
                {"parameter": k, "value": v, "fixed_or_estimated": "estimated",
                 "rationale": "edit applied to snapshot"}
                for k, v in applied.get("parameters", {}).items()
            ],
            "predicted_profiles": predicted_profiles,
        }
    }
    with open(args.submission, "w", encoding="utf-8") as fh:
        json.dump(submission, fh, ensure_ascii=False, indent=2)

    # grade (reuse the committed grader; --reason adds the Claude judge)
    reason = True if args.reason else False
    sc = G.grade(args.input, args.submission, args.key, reason,
                 cfg.model, cfg.effort)
    G._print_summary(sc)
    out = os.path.join(_ALF, "scorecards",
                       os.path.basename(args.submission).replace(".json", "")
                       + ".scorecard.json")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(sc, fh, indent=2, ensure_ascii=False)
    print(f"wrote {args.submission} and {out}")


if __name__ == "__main__":
    main()
