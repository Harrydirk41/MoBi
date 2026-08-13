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
import re
import sys

from pkpd_agent.config import AgentConfig
from pkpd_agent.engines.osp_cli import OSPCli

_HERE = os.path.dirname(os.path.abspath(__file__))
_ALF = os.path.normpath(os.path.join(_HERE, "..", "..",
                                     "OSP-PBPK-Model-Library", "Alfentanil"))
sys.path.insert(0, _ALF)
import grade_submission as G          # noqa: E402  (path injected above)


_MG = {"µg": 1e-3, "ug": 1e-3, "mg": 1.0, "g": 1e3}


def _norm_study(s: str | None) -> str:
    """Author+year key so 'Kharasch2012_Alfentanil_alone_IV', 'Kharasch 2012'
    collapse to the same thing while 'Kharasch 2011' vs '2011b' stay distinct."""
    m = re.search(r"([A-Za-z]+)\s*(\d{4}[a-z]?)", s or "")
    if m:
        return (m.group(1) + m.group(2)).lower()
    return re.sub(r"[^a-z0-9]", "", (s or "").lower())


def _norm_route(s) -> str | None:
    if not s:
        return None
    u = str(s).upper()
    return "PO" if u in ("PO", "ORAL") else ("IV" if "IV" in u else u)


def _dose_canon(s) -> tuple[float, bool] | None:
    """Canonicalize a dose string to (milligrams, per_kg) so that
    '20 µg/kg' and '0.02 mg/kg' compare equal."""
    if not s:
        return None
    m = re.search(r"([\d.]+)\s*(µg|ug|mg|g)\b\s*(/?\s*kg)?", str(s), re.IGNORECASE)
    if not m:
        return None
    return float(m.group(1)) * _MG[m.group(2).lower()], bool(m.group(3))


def _obs_key(o: dict):
    """(study, route, dose) for an observed dataset, falling back to parsing
    the dataset name when the metadata fields are empty (e.g. Kharasch 2012)."""
    from pkpd_agent.engines.osp_cli import OSPCli
    study, route, dose = o.get("study"), o.get("route"), o.get("dose")
    if not (study and route):
        ps, pr, pd = OSPCli._parse_sim_name(o.get("dataset", ""))
        study = study or ps
        route = route or pr
        dose = dose or pd
    return study, route, dose


def _match_score(obs, pred) -> int | None:
    """None if incompatible; higher = more specific match."""
    o_study, o_route, o_dose = _obs_key(obs)
    if _norm_study(o_study) != _norm_study(pred.study):
        return None
    score = 1
    o_r, p_r = _norm_route(o_route), _norm_route(pred.route)
    if o_r and p_r:
        if o_r != p_r:
            return None
        score += 2
    o_d, p_d = _dose_canon(o_dose), _dose_canon(pred.dose)
    if o_d and p_d:
        if o_d[1] != p_d[1] or abs(o_d[0] - p_d[0]) > 1e-6 + 0.01 * o_d[0]:
            return None
        score += 2
    return score


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
    preds = res["profiles"]
    predicted_profiles = []
    unmatched_obs, matched = [], 0
    for o in observed:
        best, best_score = None, 0
        for p in preds:
            s = _match_score(o, p)
            if s and s > best_score:
                best, best_score = p, s
        if best:
            predicted_profiles.append({"dataset": o["dataset"],
                                       "time_h": best.time_h,
                                       "pred_conc_mg_L": best.conc_mg_L,
                                       "_from_simulation": best.simulation})
            matched += 1
        else:
            unmatched_obs.append(o["dataset"])
    print(f"  matched {matched}/{len(observed)} observed datasets to simulations")
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
