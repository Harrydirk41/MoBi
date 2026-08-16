r"""Generate a biologic (protein / mAb) benchmark from an OSP snapshot.

The large-molecule analogue of build_benchmark.py. A biologic model's ANSWER is
its protein-disposition parameter(s) - the FcRn recycling affinity
``Kd (FcRn) in endosomal space`` and/or the renal ``GFR fraction`` - identified
from a biodistribution (concentration in whole blood and tissues). There are no
enzymes to fit; the structure is PK-Sim's large-molecule disposition and it is
GIVEN. This blanks the fitted disposition parameter(s) to a no-leak prior and
writes the task/answer files. The benchmark quantity is the fit across every
measured matrix.

Writes:
    <Compound>/json_input/<stem>.bio_input.json        agent task
    <Compound>/benchmark/<stem>.bio_blanked.json       value-blanked start
    <Compound>/answer_key/<stem>.bio_answer_edits.json fitted disposition params
    <Compound>/answer_key/<stem>.bio_answer_key.json   + observed matrices

    python -m examples.build_biologic_benchmark ^
        ..\OSP-PBPK-Model-Library\BAY794620\json\BAY794620.json --write
    python -m examples.build_biologic_benchmark --audit ..\OSP-PBPK-Model-Library
"""

from __future__ import annotations

import argparse
import copy
import glob
import json
import os

from pkpd_agent.engines import osp_biologic as B

# no-leak naive priors for the biologic disposition parameters
BIOLOGIC_DEFAULTS = {"Kd (FcRn) in endosomal space": 1.0, "GFR fraction": 1.0,
                     "Radius (solute)": 0.005, "Kd": 1.0, "koff": 1.0}


def _blank_default(name: str) -> float:
    return BIOLOGIC_DEFAULTS.get(name, 1.0)


def _blank(snapshot: dict, param_names: set[str]) -> dict:
    snap = copy.deepcopy(snapshot)

    def w(o):
        if isinstance(o, dict):
            nm = o.get("Name")
            vo = o.get("ValueOrigin") or {}
            if isinstance(nm, str) and isinstance(o.get("Value"), (int, float)) \
                    and vo.get("Source") == "ParameterIdentification":
                o["Value"] = _blank_default(nm)
                o["ValueOrigin"] = {"Source": "Unknown",
                                    "Description": "benchmark naive prior (blanked)"}
            for v in o.values():
                w(v)
        elif isinstance(o, list):
            for v in o:
                w(v)

    for c in snap.get("Compounds") or []:
        w(c)
    return snap


def build(path: str) -> dict:
    with open(path, encoding="utf-8") as fh:
        snap = json.load(fh)
    stem = os.path.splitext(os.path.basename(path))[0]
    bstruct = B.analyze_biologic(snap)
    if not bstruct:
        return {"skip": True, "stem": stem, "reason": "not a large molecule"}
    fitted = bstruct.get("disposition_parameters") or []
    if not fitted:
        return {"skip": True, "stem": stem,
                "reason": "no fitted disposition parameters (prediction/scaling, "
                          "not a parameter-recovery task)"}
    observed = B.biologic_observed(snap, bstruct["molecule"])
    if not observed:
        return {"skip": True, "stem": stem, "reason": "no observed biodistribution"}

    blanked = _blank(snap, {p["name"] for p in fitted})
    agent_input = {
        "schema": "osp-biologic-task/v1",
        "task_id": stem, "source_snapshot": os.path.basename(path),
        "objective": (
            f"Model the disposition of {bstruct['molecule']} (a large molecule) so "
            "it reproduces the observed concentration in whole blood and tissues. "
            "The large-molecule distribution structure is given; recover the "
            "protein-disposition parameter(s) - FcRn recycling affinity and/or the "
            "renal GFR fraction - from the biodistribution."),
        "biologic_structure": {
            "molecule": bstruct["molecule"],
            "molecular_weight": bstruct["molecular_weight"],
            "radius": bstruct["radius"],
            "disposition_note": bstruct["disposition_note"],
            "has_metabolizing_processes": bstruct["has_metabolizing_processes"],
            "observed_matrices": bstruct["observed_matrices"],
        },
        "given_data": {"biodistribution_by_matrix": observed},
        "unknowns_guidance": (
            "This is NOT a small molecule: do not add enzyme or transporter "
            "clearance. The unknowns are the protein-disposition parameter(s) "
            "(FcRn Kd and/or GFR fraction). Estimate them so the model reproduces "
            "the whole-blood and tissue concentrations."),
        "how_scored": ("Per-matrix GMFE (whole blood and each tissue) and an "
                       "overall biodistribution GMFE."),
    }
    answer_edits = {"disposition_parameters": [
        {"compound": bstruct["molecule"], "parameters": fitted}]}
    answer_key = {
        "schema": "osp-biologic-answer-key/v1",
        "warning": "REFERENCE ANSWERS - grading only, never shown to the agent.",
        "molecule": bstruct["molecule"],
        "disposition_note": bstruct["disposition_note"],
        "estimated_parameters": fitted,
        "observed_matrices": bstruct["observed_matrices"],
    }
    return {"skip": False, "stem": stem, "input": agent_input, "blanked": blanked,
            "answer_edits": answer_edits, "answer_key": answer_key,
            "molecule": bstruct["molecule"], "fitted": fitted,
            "n_matrices": len(bstruct["observed_matrices"])}


def _write(base_dir: str, stem: str, res: dict) -> None:
    for sub in ("json_input", "benchmark", "answer_key"):
        os.makedirs(os.path.join(base_dir, sub), exist_ok=True)
    files = {
        os.path.join(base_dir, "json_input", f"{stem}.bio_input.json"): res["input"],
        os.path.join(base_dir, "benchmark", f"{stem}.bio_blanked.json"): res["blanked"],
        os.path.join(base_dir, "answer_key", f"{stem}.bio_answer_edits.json"): res["answer_edits"],
        os.path.join(base_dir, "answer_key", f"{stem}.bio_answer_key.json"): res["answer_key"],
    }
    for p, obj in files.items():
        with open(p, "w", encoding="utf-8") as fh:
            json.dump(obj, fh, indent=2, ensure_ascii=False)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("target", help="a biologic snapshot .json, or a directory")
    ap.add_argument("--write", action="store_true")
    ap.add_argument("--audit", action="store_true")
    args = ap.parse_args()

    files = ([args.target] if os.path.isfile(args.target)
             else sorted(glob.glob(os.path.join(args.target, "*/json/*.json"))))
    for f in files:
        res = build(f)
        stem = res["stem"]
        if res["skip"]:
            if args.audit or os.path.isfile(args.target):
                print(f"SKIP {stem}: {res['reason']}")
            continue
        print(f"{'WROTE' if (args.write and not args.audit) else 'BIO  '} {stem}: "
              f"{res['molecule']}, recover {[p['name'] for p in res['fitted']]}, "
              f"{res['n_matrices']} matrix(es)")
        if args.write and not args.audit:
            _write(os.path.dirname(os.path.dirname(f)), stem, res)


if __name__ == "__main__":
    main()
