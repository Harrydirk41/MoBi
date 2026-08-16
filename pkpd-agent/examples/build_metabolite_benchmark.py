r"""Generate a metabolite-cascade benchmark from an OSP snapshot.

The metabolite analogue of build_benchmark.py / build_ddi_benchmark.py. A
cascade model's ANSWER is the parent AND daughter disposition parameters that
were identified, and - in HARD mode - the cascade STRUCTURE itself (which
metabolization step produces which metabolite). The benchmark quantity is the
plasma of the PARENT and every measured METABOLITE, so a parent-only fit that
gets a metabolite exposure wrong is caught.

Writes (normal):
    <Compound>/json_input/<stem>.metab_input.json        agent task
    <Compound>/benchmark/<stem>.metab_blanked.json       value-blanked start
    <Compound>/answer_key/<stem>.metab_answer_edits.json fitted params (cascade)
    <Compound>/answer_key/<stem>.metab_answer_key.json   + cascade structure

Writes (--hard): the STRUCTURE is stripped too - the metabolization processes
that name each product are removed, so the agent must rebuild the cascade from
the biology.
    <Compound>/json_input/<stem>.metab.hard.input.json
    <Compound>/benchmark/<stem>.metab.hard_blanked.json

    python -m examples.build_metabolite_benchmark ^
        ..\OSP-PBPK-Model-Library\Itraconazole\json\Itraconazole-Model.json --write
    python -m examples.build_metabolite_benchmark --audit ..\OSP-PBPK-Model-Library
"""

from __future__ import annotations

import argparse
import copy
import glob
import json
import os

from pkpd_agent.engines import osp_metabolite as M

# no-leak naive priors, shared with the single-compound generator's spirit
BLANK_DEFAULTS = {"Lipophilicity": 2.0, "Solubility at reference pH": 1.0,
                  "Specific intestinal permeability (transcellular)": 1e-6,
                  "Specific clearance": 1.0, "Ki": 1.0, "kcat": 1.0,
                  "Fraction unbound (plasma, reference value)": 0.1,
                  "GFR fraction": 1.0}


def _blank_default(name: str) -> float:
    if name in BLANK_DEFAULTS:
        return BLANK_DEFAULTS[name]
    n = name.lower()
    if "lipophil" in n:
        return 2.0
    if "clearance" in n or "kcat" in n or "vmax" in n:
        return 1.0
    if "unbound" in n:
        return 0.1
    if "permeab" in n:
        return 1e-6
    return 1.0


def _fitted(compound: dict) -> list[tuple[str, float, str]]:
    out = []

    def w(o):
        if isinstance(o, dict):
            nm = o.get("Name")
            if isinstance(nm, str) and isinstance(o.get("Value"), (int, float)) \
                    and (o.get("ValueOrigin") or {}).get("Source") == "ParameterIdentification":
                out.append((nm, o["Value"], o.get("Unit", "")))
            for v in o.values():
                w(v)
        elif isinstance(o, list):
            for v in o:
                w(v)

    w(compound)
    return out


def fitted_cascade_params(snapshot: dict) -> list[dict]:
    """The identified disposition parameters, per compound in the cascade."""
    specs = []
    for c in snapshot.get("Compounds") or []:
        fit = _fitted(c)
        if fit:
            specs.append({"compound": c.get("Name"),
                          "parameters": [{"name": n, "value": v, "unit": u}
                                         for (n, v, u) in fit]})
    return specs


def _blank_values(snapshot: dict) -> dict:
    """Reset every ParameterIdentification value (across all cascade compounds)
    to a no-leak naive prior."""
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


def _strip_cascade_structure(blanked: dict) -> tuple[dict, list[dict]]:
    """HARD mode: remove the metabolization processes that NAME each product, so
    the agent must rebuild the cascade from biology (which enzyme turns the parent
    into which metabolite). Keeps the daughter compounds and their physchem - the
    unknown is the topology + kinetics, not the existence of the metabolites."""
    snap = copy.deepcopy(blanked)
    removed = []
    for c in snap.get("Compounds") or []:
        kept = []
        for p in c.get("Processes") or []:
            if p.get("Metabolite"):
                removed.append({"parent": c.get("Name"),
                                "enzyme": p.get("Molecule"),
                                "internal": p.get("InternalName"),
                                "metabolite": p.get("Metabolite")})
            else:
                kept.append(p)
        c["Processes"] = kept
    # drop the mirrored process refs in simulations (no dangling references)
    for sim in snap.get("Simulations") or []:
        for sc in sim.get("Compounds") or []:
            if sc.get("Processes"):
                sc["Processes"] = [p for p in sc["Processes"]
                                   if not p.get("Metabolite")]
    return snap, removed


def cascade_biology_facts(mstruct: dict) -> list[str]:
    """The metabolic pathway stated as biology - the ONLY structure source in
    hard mode. E.g. 'Itraconazole is metabolized by CYP3A4 to Hydroxy-...'."""
    facts = []
    for e in mstruct.get("edges") or []:
        facts.append(f"{e['parent']} is metabolized by {e['enzyme']} to "
                     f"{e['metabolite']}.")
    scor = mstruct.get("scorable_molecules") or []
    if scor:
        facts.append("Measured plasma is available for: " + ", ".join(scor) + ".")
    return facts


def build(path: str) -> dict:
    with open(path, encoding="utf-8") as fh:
        snap = json.load(fh)
    stem = os.path.splitext(os.path.basename(path))[0]
    mstruct = M.analyze_metabolites(snap)
    if not mstruct:
        return {"skip": True, "stem": stem, "reason": "not a metabolite cascade "
                "(no process names a product Metabolite)"}
    specs = fitted_cascade_params(snap)
    if not specs:
        return {"skip": True, "stem": stem,
                "reason": "no fitted parameters to recover"}
    scorable = mstruct.get("scorable_molecules") or []
    if not scorable:
        return {"skip": True, "stem": stem,
                "reason": "no measured plasma for parent or any metabolite"}

    obs_by_mol = {m: M.plasma_observed(snap, m) for m in scorable}
    blanked = _blank_values(snap)
    hard_blanked, removed = _strip_cascade_structure(blanked)

    root = mstruct["root"]
    agent_input = {
        "schema": "osp-metabolite-task/v1",
        "task_id": stem, "source_snapshot": os.path.basename(path),
        "objective": (
            f"Build a PBPK model for {root} AND its metabolite cascade that "
            f"reproduces the observed plasma of {root} and its measured "
            "metabolites. Decide the disposition parameters you cannot read from "
            "the data and justify them. The model must fit the whole cascade, not "
            "just the parent."),
        "cascade_structure": {
            "parent": root,
            "chain": mstruct["chain"],
            "metabolites": mstruct["metabolites"],
            "metabolic_steps": [{"parent": e["parent"], "enzyme": e["enzyme"],
                                 "metabolite": e["metabolite"]}
                                for e in mstruct["edges"]],
            "scorable_molecules": scorable,
        },
        "background": {"metabolic_pathway": cascade_biology_facts(mstruct)},
        "given_data": {"plasma_by_molecule": obs_by_mol},
        "unknowns_guidance": (
            "The metabolization kinetics and each compound's disposition are the "
            "unknowns. Fit them so the model reproduces the parent AND metabolite "
            "plasma. A metabolite with sparse/absent data is structurally weak - "
            "hold its disposition at literature values rather than fitting both "
            "its formation and its clearance."),
        "how_scored": ("Per-molecule GMFE on plasma (parent and each measured "
                       "metabolite) and a cascade GMFE across all molecules."),
    }
    answer_edits = {"cascade_parameters": specs}
    answer_key = {
        "schema": "osp-metabolite-answer-key/v1",
        "warning": "REFERENCE ANSWERS - grading only, never shown to the agent.",
        "parent": root,
        "cascade_structure": [{"parent": e["parent"], "type": e["internal_name"],
                               "enzyme": e["enzyme"], "metabolite": e["metabolite"]}
                              for e in mstruct["edges"]],
        "estimated_parameters": specs,
    }
    hard_input = copy.deepcopy(agent_input)
    hard_input["mode"] = "hard"
    hard_input["unknowns_guidance"] = (
        "HARD mode - the cascade STRUCTURE is NOT in the snapshot: the "
        "metabolization processes that turn each compound into the next have been "
        "removed. Using the metabolic pathway in the background, ADD the "
        "metabolization processes (each producing the correct metabolite) on the "
        "expressed enzymes, then fit the kinetics so the model reproduces the "
        "parent and metabolite plasma.")
    # in hard mode the cascade topology is the answer, so it must not be echoed
    # back as a given structure - only the biology facts remain.
    hard_input.pop("cascade_structure", None)

    return {"skip": False, "stem": stem, "input": agent_input, "blanked": blanked,
            "hard_input": hard_input, "hard_blanked": hard_blanked,
            "removed_processes": removed,
            "answer_edits": answer_edits, "answer_key": answer_key,
            "n_molecules": len(scorable), "specs": specs,
            "scorable": scorable, "chain": mstruct["chain"]}


def _write(base_dir: str, stem: str, res: dict, hard: bool) -> None:
    subs = ("json_input", "benchmark") if hard else ("json_input", "benchmark", "answer_key")
    for sub in subs:
        os.makedirs(os.path.join(base_dir, sub), exist_ok=True)
    if hard:
        files = {
            os.path.join(base_dir, "json_input", f"{stem}.metab.hard.input.json"): res["hard_input"],
            os.path.join(base_dir, "benchmark", f"{stem}.metab.hard_blanked.json"): res["hard_blanked"],
        }
    else:
        files = {
            os.path.join(base_dir, "json_input", f"{stem}.metab_input.json"): res["input"],
            os.path.join(base_dir, "benchmark", f"{stem}.metab_blanked.json"): res["blanked"],
            os.path.join(base_dir, "answer_key", f"{stem}.metab_answer_edits.json"): res["answer_edits"],
            os.path.join(base_dir, "answer_key", f"{stem}.metab_answer_key.json"): res["answer_key"],
        }
    for p, obj in files.items():
        with open(p, "w", encoding="utf-8") as fh:
            json.dump(obj, fh, indent=2, ensure_ascii=False)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("target", help="a metabolite snapshot .json, or a directory")
    ap.add_argument("--write", action="store_true")
    ap.add_argument("--hard", action="store_true",
                    help="write ONLY the hard-mode artifacts (cascade structure removed)")
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
        tag = "[HARD]" if args.hard else ""
        print(f"{'WROTE' if (args.write and not args.audit) else 'METAB'} {stem} {tag}: "
              f"cascade {' -> '.join(res['chain'])}; {res['n_molecules']} scorable "
              f"molecule(s), {sum(len(s['parameters']) for s in res['specs'])} params")
        if args.hard:
            for r in res["removed_processes"]:
                print(f"      removed: {r['parent']} --{r['enzyme']}--> {r['metabolite']}")
        if args.write and not args.audit:
            _write(os.path.dirname(os.path.dirname(f)), stem, res, args.hard)


if __name__ == "__main__":
    main()
