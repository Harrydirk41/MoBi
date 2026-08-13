"""Build an agent-facing INPUT file (and a separate answer key) from an OSP snapshot.

An OSP (PK-Sim / MoBi) snapshot mixes three things: the INPUTS a modeler was
handed, the MODELING STRATEGY they chose, and the RESULTS they obtained. For a
blind modeling benchmark the agent must see only the inputs, be told how it will
be graded and what to produce, and never see the reference answers.

This script emits two files per snapshot:

  json_input/<stem>.input.json   -- AGENT-FACING. Objective, public background,
      literature physicochemical data (ranges, no fitted values), the clinical
      observed data to reproduce, study designs, demographics, a generic hint
      about unknowns, an output/submission template, and an evaluation rubric.
      No fitted values, no method choices, no reference metrics.

  answer_key/<stem>.answer_key.json  -- GRADING ONLY. The parameters the
      reference model estimated (with values), the distribution/permeability
      method choices, the metabolizing processes, and any reference metrics
      (e.g. GMFE) carried in the problem card. NOT shown to the agent.

Public, curated per-compound context (objective + literature physchem + known
biology) comes from a problem-cards file (see problem_cards.json), keyed by
snapshot stem, so no compound-specific text is invented by the generic code.

Usage:
    python build_clean_input.py json/Alfentanil-Model.json
    python build_clean_input.py json/Alfentanil-Pediatrics.json
    # options: --cards problem_cards.json  --outdir json_input  --keydir answer_key
"""

from __future__ import annotations

import argparse
import json
import os
from typing import Any


# --------------------------------------------------------------------------- #
# unit helpers
# --------------------------------------------------------------------------- #

def _num(v: Any) -> float | None:
    """Coerce to float; missing markers ('', 'NaN', None, lists) -> None."""
    if isinstance(v, bool) or isinstance(v, (list, dict)) or v is None:
        return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return None if f != f else f                          # drop NaN


def _to_mg_L(values: list, unit: str, mol_weight: float | None) -> list:
    """Convert a concentration column to mg/L (non-numeric entries -> None)."""
    nums = [_num(v) for v in values]
    u = (unit or "").strip().lower()
    mass = {"mg/l": 1.0, "µg/l": 1e-3, "ug/l": 1e-3, "ng/ml": 1e-3,
            "ng/l": 1e-6, "g/l": 1e3}
    if u in mass:
        f = mass[u]
    elif "mol" in u and mol_weight:                      # molar -> mass
        f = mol_weight / 1000.0                           # µmol/L * g/mol /1000 = mg/L
        if u.startswith("nmol"):
            f /= 1000.0
    else:
        return nums                                       # unknown unit: pass through
    return [None if v is None else v * f for v in nums]


def _time_to_h(values: list, unit: str) -> list:
    nums = [_num(v) for v in values]
    u = (unit or "").strip().lower()
    if u.startswith("min"):
        return [None if v is None else v / 60.0 for v in nums]
    if u.startswith("day"):
        return [None if v is None else v * 24.0 for v in nums]
    return nums


def _ext_props(od: dict[str, Any]) -> dict[str, Any]:
    return {e.get("Name"): e.get("Value") for e in od.get("ExtendedProperties") or []}


# --------------------------------------------------------------------------- #
# INPUT extractors (safe to show the agent)
# --------------------------------------------------------------------------- #

def clinical_observed(data: dict[str, Any]) -> list[dict[str, Any]]:
    """Tidy every observed profile: metadata + time(h) + conc(mg/L) + SD."""
    out = []
    for od in data.get("ObservedData") or []:
        base = od.get("BaseGrid") or {}
        time_h = _time_to_h(list(base.get("Values") or []), base.get("Unit", ""))
        cols = od.get("Columns") or []
        if not cols or not time_h:
            continue
        col = cols[0]
        mw = (col.get("DataInfo") or {}).get("MolWeight")
        conc_native = list(col.get("Values") or [])
        conc_mg_L = _to_mg_L(conc_native, col.get("Unit", ""), mw)

        sd_mg_L = None
        for rc in col.get("RelatedColumns") or []:
            if (rc.get("DataInfo") or {}).get("AuxiliaryType") == "ArithmeticStdDev":
                sd_mg_L = _to_mg_L(list(rc.get("Values") or []), rc.get("Unit", ""), mw)
                break

        ep = _ext_props(od)
        n = min(len(time_h), len(conc_native))
        rec = {
            "dataset": od.get("Name"),
            "study": ep.get("Study Id"),
            "molecule": ep.get("Molecule"),
            "route": ep.get("Route"),
            "dose": ep.get("Dose"),
            "matrix": "/".join(x for x in (ep.get("Organ"), ep.get("Compartment")) if x),
            "n_points": n,
            "time_unit": "h",
            "conc_unit": "mg/L",
            "mol_weight_g_per_mol": mw,
            "time_h": [None if t is None else round(t, 5) for t in time_h[:n]],
            "conc_mg_L": [None if c is None else round(c, 8) for c in conc_mg_L[:n]],
        }
        if sd_mg_L is not None:
            rec["sd_mg_L"] = [None if s is None else round(s, 8) for s in sd_mg_L[:n]]
        out.append(rec)
    return out


def study_designs(data: dict[str, Any]) -> list[dict[str, Any]]:
    out = []
    for pr in data.get("Protocols") or []:
        params = {p.get("Name"): (p.get("Value"), p.get("Unit"))
                  for p in pr.get("Parameters") or []}
        dose, dunit = params.get("InputDose", (None, None))
        start, sunit = params.get("Start time", (None, None))
        out.append({
            "name": pr.get("Name"),
            "application_type": pr.get("ApplicationType"),
            "dosing_interval": pr.get("DosingInterval"),
            "dose": dose, "dose_unit": dunit,
            "start_time": start, "start_time_unit": sunit,
        })
    return out


def demographics(data: dict[str, Any]) -> list[dict[str, Any]]:
    out = []
    for ind in data.get("Individuals") or []:
        od = ind.get("OriginData") or {}
        age = od.get("Age") or {}
        wt = od.get("Weight") or {}
        out.append({
            "name": ind.get("Name"), "kind": "individual",
            "species": od.get("Species"), "population": od.get("Population"),
            "gender": od.get("Gender"),
            "age": age.get("Value"), "age_unit": age.get("Unit"),
            "weight": wt.get("Value"), "weight_unit": wt.get("Unit"),
        })
    for pop in data.get("Populations") or []:
        settings = pop.get("Settings") or {}
        out.append({
            "name": pop.get("Name"), "kind": "population",
            "number_of_individuals": settings.get("NumberOfIndividuals")
            or pop.get("NumberOfIndividuals"),
        })
    return out


# --------------------------------------------------------------------------- #
# ANSWER-KEY extractors (must NOT be shown to the agent)
# --------------------------------------------------------------------------- #

def answer_key(data: dict[str, Any], card: dict[str, Any]) -> dict[str, Any]:
    comp = (data.get("Compounds") or [{}])[0]
    fitted: list[dict[str, Any]] = []
    seen: set[tuple] = set()

    def walk(obj: Any) -> None:
        if isinstance(obj, dict):
            if "Name" in obj and isinstance(obj.get("Value"), (int, float)):
                vo = obj.get("ValueOrigin") or {}
                if vo.get("Source") == "ParameterIdentification":
                    key = (obj["Name"], obj.get("Unit"))
                    if key not in seen:
                        seen.add(key)
                        fitted.append({"parameter": obj["Name"], "value": obj["Value"],
                                       "unit": obj.get("Unit", "")})
            for v in obj.values():
                walk(v)
        elif isinstance(obj, list):
            for v in obj:
                walk(v)

    walk(comp)

    methods = [m if isinstance(m, str) else m.get("Name")
               for m in comp.get("CalculationMethods") or []]
    processes = [{"molecule": p.get("Molecule"), "type": p.get("InternalName")}
                 for p in comp.get("Processes") or [] if p.get("Molecule")]

    return {
        "schema": "osp-answer-key/v1",
        "warning": "REFERENCE ANSWERS - do NOT provide to the agent. Grading only.",
        "compound": comp.get("Name"),
        "structural_choices": {
            "calculation_methods": methods,
            "metabolizing_processes": processes,
        },
        "estimated_parameters": fitted,
        "reference_metrics": card.get("reference_metrics", {}),
        "notes": card.get("answer_notes",
                          "Reference goodness-of-fit (GMFE) and diagnostic plots "
                          "are in Alfentanil_evaluation_report.md sections 3.1-3.3."),
    }


# --------------------------------------------------------------------------- #
# fixed agent-facing scaffolding
# --------------------------------------------------------------------------- #

def submission_template() -> dict[str, Any]:
    return {
        "instructions": (
            "This is what you must produce. Return an object with this exact "
            "shape (or add it back into this file under a top-level 'submission' "
            "key). predicted_profiles MUST align to the datasets and time points "
            "in clinical_observed_data (same dataset names, same time_h grid). "
            "List EVERY parameter you fixed or estimated beyond the given data, "
            "with units and a short rationale."),
        "submission": {
            "structural_model": {
                "distribution_model": None,
                "absorption_model": None,
                "elimination_pathways": [],
                "engine": None,
                "notes": None,
            },
            "parameters": [
                {"parameter": None, "value": None, "unit": None,
                 "fixed_or_estimated": None, "rationale": None},
            ],
            "predicted_profiles": [
                {"dataset": None, "time_h": [], "pred_conc_mg_L": []},
            ],
            "self_assessment": {
                "data_fit_gmfe_overall": None,
                "pct_within_2fold": None,
                "parameter_plausibility_notes": None,
                "output_plausibility_notes": None,
            },
        },
    }


def evaluation_rubric() -> dict[str, Any]:
    return {
        "how_scored": (
            "Your submission is graded on three independent dimensions. You are "
            "NOT shown the reference answers; aim for a model that is both "
            "well-fitting and physically sensible."),
        "dimensions": [
            {
                "id": "data_fit",
                "goal": "Predicted plasma concentrations match the observed data.",
                "metrics": [
                    {"name": "GMFE",
                     "definition": "geometric mean fold error = exp(mean(|ln(pred/obs)|)) over paired points",
                     "target": "<= 2 (good <= 1.5)"},
                    {"name": "pct_within_2fold",
                     "definition": "fraction of predicted points within 2x of observed",
                     "target": "majority; higher is better"},
                ],
                "report_by": ["overall", "route (IV vs PO)", "study"],
            },
            {
                "id": "parameter_plausibility",
                "goal": "Every parameter you set is physically/physiologically sensible.",
                "checks": [
                    "fraction unbound in (0, 1]",
                    "lipophilicity (logP/logD) within about [-2, 7]",
                    "clearance does not exceed the relevant organ blood flow (e.g. hepatic ~1.5 L/min in an adult)",
                    "permeabilities > 0 and within physiological range (cm/min)",
                    "volume of distribution physiologically bounded (~0.05-50 L/kg)",
                    "all rate constants and doses positive; units consistent",
                ],
            },
            {
                "id": "output_plausibility",
                "goal": "Simulated profiles behave like real pharmacokinetics.",
                "checks": [
                    "concentrations >= 0 at all times; no numerical blow-up",
                    "IV bolus declines monotonically from t0; oral rises to a Cmax then declines",
                    "terminal half-life consistent (within ~2x) across dose levels unless nonlinearity is justified",
                    "dose-normalized AUC roughly constant across doses, or nonlinearity explicitly explained",
                    "Cmax and AUC increase with dose",
                ],
            },
        ],
        "self_assessment_required": True,
    }


UNKNOWNS_GUIDANCE = (
    "Some compound parameters cannot be read off the given literature/in-vitro "
    "data and must be chosen or estimated so the model reproduces the observed "
    "data. Decide which parameters to fix (from the given values) and which to "
    "estimate, choose a distribution/permeability approach and an absorption "
    "model, and justify each choice. The reference model's parameter set, method "
    "choices, and results are withheld.")


# --------------------------------------------------------------------------- #
# build
# --------------------------------------------------------------------------- #

def build(path: str, card: dict[str, Any]) -> tuple[dict, dict]:
    with open(path, encoding="utf-8") as fh:
        data = json.load(fh)
    stem = os.path.splitext(os.path.basename(path))[0]

    comp = (data.get("Compounds") or [{}])[0]
    obs = clinical_observed(data)
    routes = sorted({o["route"] for o in obs if o.get("route")})
    studies = sorted({o["study"] for o in obs if o.get("study")})

    agent_input = {
        "schema": "osp-agent-task/v2",
        "task_id": card.get("task_id", stem),
        "source_snapshot": os.path.basename(path),
        "compound": comp.get("Name"),
        "objective": card.get("objective",
                              f"Build a PBPK model for {comp.get('Name')} that "
                              f"reproduces the observed plasma concentration-time data."),
        "background": card.get("background", {}),
        "given_data": {
            "compound_identity": {
                "name": comp.get("Name"),
                "is_small_molecule": comp.get("IsSmallMolecule"),
            },
            "literature_physicochemical": card.get("literature_physicochemical", []),
            "clinical_observed_data": obs,
            "study_designs": study_designs(data),
            "demographics": demographics(data),
        },
        "data_overview": {
            "n_observed_datasets": len(obs),
            "routes": routes,
            "studies": studies,
            "n_study_designs": len(data.get("Protocols") or []),
            "n_demographics": len(data.get("Individuals") or [])
            + len(data.get("Populations") or []),
        },
        "unknowns_guidance": UNKNOWNS_GUIDANCE,
        "what_you_must_produce": submission_template(),
        "evaluation_rubric": evaluation_rubric(),
    }
    return agent_input, answer_key(data, card)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("snapshot", help="path to the OSP snapshot .json")
    ap.add_argument("--cards", default="problem_cards.json",
                    help="problem-cards JSON keyed by snapshot stem (default: problem_cards.json)")
    ap.add_argument("--outdir", default="json_input", help="agent-input folder")
    ap.add_argument("--keydir", default="answer_key", help="answer-key folder")
    args = ap.parse_args()

    stem = os.path.splitext(os.path.basename(args.snapshot))[0]
    cards = {}
    if args.cards and os.path.exists(args.cards):
        with open(args.cards, encoding="utf-8") as fh:
            cards = json.load(fh)
    card = cards.get(stem, {})
    if not card:
        print(f"warning: no problem card for '{stem}' in {args.cards}; "
              f"input will have no curated objective/background/literature.")

    agent_input, key = build(args.snapshot, card)

    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.keydir, exist_ok=True)
    in_path = os.path.join(args.outdir, f"{stem}.input.json")
    key_path = os.path.join(args.keydir, f"{stem}.answer_key.json")
    with open(in_path, "w", encoding="utf-8") as fh:
        json.dump(agent_input, fh, indent=2, ensure_ascii=False)
    with open(key_path, "w", encoding="utf-8") as fh:
        json.dump(key, fh, indent=2, ensure_ascii=False)

    ov = agent_input["data_overview"]
    print(f"AGENT INPUT -> {in_path}")
    print(f"  compound          : {agent_input['compound']}")
    print(f"  observed datasets : {ov['n_observed_datasets']}  routes={ov['routes']}")
    print(f"  given physchem    : {len(agent_input['given_data']['literature_physicochemical'])} literature entries")
    print(f"  leaks             : none (no fitted values / methods / metrics)")
    print(f"ANSWER KEY  -> {key_path}")
    print(f"  estimated params  : {len(key['estimated_parameters'])} (values hidden from agent)")
    print(f"  method choices    : {key['structural_choices']['calculation_methods']}")


if __name__ == "__main__":
    main()
