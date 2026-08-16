"""Understand and score a large-molecule (protein / mAb) PBPK snapshot.

A biologic model is a different class from a small molecule:

  * ``IsSmallMolecule = False`` - the disposition is size-limited protein
    distribution (governed by the hydrodynamic Radius), FcRn-mediated recycling
    (``Kd (FcRn) in endosomal space`` drives the long half-life), and renal
    filtration (``GFR fraction``) - NOT enzymes or transporters. These models
    carry no metabolization processes at all.
  * the data are a BIODISTRIBUTION: concentration-time in whole blood AND many
    tissues, not a single plasma curve. The benchmark quantity is the fit across
    every measured matrix (each organ/compartment), so a model that matches blood
    but not the tissues is caught.

This module reads the biologic structure, extracts observed concentration per
(organ, compartment) matrix (which the small-molecule plasma extractor drops -
it is tissue, not plasma), runs the model once pulling every matrix from the same
run, and scores each matrix with the shared GMFE metric. Structures verified
against the OSP library BAY794620 (mAb, 17 tissue + whole-blood matrices, recovers
Kd(FcRn)), dAb2 (recovers GFR fraction), and Tefibazumab (plasma, recovers Kd).
"""

from __future__ import annotations

import math
from typing import Any

from . import osp_score

# disposition parameters that ARE the biologic unknown (not enzyme kinetics)
_BIOLOGIC_PARAMS = ("Kd (FcRn) in endosomal space", "GFR fraction",
                    "Radius (solute)", "Kd", "koff", "kass (FcRn)",
                    "Reference concentration")

_MASS = {"mg/l": 1.0, "µg/l": 1e-3, "ug/l": 1e-3, "µg/ml": 1.0, "ug/ml": 1.0,
         "mg/ml": 1e3, "ng/ml": 1e-3, "ng/l": 1e-6, "g/l": 1e3,
         "µg/g": 1.0, "ug/g": 1.0, "% id/g": None}
_TIME = {"h": 1.0, "hr": 1.0, "hour": 1.0, "hours": 1.0, "min": 1 / 60.0,
         "minute": 1 / 60.0, "minutes": 1 / 60.0, "s": 1 / 3600.0,
         "day": 24.0, "days": 24.0, "d": 24.0, "week": 168.0, "weeks": 168.0}


def _to_mg_L(vals, unit, mw):
    u = (unit or "").strip().lower()
    out = []
    for v in vals:
        if not isinstance(v, (int, float)):
            out.append(None)
            continue
        if u in _MASS and _MASS[u] is not None:
            out.append(v * _MASS[u])
        elif "mol" in u and mw:
            f = mw / 1000.0 * (1e-3 if u.startswith("nmol")
                               else (1e-6 if u.startswith("pmol") else 1.0))
            out.append(v * f)
        else:
            out.append(v)     # unknown unit (e.g. %ID/g): keep raw, scored in-kind
    return out


def _to_h(vals, unit):
    f = _TIME.get((unit or "").strip().lower(), 1.0)
    return [v * f if isinstance(v, (int, float)) else None for v in vals]


# --------------------------------------------------------------------------- #
# structure
# --------------------------------------------------------------------------- #

def _fitted(compound: dict) -> list[dict]:
    out = []

    def w(o):
        if isinstance(o, dict):
            nm = o.get("Name")
            if isinstance(nm, str) and isinstance(o.get("Value"), (int, float)) \
                    and (o.get("ValueOrigin") or {}).get("Source") == "ParameterIdentification":
                out.append({"name": nm, "value": o["Value"], "unit": o.get("Unit", "")})
            for v in o.values():
                w(v)
        elif isinstance(o, list):
            for v in o:
                w(v)

    w(compound)
    return out


def observed_matrices(snapshot: dict, molecule: str) -> list[tuple[str, str]]:
    seen = []
    for od in snapshot.get("ObservedData") or []:
        ep = {e.get("Name"): e.get("Value") for e in od.get("ExtendedProperties") or []}
        if ep.get("Molecule") != molecule:
            continue
        key = (ep.get("Organ") or "", ep.get("Compartment") or "")
        if key not in seen:
            seen.append(key)
    return seen


def analyze_biologic(snapshot: dict) -> dict[str, Any]:
    """Return the biologic structure (empty dict if not a large molecule)."""
    compounds = snapshot.get("Compounds") or []
    large = [c for c in compounds if c.get("IsSmallMolecule") is False]
    if not large:
        return {}
    primary = large[0]
    name = primary.get("Name")
    fitted = _fitted(primary)
    # physchem the modeler is given
    physchem = {}
    for p in primary.get("Parameters") or []:
        if p.get("Name") in ("Molecular weight", "Radius (solute)"):
            physchem[p["Name"]] = {"value": p.get("Value"), "unit": p.get("Unit", "")}
    matrices = observed_matrices(snapshot, name)
    has_enzyme_process = any(p.get("Molecule") for c in compounds
                             for p in c.get("Processes") or [])
    return {
        "is_biologic": True,
        "molecule": name,
        "compounds": [c.get("Name") for c in compounds],
        "molecular_weight": physchem.get("Molecular weight"),
        "radius": physchem.get("Radius (solute)"),
        "disposition_parameters": fitted,     # the unknown(s) to recover
        "observed_matrices": [{"organ": o, "compartment": c} for (o, c) in matrices],
        "has_metabolizing_processes": has_enzyme_process,
        "disposition_note": ("large-molecule disposition: size-limited tissue "
                             "distribution + FcRn recycling + renal filtration; "
                             "no enzymatic clearance"),
    }


def biologic_identifiability(bstruct: dict) -> list[dict]:
    """A priori note: with only whole-blood data an FcRn Kd and the renal GFR
    fraction both change the terminal slope and are hard to separate; tissue data
    resolves distribution but not necessarily recycling vs filtration."""
    acts = []
    params = {p["name"] for p in bstruct.get("disposition_parameters") or []}
    matrices = bstruct.get("observed_matrices") or []
    only_blood = all(("blood" in (m["organ"] or "").lower()
                      or "plasma" in (m["compartment"] or "").lower())
                     for m in matrices) if matrices else True
    if {"Kd (FcRn) in endosomal space", "GFR fraction"} <= params and only_blood:
        acts.append({
            "issue": "FcRn recycling (Kd) and renal filtration (GFR fraction) both "
                     "shape the terminal decline; blood-only data cannot fully "
                     "separate them",
            "action": "fix the better-known one (GFR fraction from protein size) "
                      "and estimate the other",
            "severity": "medium"})
    return acts


# --------------------------------------------------------------------------- #
# observed extraction (per matrix: whole blood + each tissue)
# --------------------------------------------------------------------------- #

def biologic_observed(snapshot: dict, molecule: str) -> list[dict]:
    """Concentration-time datasets for the biologic across ALL matrices - whole
    blood, plasma, and each tissue - keyed by (organ, compartment). Unlike the
    small-molecule extractor (plasma only), tissue matrices are kept: a mAb model
    is graded on its biodistribution."""
    out = []
    for od in snapshot.get("ObservedData") or []:
        ep = {e.get("Name"): e.get("Value") for e in od.get("ExtendedProperties") or []}
        if ep.get("Molecule") != molecule:
            continue
        organ = ep.get("Organ") or ""
        compartment = ep.get("Compartment") or ""
        base = od.get("BaseGrid") or {}
        time_h = _to_h(list(base.get("Values") or []), base.get("Unit", ""))
        cols = od.get("Columns") or []
        if not cols or not time_h:
            continue
        col = cols[0]
        mw = (col.get("DataInfo") or {}).get("MolWeight")
        conc = _to_mg_L(list(col.get("Values") or []), col.get("Unit", ""), mw)
        n = min(len(time_h), len(conc))
        if n < 2:
            continue
        key = f"{organ}|{compartment}"
        out.append({"dataset": od.get("Name") or key, "matrix": key,
                    "organ": organ, "compartment": compartment,
                    "molecule": molecule, "study": ep.get("Study Id"),
                    "n_points": n,
                    "time_h": [None if t is None else round(t, 5) for t in time_h[:n]],
                    "conc_mg_L": [None if c is None else round(c, 8) for c in conc[:n]]})
    return out


def matrix_specs(observed: list[dict], molecule: str) -> list[dict]:
    """Column selectors for build_and_run(target_matrices=...): one per distinct
    matrix. Carries molecule + organ + compartment so the CLI matches the organ as
    an exact pipe segment (with aliases) - a substring organ name cannot match a
    different organ's column."""
    specs, seen = [], set()
    for o in observed:
        key = o["matrix"]
        if key in seen:
            continue
        seen.add(key)
        specs.append({"key": key, "molecule": molecule,
                      "organ": o.get("organ"), "compartment": o.get("compartment")})
    return specs


# --------------------------------------------------------------------------- #
# scoring: per matrix, then an overall biodistribution GMFE
# --------------------------------------------------------------------------- #

def score_biologic(observed: list[dict],
                   profiles_by_matrix: dict[str, list]) -> dict[str, Any]:
    """Score each matrix by interpolating the predicted profile onto the observed
    times (fold error per point), then a biodistribution-wide GMFE. Matching is by
    matrix key (organ|compartment), not study/route - a mAb biodistribution is
    one dose, many matrices."""
    obs_by_matrix: dict[str, list[dict]] = {}
    for o in observed:
        obs_by_matrix.setdefault(o["matrix"], []).append(o)
    per_matrix, all_fe = {}, []
    for key, obs_list in obs_by_matrix.items():
        profs = profiles_by_matrix.get(key) or []
        pred = profs[0] if profs else None
        fe = []
        if pred is not None:
            for o in obs_list:
                for t, c in zip(o["time_h"], o["conc_mg_L"]):
                    of = osp_score._finite(c)
                    tf = osp_score._finite(t)
                    if of is None or of <= 0 or tf is None:
                        continue
                    p = osp_score._interp(pred.time_h, pred.conc_mg_L, tf)
                    if p is None or p <= 0:
                        continue
                    fe.append(p / of)
        per_matrix[key] = {**osp_score._metrics(fe), "matched": pred is not None}
        all_fe.extend(fe)
    overall = osp_score._metrics(all_fe)
    return {"overall": overall, "per_matrix": per_matrix,
            "n_matrices": len(obs_by_matrix),
            "n_matched": sum(1 for v in per_matrix.values() if v["matched"])}


# --------------------------------------------------------------------------- #
# orchestration
# --------------------------------------------------------------------------- #

def run_biologic_prediction(cli, snapshot_path: str, bstruct: dict,
                            snapshot: dict | None = None,
                            edits: dict | None = None,
                            simulations: list[str] | None = None) -> dict[str, Any]:
    """Build+run the biologic ONCE, pull every measured matrix from that run, and
    score the biodistribution."""
    if snapshot is None:
        import json
        with open(snapshot_path, encoding="utf-8") as fh:
            snapshot = json.load(fh)
    molecule = bstruct.get("molecule")
    observed = biologic_observed(snapshot, molecule)
    if not observed:
        return {"ok": False, "message": f"no observed matrices for {molecule}"}
    specs = matrix_specs(observed, molecule)
    res = cli.build_and_run(snapshot_path, edits=edits, simulations=simulations,
                            prune_simulations=bool(simulations),
                            target_matrices=specs)
    if not res.get("ok"):
        return {"ok": False, "message": res.get("message"),
                "edits_applied": res.get("edits_applied")}
    by_matrix = res.get("profiles_by_matrix") or {}
    score = score_biologic(observed, by_matrix)
    out = {"ok": True, "molecule": molecule, "n_matrices": len(specs),
           "n_matched": score["n_matched"], "score": score}
    recs = biologic_identifiability(bstruct)
    if recs:
        out["recommendations"] = recs
    return out
