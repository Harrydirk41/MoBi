"""Understand a multi-compound DDI snapshot (perpetrator -> victim interaction).

A drug-drug-interaction PBPK model is structurally different from a single-
compound one:

  * it contains >1 Compound - a PERPETRATOR (carries inhibition/induction
    processes on an enzyme/transporter) and a VICTIM (cleared by that enzyme);
  * each treatment Simulation carries an ``Interactions`` block
    ``[{Name, MoleculeName, CompoundName}]`` that switches the perpetrator's
    effect on into that run;
  * CONTROL simulations (victim alone, no Interactions) and TREATMENT
    simulations (victim + perpetrator) come in pairs, and the benchmark quantity
    is the DDI RATIO - how much the victim's exposure changes (AUC_treatment /
    AUC_control), not the absolute concentration.

This module reads the snapshot and returns that structure so the agent can be
handed a DDI task instead of a single-compound fit. Structures verified against
the OSP library Erythromycin (Midazolam victim) snapshot.
"""

from __future__ import annotations

import re
from typing import Any

from .osp_cli import OSPCli
from .osp_catalog import INTERACTION_PROCESS_TYPES

# InternalNames that mark a compound as a perpetrator (inhibition / induction)
_INTERACTION_INTERNAL_NAMES = {s["internal_name"]
                               for s in INTERACTION_PROCESS_TYPES.values()}


def _perpetrator_mechanisms(compound: dict) -> list[dict[str, Any]]:
    out = []
    for p in compound.get("Processes") or []:
        inm = p.get("InternalName")
        if inm in _INTERACTION_INTERNAL_NAMES:
            out.append({
                "internal_name": inm,
                "target": p.get("Molecule") or p.get("MoleculeName"),
                "data_source": p.get("DataSource"),
                "parameters": {par.get("Name"): par.get("Value")
                               for par in p.get("Parameters") or []},
            })
    return out


def _pair_key(sim_name: str) -> tuple[str, str | None]:
    """(study, route) key so a control and ITS treatment arm match while unrelated
    arms do not. The study is the author+year token (e.g. 'olkkola1993'); the
    route folds Oral->PO. A control (victim alone) and its treatment (victim +
    perpetrator) share the study author+year; a perpetrator's own PK arm from a
    different study has a different author+year, so it will not pair with the
    victim's control."""
    _, route, _ = OSPCli._parse_sim_name(sim_name)
    m = re.search(r"([A-Za-z]{3,})[ _]?((?:19|20)\d{2})", sim_name)
    study = (m.group(1) + m.group(2)).lower() if m else \
        re.sub(r"(?i)(control|treatment|baseline|alone|_|\W)", "", sim_name).lower()
    return study, route


def analyze_ddi(snapshot: dict) -> dict[str, Any]:
    """Return the DDI structure of a snapshot (empty dict if not a DDI model)."""
    compounds = snapshot.get("Compounds") or []
    sims = snapshot.get("Simulations") or []

    perpetrators, victims = [], []
    for c in compounds:
        mechs = _perpetrator_mechanisms(c)
        if mechs:
            perpetrators.append({"name": c.get("Name"), "mechanisms": mechs})
        else:
            victims.append(c.get("Name"))

    treatment_sims = [s for s in sims if s.get("Interactions")]
    control_sims = [s for s in sims if not s.get("Interactions")]
    if not treatment_sims or not perpetrators:
        return {}

    # unique interaction links used across treatment simulations
    seen, interactions = set(), []
    for s in treatment_sims:
        for it in s.get("Interactions") or []:
            key = (it.get("Name"), it.get("MoleculeName"), it.get("CompoundName"))
            if key not in seen:
                seen.add(key)
                interactions.append({"name": it.get("Name"),
                                     "molecule": it.get("MoleculeName"),
                                     "perpetrator": it.get("CompoundName")})

    # pair each control (victim alone) with its treatment arm (victim +
    # perpetrator) on shared study author+year AND route. Treatment arms with no
    # matching control (e.g. the perpetrator's own PK) are left UNPAIRED - they
    # are not victim-DDI pairs.
    ctrl_by_key = {}
    for s in control_sims:
        ctrl_by_key.setdefault(_pair_key(s["Name"]), s["Name"])
    pairs, unpaired = [], []
    for s in treatment_sims:
        ctrl = ctrl_by_key.get(_pair_key(s["Name"]))
        if ctrl:
            pairs.append({"control": ctrl, "treatment": s["Name"],
                          "route": _pair_key(s["Name"])[1]})
        else:
            unpaired.append(s["Name"])

    return {
        "is_ddi": True,
        "compounds": [c.get("Name") for c in compounds],
        "perpetrators": perpetrators,
        "victims": victims,
        "n_control": len(control_sims),
        "n_treatment": len(treatment_sims),
        "control_simulations": [s["Name"] for s in control_sims],
        "treatment_simulations": [s["Name"] for s in treatment_sims],
        "interactions": interactions,
        "pairs": pairs,
        "unpaired_treatments": unpaired,
    }


# --------------------------------------------------------------------------- #
# DDI metric: the interaction RATIO (victim exposure with vs without perpetrator)
# --------------------------------------------------------------------------- #
# The benchmark quantity for a victim DDI is not the absolute concentration but
# how much the perpetrator changes the victim's exposure: AUCR = AUC_treatment /
# AUC_control (and analogously Cmax ratio). A good DDI model reproduces that
# ratio - e.g. erythromycin roughly quadruples midazolam AUC.

def _auc(time_h: list, conc: list) -> float:
    pts = sorted((t, c) for t, c in zip(time_h, conc)
                 if isinstance(t, (int, float)) and isinstance(c, (int, float))
                 and c is not None and c >= 0)
    return sum((pts[i][0] - pts[i-1][0]) * (pts[i][1] + pts[i-1][1]) / 2.0
               for i in range(1, len(pts)))


def _cmax(conc: list) -> float:
    vals = [c for c in conc if isinstance(c, (int, float)) and c is not None]
    return max(vals) if vals else 0.0


def interaction_ratios(pairs: list, profiles_by_sim: dict) -> list[dict[str, Any]]:
    """For each control/treatment pair, the victim AUC and Cmax ratios.

    ``profiles_by_sim``: {sim_name: {"time_h": [...], "conc": [...]}} for the
    VICTIM compound (the measured plasma), control and treatment alike."""
    out = []
    for pr in pairs:
        c, t = profiles_by_sim.get(pr["control"]), profiles_by_sim.get(pr["treatment"])
        if not c or not t:
            continue
        ac, at = _auc(c["time_h"], c["conc"]), _auc(t["time_h"], t["conc"])
        cc, ct = _cmax(c["conc"]), _cmax(t["conc"])
        out.append({
            "treatment": pr["treatment"], "control": pr["control"],
            "route": pr.get("route"),
            "auc_ratio": (at / ac) if ac > 0 else None,
            "cmax_ratio": (ct / cc) if cc > 0 else None})
    return out


def score_ddi(predicted: list[dict], observed: list[dict]) -> dict[str, Any]:
    """Compare predicted vs observed interaction ratios. observed/predicted are
    lists of {treatment, auc_ratio, cmax_ratio}. Reports per-arm fold error and
    the geometric-mean fold error (GMFE) on the AUC ratio - the DDI accuracy."""
    obs = {o["treatment"]: o for o in observed}
    import math
    rows, logfolds = [], []
    for p in predicted:
        o = obs.get(p["treatment"])
        pr, orr = p.get("auc_ratio"), (o or {}).get("auc_ratio")
        fold = None
        if isinstance(pr, (int, float)) and isinstance(orr, (int, float)) \
                and pr > 0 and orr > 0:
            fold = pr / orr
            logfolds.append(abs(math.log(fold)))
        rows.append({"treatment": p["treatment"], "route": p.get("route"),
                     "predicted_aucr": round(pr, 3) if isinstance(pr, (int, float)) else None,
                     "observed_aucr": round(orr, 3) if isinstance(orr, (int, float)) else None,
                     "fold_error": round(fold, 3) if fold else None})
    gmfe = round(math.exp(sum(logfolds) / len(logfolds)), 4) if logfolds else None
    within2 = (round(100.0 * sum(1 for lf in logfolds if lf <= math.log(2))
                     / len(logfolds), 1) if logfolds else None)
    return {"gmfe_aucr": gmfe, "within_2fold_pct": within2, "per_arm": rows,
            "n": len(logfolds)}


# --------------------------------------------------------------------------- #
# orchestration: run the control/treatment pairs and predict the DDI ratio
# --------------------------------------------------------------------------- #

def run_ddi_prediction(cli, snapshot_path: str, ddi: dict, victim: str,
                       edits: dict | None = None,
                       observed_ratios: list | None = None) -> dict[str, Any]:
    """Build+run the control and treatment simulations, extract the VICTIM's
    plasma from each, and compute the predicted interaction ratios (optionally
    scored against observed). ``edits`` may tune the perpetrator's interaction
    parameters (interaction_parameters). This is the DDI analogue of the
    single-compound build-and-score path; the model structure (which enzyme, the
    mechanism) is the ground-truth model's - the agent tunes/predicts the
    interaction, then checks the ratio."""
    pairs = ddi.get("pairs") or []
    sims = sorted({p["control"] for p in pairs} | {p["treatment"] for p in pairs})
    if not sims:
        return {"ok": False, "message": "no control/treatment pairs to run"}
    res = cli.build_and_run(snapshot_path, edits=edits, simulations=sims,
                            prune_simulations=True, target_molecule=victim)
    if not res.get("ok"):
        return {"ok": False, "message": res.get("message"), "edits_applied": res.get("edits_applied")}
    profiles_by_sim = {p.simulation: {"time_h": p.time_h, "conc": p.conc_mg_L}
                       for p in res["profiles"]}
    predicted = interaction_ratios(pairs, profiles_by_sim)
    out = {"ok": True, "victim": victim, "predicted_ratios": predicted,
           "n_pairs": len(pairs), "n_ran": len(profiles_by_sim)}
    if observed_ratios:
        out["score"] = score_ddi(predicted, observed_ratios)
    return out
