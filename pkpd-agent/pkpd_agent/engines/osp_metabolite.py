"""Understand and score a metabolite-cascade PBPK snapshot (parent -> daughter).

A metabolite model is structurally different from a single-compound one:

  * it has >1 Compound, linked by a metabolization process that NAMES its
    product - ``process.Metabolite = "<daughter compound>"`` (e.g. Itraconazole's
    rCYP450_MM on CYP3A4 produces "Hydroxy-Itraconazole");
  * the benchmark quantity is not only the PARENT plasma but each daughter
    METABOLITE's plasma. A good model reproduces the whole cascade, so a parent
    fit that happens to get the metabolite exposure wrong is CAUGHT - which a
    parent-only score cannot do.

This module reads the cascade structure, extracts the observed plasma for each
molecule (parent and daughters alike - the single-compound extractor drops
metabolites on purpose), runs the model once pulling every molecule's plasma from
the same run, and scores each molecule with the shared GMFE metric. Structure
verified against the OSP library Itraconazole snapshot (Itraconazole ->
Hydroxy- -> Keto- -> N-desalkyl-Itraconazole, all via rCYP450_MM on CYP3A4).
"""

from __future__ import annotations

from typing import Any

from . import osp_score

# unit handling (mirrors the observed-data conversion used elsewhere)
_MASS = {"mg/l": 1.0, "µg/l": 1e-3, "ug/l": 1e-3, "µg/ml": 1.0, "ug/ml": 1.0,
         "mg/ml": 1e3, "ng/ml": 1e-3, "ng/l": 1e-6, "g/l": 1e3, "µg/dl": 1e-2,
         "ug/dl": 1e-2}
_TIME = {"h": 1.0, "hr": 1.0, "hour": 1.0, "hours": 1.0, "min": 1 / 60.0,
         "minute": 1 / 60.0, "minutes": 1 / 60.0, "s": 1 / 3600.0,
         "sec": 1 / 3600.0, "day": 24.0, "d": 24.0}


def _to_mg_L(vals, unit, mw):
    u = (unit or "").strip().lower()
    out = []
    for v in vals:
        if not isinstance(v, (int, float)):
            out.append(None)
            continue
        if u in _MASS:
            out.append(v * _MASS[u])
        elif "mol" in u and mw:
            f = mw / 1000.0 * (1e-3 if u.startswith("nmol")
                               else (1e-6 if u.startswith("pmol") else 1.0))
            out.append(v * f)
        else:
            out.append(v)
    return out


def _to_h(vals, unit):
    f = _TIME.get((unit or "").strip().lower(), 1.0)
    return [v * f if isinstance(v, (int, float)) else None for v in vals]


# --------------------------------------------------------------------------- #
# structure
# --------------------------------------------------------------------------- #

def metabolite_edges(snapshot: dict) -> list[dict[str, Any]]:
    """Parent -> daughter edges declared by ``process.Metabolite``."""
    edges = []
    for c in snapshot.get("Compounds") or []:
        for p in c.get("Processes") or []:
            met = p.get("Metabolite")
            if met:
                edges.append({"parent": c.get("Name"),
                              "enzyme": p.get("Molecule"),
                              "internal_name": p.get("InternalName"),
                              "metabolite": met})
    return edges


def _order_chain(root: str | None, edges: list[dict]) -> list[str]:
    """Molecules in cascade order from the root (breadth-first over the edges)."""
    children: dict[str, list[str]] = {}
    for e in edges:
        children.setdefault(e["parent"], []).append(e["metabolite"])
    order, seen = [], set()
    queue = [root] if root else []
    while queue:
        n = queue.pop(0)
        if not n or n in seen:
            continue
        seen.add(n)
        order.append(n)
        queue.extend(children.get(n, []))
    # append any metabolite not reachable from the declared root (defensive)
    for e in edges:
        for m in (e["parent"], e["metabolite"]):
            if m and m not in seen:
                seen.add(m)
                order.append(m)
    return order


def observed_molecules(snapshot: dict) -> list[str]:
    out = set()
    for od in snapshot.get("ObservedData") or []:
        ep = {e.get("Name"): e.get("Value") for e in od.get("ExtendedProperties") or []}
        mol = ep.get("Molecule")
        if mol:
            out.add(mol)
    return sorted(out)


def analyze_metabolites(snapshot: dict) -> dict[str, Any]:
    """Return the metabolite-cascade structure (empty dict if not a cascade)."""
    edges = metabolite_edges(snapshot)
    if not edges:
        return {}
    compounds = [c.get("Name") for c in snapshot.get("Compounds") or []]
    root = compounds[0] if compounds else None
    chain = _order_chain(root, edges)
    obs_mols = observed_molecules(snapshot)
    scorable = [m for m in chain if m in obs_mols]
    metabolites = [m for m in chain if m != root]
    return {
        "is_metabolite": True,
        "root": root,
        "compounds": compounds,
        "metabolites": metabolites,
        "edges": edges,
        "chain": chain,
        "observed_molecules": obs_mols,
        # molecules we can actually grade against data (parent + any daughter
        # with measured plasma); a daughter with no data is a structure we build
        # but cannot score directly.
        "scorable_molecules": scorable,
        "unscorable_metabolites": [m for m in metabolites if m not in obs_mols],
    }


def analyze_multicompound(snapshot: dict) -> dict[str, Any]:
    """Superset of analyze_metabolites: any small-molecule model with >=2
    co-modelled compounds that each carry measured plasma - whether linked by a
    metabolization CASCADE (process.Metabolite edges, e.g. Itraconazole) or run in
    PARALLEL with no cascade edge (e.g. Omeprazole's R/S enantiomers, each with its
    own clearance). Both are scored the same way (per-molecule GMFE + roll-up); the
    only difference is the structure, reported as ``kind``.

    Returns {} for a single-compound model (that is the single-compound task)."""
    compounds = [c for c in snapshot.get("Compounds") or []
                 if c.get("IsSmallMolecule") is not False]
    names = [c.get("Name") for c in compounds]
    if len(names) < 2:
        return {}
    obs_mols = observed_molecules(snapshot)
    scorable = [n for n in names if n in obs_mols]
    if len(scorable) < 2:
        return {}
    edges = metabolite_edges(snapshot)
    root = names[0]
    if edges:
        chain = _order_chain(root, edges)
        # keep only compounds actually present (edges may name the same set)
        chain = [m for m in chain if m in names] or names
        kind = "cascade"
        metabolites = [m for m in chain if m != root]
    else:
        # no cascade edge. A parallel multi-compound fit (e.g. R/S enantiomers) is a
        # metabolite-style task ONLY when it is not really a victim DDI. A
        # perpetrator->victim model (distinct compounds, control/treatment pairs) is
        # the DDI task. Auto-inhibition (every compound inhibits an enzyme that
        # clears itself, e.g. omeprazole on CYP2C19) has no victim and stays here.
        from . import osp_ddi
        d = osp_ddi.analyze_ddi(snapshot)
        if d.get("victims") and d.get("pairs"):
            return {}
        chain = list(names)
        kind = "parallel"
        metabolites = []
    return {
        "is_metabolite": bool(edges),
        "is_multicompound": True,
        "kind": kind,
        "root": root,
        "compounds": names,
        "metabolites": metabolites,
        "edges": edges,
        "chain": chain,
        "observed_molecules": obs_mols,
        "scorable_molecules": [m for m in chain if m in obs_mols],
        "unscorable_metabolites": [m for m in metabolites if m not in obs_mols],
    }


def metabolite_identifiability(mstruct: dict) -> list[dict]:
    """A priori identifiability notes for a cascade, known before any run.

    A daughter metabolite with NO measured plasma cannot pin the fraction of
    parent turned into it vs the rate it is itself cleared - only the product is
    constrained by the parent's disappearance. Surfaced up front, like the DDI
    and single-compound identifiability guidance."""
    acts = []
    for m in mstruct.get("unscorable_metabolites") or []:
        acts.append({
            "molecule": m,
            "issue": (f"{m} has no measured plasma - its formation fraction and "
                      "its own clearance trade off and are not separately "
                      "identifiable from the parent data alone"),
            "action": ("hold this metabolite's disposition at literature/in-vitro "
                       "values; do not free both its formation and elimination"),
            "severity": "medium",
        })
    return acts


# --------------------------------------------------------------------------- #
# observed extraction (per molecule - parent AND daughters)
# --------------------------------------------------------------------------- #

def plasma_observed(snapshot: dict, molecule: str) -> list[dict]:
    """Plasma concentration-time datasets for ONE molecule (parent or daughter).

    Unlike the single-compound extractor - which deliberately drops metabolite
    datasets - this keeps exactly the requested molecule's plasma, so a daughter
    metabolite can be scored."""
    out = []
    for od in snapshot.get("ObservedData") or []:
        ep = {e.get("Name"): e.get("Value") for e in od.get("ExtendedProperties") or []}
        if ep.get("Molecule") != molecule:
            continue
        organ = (ep.get("Organ") or "")
        comp_meas = (ep.get("Compartment") or "")
        is_plasma = ("plasma" in comp_meas.lower()
                     or "venous blood" in organ.lower()
                     or "plasma" in organ.lower())
        is_excretion = any(k in (comp_meas + organ).lower()
                           for k in ("feces", "urine", "bile", "lumen"))
        if is_excretion or (comp_meas and not is_plasma):
            continue
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
        out.append({"dataset": od.get("Name"), "study": ep.get("Study Id"),
                    "molecule": molecule, "route": ep.get("Route"),
                    "dose": ep.get("Dose"), "n_points": n,
                    "time_h": [None if t is None else round(t, 5) for t in time_h[:n]],
                    "conc_mg_L": [None if c is None else round(c, 8) for c in conc[:n]]})
    return out


# --------------------------------------------------------------------------- #
# scoring: score EACH molecule with the shared GMFE metric, then combine
# --------------------------------------------------------------------------- #

def score_metabolites(observed_by_mol: dict[str, list],
                      profiles_by_mol: dict[str, list]) -> dict[str, Any]:
    """Per-molecule fit + a cascade-level roll-up.

    ``observed_by_mol`` : {molecule: [observed dataset dicts]}
    ``profiles_by_mol`` : {molecule: [PredictedProfile]} from ONE run.
    Each molecule is mapped and scored with the same metric the single-compound
    path uses, so a metabolite GMFE is directly comparable to a parent GMFE.
    The cascade GMFE is the geometric mean of per-fold-errors across ALL
    molecules - a parent that fits but whose metabolite is 5x off scores badly."""
    import math
    per_mol, wsum, nsum = {}, 0.0, 0
    for mol, obs in observed_by_mol.items():
        profs = profiles_by_mol.get(mol) or []
        pred, unmatched = osp_score.map_predictions(profs, obs)
        sc = osp_score.score_fit(obs, pred)
        per_mol[mol] = {"overall": sc["overall"], "by_route": sc["by_route"],
                        "n_datasets": len(obs), "n_matched": len(pred),
                        "unmatched": unmatched}
        g, n = sc["overall"].get("gmfe"), sc["overall"].get("n") or 0
        if g and n:                      # weighted geometric mean of molecule GMFEs
            wsum += n * math.log(g)
            nsum += n
    cascade = None
    if nsum:
        cascade = {
            "n": nsum,
            "gmfe": round(math.exp(wsum / nsum), 3),
            "molecules_scored": [m for m, v in per_mol.items()
                                 if v["overall"]["n"]],
        }
    return {"per_molecule": per_mol, "cascade": cascade}


# --------------------------------------------------------------------------- #
# orchestration: run the cascade once and score every molecule
# --------------------------------------------------------------------------- #

def run_metabolite_prediction(cli, snapshot_path: str, mstruct: dict,
                              snapshot: dict | None = None,
                              edits: dict | None = None,
                              simulations: list[str] | None = None) -> dict[str, Any]:
    """Build+run the model ONCE, pull the plasma of the parent and every scorable
    daughter from that single run, and score each molecule. ``snapshot`` (the
    loaded dict) is used to read the observed data; if omitted it is loaded from
    ``snapshot_path``."""
    if snapshot is None:
        import json
        with open(snapshot_path, encoding="utf-8") as fh:
            snapshot = json.load(fh)
    molecules = mstruct.get("scorable_molecules") or []
    if not molecules:
        return {"ok": False, "message": "no scorable molecules (no plasma data "
                "for the parent or any metabolite)"}
    res = cli.build_and_run(snapshot_path, edits=edits, simulations=simulations,
                            prune_simulations=bool(simulations),
                            target_molecules=molecules)
    if not res.get("ok"):
        return {"ok": False, "message": res.get("message"),
                "edits_applied": res.get("edits_applied")}
    profiles_by_mol = res.get("profiles_by_molecule") or {}
    observed_by_mol = {m: plasma_observed(snapshot, m) for m in molecules}
    score = score_metabolites(observed_by_mol, profiles_by_mol)
    out = {"ok": True, "root": mstruct.get("root"),
           "scored_molecules": molecules,
           "n_ran": sum(len(v) for v in profiles_by_mol.values()),
           "score": score}
    recs = metabolite_identifiability(mstruct)
    if recs:
        out["recommendations"] = recs
    return out
