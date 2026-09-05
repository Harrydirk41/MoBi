"""Map PK-Sim simulations to observed datasets and score the fit.

Shared by the standalone runner (examples/osp_run) and the LLM loop tools so
they agree on mapping and metrics. Beyond GMFE and %-within-2-fold this also
reports the geometric-mean BIAS per route (>1 = the model over-predicts), which
is the signal the agent needs to know *which way* to move a parameter:

    IV over-predicts  -> clearance too low        -> raise intrinsic clearance
    PO under-predicts -> absorption too low        -> raise intestinal permeability
"""

from __future__ import annotations

import math
import re
from typing import Any

from .osp_cli import OSPCli

_MG = {"µg": 1e-3, "ug": 1e-3, "mg": 1.0, "g": 1e3}


# --------------------------------------------------------------------------- #
# simulation <-> observed matching (study, route, dose)
# --------------------------------------------------------------------------- #

def _norm_study(s: str | None) -> str:
    m = re.search(r"([A-Za-z]+)\s*(\d{4}[a-z]?)", s or "")
    if m:
        return (m.group(1) + m.group(2)).lower()
    return re.sub(r"[^a-z0-9]", "", (s or "").lower())


def _norm_route(s) -> str | None:
    if not s:
        return None
    u = str(s).upper()
    return "PO" if u in ("PO", "ORAL") else ("IV" if "IV" in u else u)


def _dose_canon(s):
    if not s:
        return None
    m = re.search(r"([\d.]+)\s*(µg|ug|mg|g)\b\s*(/?\s*kg)?", str(s), re.IGNORECASE)
    if not m:
        return None
    return float(m.group(1)) * _MG[m.group(2).lower()], bool(m.group(3))


def _obs_key(o: dict):
    study, route, dose = o.get("study"), o.get("route"), o.get("dose")
    if not (study and route):
        ps, pr, pd = OSPCli._parse_sim_name(o.get("dataset", ""))
        study = study or ps
        route = route or pr
        dose = dose or pd
    return study, route, dose


def _match_score(obs: dict, pred) -> int | None:
    """Score a candidate (observed, predicted) pairing; None = incompatible.

    Route and dose are HARD constraints (a mismatch rules the pairing out).
    Study is a soft PREFERENCE: an exact-study match scores higher so it wins the
    best-match search, but a simulation named only by route+dose (e.g.
    'PO SD 50 mg', no author) still matches the right observed arm by route+dose.
    Requiring at least route or dose to positively match prevents a study-less
    or route-less simulation from matching everything."""
    o_study, o_route, o_dose = _obs_key(obs)
    score, hard = 0, False
    if _norm_study(o_study) and _norm_study(pred.study) \
            and _norm_study(o_study) == _norm_study(pred.study):
        score += 1                                   # study bonus (soft)
    o_r, p_r = _norm_route(o_route), _norm_route(pred.route)
    if o_r and p_r:
        if o_r != p_r:
            return None                              # route mismatch -> out
        score += 2; hard = True
    o_d, p_d = _dose_canon(o_dose), _dose_canon(pred.dose)
    if o_d and p_d:
        if o_d[1] != p_d[1] or abs(o_d[0] - p_d[0]) > 1e-6 + 0.01 * o_d[0]:
            return None                              # dose mismatch -> out
        score += 2; hard = True
    # need a positive route/dose match; study alone is too weak to map on.
    return score if hard else None


def map_predictions(profiles: list, observed: list[dict]):
    """-> (predicted_profiles[list of dicts], unmatched[list of dataset names])."""
    predicted_profiles, unmatched = [], []
    for o in observed:
        best, best_score = None, 0
        for p in profiles:
            s = _match_score(o, p)
            if s and s > best_score:
                best, best_score = p, s
        if best:
            predicted_profiles.append({
                "dataset": o["dataset"], "time_h": best.time_h,
                "pred_conc_mg_L": best.conc_mg_L,
                "_from_simulation": best.simulation})
        else:
            unmatched.append(o["dataset"])
    return predicted_profiles, unmatched


# --------------------------------------------------------------------------- #
# metrics
# --------------------------------------------------------------------------- #

def _finite(x):
    try:
        f = float(x)
    except (TypeError, ValueError):
        return None
    return None if f != f else f


def _interp(xp, fp, x):
    pts = sorted((a, b) for a, b in zip(xp, fp)
                 if _finite(a) is not None and _finite(b) is not None)
    if len(pts) < 2 or x < pts[0][0] or x > pts[-1][0]:
        for a, b in pts:
            if abs(a - x) < 1e-9:
                return b
        return None
    for i in range(1, len(pts)):
        x0, y0 = pts[i - 1]
        x1, y1 = pts[i]
        if x0 <= x <= x1:
            return y0 if x1 == x0 else y0 + (x - x0) / (x1 - x0) * (y1 - y0)
    return None


def _metrics(fe: list[float]) -> dict[str, Any]:
    if not fe:
        return {"n": 0, "gmfe": None, "within_2fold_pct": None, "bias": None}
    logs = [math.log(f) for f in fe]
    return {
        "n": len(fe),
        "gmfe": round(math.exp(sum(abs(l) for l in logs) / len(logs)), 3),
        "within_2fold_pct": round(100 * sum(1 for f in fe if 0.5 <= f <= 2) / len(fe), 1),
        "bias": round(math.exp(sum(logs) / len(logs)), 3),   # >1 = over-predicts
    }


def score_fit(observed: list[dict], predicted_profiles: list[dict]) -> dict[str, Any]:
    preds = {p["dataset"]: p for p in predicted_profiles}
    obs_by = {o["dataset"]: o for o in observed}
    all_fe, by_route, per_ds = [], {}, []
    for name, o in obs_by.items():
        pr = preds.get(name)
        if not pr:
            continue
        route = (_norm_route(_obs_key(o)[1]) or "NA")
        fe = []
        for t, c in zip(o["time_h"], o["conc_mg_L"]):
            of, tf = _finite(c), _finite(t)
            if of is None or of <= 0 or tf is None:
                continue
            p = _interp(pr["time_h"], pr["pred_conc_mg_L"], tf)
            if p is None or p <= 0:
                continue
            fe.append(p / of)
        if fe:
            all_fe.extend(fe)
            by_route.setdefault(route, []).extend(fe)
            m = _metrics(fe)
            per_ds.append({"dataset": name, "study": o.get("study"),
                           "route": route, **m})
    return {
        "overall": _metrics(all_fe),
        "by_route": {r: _metrics(fe) for r, fe in by_route.items()},
        "per_dataset": sorted(per_ds, key=lambda d: (d["gmfe"] or 0), reverse=True),
    }


# --------------------------------------------------------------------------- #
# parameter plausibility (compact; shares the rubric's bounds)
# --------------------------------------------------------------------------- #

def physical_bounds(name: str) -> "tuple[float, float] | None":
    """The GENERAL physical plausibility range for a parameter, by kind (not drug-specific). Used to
    widen a FREE (non-given) fit-target that railed to a too-tight self-imposed bound - e.g. the
    effective Lipophilicity that PK-Sim fits for distribution, which can sit far from a measured logP.
    Returns None for a parameter with no defined physical range."""
    n = (name or "").lower()
    if "lipophil" in n or n in ("logp", "logd"):
        return (-2.0, 7.0)                         # same range the plausibility guard enforces
    if "unbound" in n or n in ("fu", "fup"):
        return (0.0, 1.0)
    if "gfr" in n:
        return (0.0, 1.0)
    return None


def plausibility(params: list[dict]) -> list[dict]:
    flags = []
    for p in params or []:
        name = (p.get("parameter") or "").lower()
        v = _finite(p.get("value"))
        unit = (p.get("unit") or "").lower()
        if v is None:
            continue
        bad = None
        if "unbound" in name or name in ("fu", "fup"):
            if not 0 < v <= 1:
                bad = "fraction unbound must be in (0, 1]"
        elif "lipophil" in name or name in ("logp", "logd"):
            if not -2 <= v <= 7:
                bad = "lipophilicity expected within [-2, 7]"
        elif "clearance" in name or re.search(r"\bcl\b", name):
            if v <= 0:
                bad = "clearance must be > 0"
            elif "l/min" in unit and v > 1.5:
                bad = "clearance > ~1.5 L/min exceeds adult hepatic blood flow"
        elif "permeab" in name and v <= 0:
            bad = "permeability must be > 0"
        elif "gfr" in name:
            # GFR fraction is a physiological fraction (unbound drug filtered at GFR); outside [0,1]
            # is unphysical, and a value pushed well below 1 to rescue a fit implies net tubular
            # reabsorption that needs a mechanistic justification, not a free knob.
            if not 0 <= v <= 1:
                bad = "GFR fraction is a physiological fraction, must be in [0, 1]"
            elif v < 0.5:
                bad = ("GFR fraction < 0.5 implies substantial net tubular reabsorption - justify it "
                       "mechanistically or check the distribution method before fitting it this low")
        if bad:
            flags.append({"parameter": p.get("parameter"), "value": v, "message": bad})
    return flags
