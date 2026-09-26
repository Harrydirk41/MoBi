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

_STUDY_RE = re.compile(r"([A-Za-z]+)[\s_-]*(\d{4}[a-z]?)")   # author (space/_/- ) year


def _norm_study(s: str | None) -> str:
    m = _STUDY_RE.search(s or "")
    if m:
        return (m.group(1) + m.group(2)).lower()
    return re.sub(r"[^a-z0-9]", "", (s or "").lower())


def _study_token(s: str | None) -> str | None:
    """A REAL author+year study id (e.g. 'bovill1984', 'gepts1987' from 'Gepts-1987')
    or None. Unlike _norm_study this does NOT fall back to squashing the whole string,
    so a simulation named only by route/dose ('Digoxin po, 0.5 mg') yields None and
    never falsely constrains."""
    m = _STUDY_RE.search(s or "")
    return (m.group(1) + m.group(2)).lower() if m else None


def _dose_from_name(s: str | None) -> "tuple | None":
    """Dose parsed from a raw simulation/dataset name, tolerant of the COMPACT forms
    OSP uses ('0.5gPer6h', '1gBID', '15ug', '5mg_kg') where the unit is not followed by
    a space. The unit must not run into a lowercase word (so 'g' in 'gemfibrozil' is not
    a gram). Used as a fallback when the structured dose field is absent."""
    if not s:
        return None
    m = re.search(r"(\d+\.?\d*)\s*(µg|ug|mcg|mg|g)(?=$|[^a-z]|[A-Z])", s)
    if not m:
        return None
    unit = m.group(2).lower().replace("mcg", "ug")
    perkg = bool(re.search(r"(µg|ug|mcg|mg|g)\s*[/_]?\s*kg", s, re.IGNORECASE))
    return float(m.group(1)) * _MG[unit], perkg


def _schedule(s: str | None) -> str | None:
    """Dosing schedule from an EXPLICIT token: 'multiple' (repeated/steady-state) or
    'single', else None. A multiple-dose arm accumulates and must not be scored against
    a single-dose simulation - but ABSENCE of a token does NOT imply single (many
    multiple-dose datasets omit it), so it returns None and does not constrain."""
    t = (s or "").lower()
    if re.search(r"\b(t\.?i\.?d|b\.?i\.?d|q\.?d|q\.?i\.?d|multiple[ -]?dose|\bm\.?d\b"
                 r"|steady|day\s*[2-9]|day\s*1[0-9])", t):
        return "multiple"
    if re.search(r"\bs\.?d\b|single[ -]?dose", t):
        return "single"
    return None


def _norm_route(s) -> str | None:
    if not s:
        return None
    u = str(s).upper()
    return "PO" if u in ("PO", "ORAL") else ("IV" if "IV" in u else u)


def _dose_canon(s):
    if s is None or s == "":
        return None
    t = str(s).strip()
    # a bare number with NO unit (e.g. observed dose "300.0") is a dose in mg -
    # the default PK unit. Without this it returned None, which DISABLED the dose
    # constraint and collapsed every arm of a study onto one simulation.
    if re.fullmatch(r"[\d.]+", t):
        try:
            return float(t) * _MG["mg"], False
        except ValueError:
            return None
    m = re.search(r"([\d.]+)\s*(µg|ug|mg|g)\b\s*(/?\s*kg)?", t, re.IGNORECASE)
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


# formulation (tablet/capsule/...) and food state (fed/fasted) are encoded in BOTH
# the observed dataset name ("...8mg.po.sd.Tab_Fed") and the simulation name
# ("Tizanidine 8mg po tablet fed"). Ignoring them lets the scorer pair a FED
# observation with a FASTED simulation (whose dissolution can differ ~10x), which
# silently inflates the reference model's GMFE. Match on them too.
_FORMS = [("tablet", "tab"), ("capsule", "cap"), ("solution", "sol"),
          ("suspension", "susp"), ("granule", "gran"), ("syrup", "syrup"),
          ("pellet", "pellet")]


def _formulation(s: str | None) -> str | None:
    t = (s or "").lower()
    for canon, key in _FORMS:
        if re.search(r"\b" + key, t):
            return canon
    return None


def _food(s: str | None) -> str | None:
    t = (s or "").lower()
    if "fasted" in t or "fasting" in t:
        return "fasted"
    if "fed" in t:
        return "fed"
    return None


def _infusion(s: str | None) -> str | None:
    """IV infusion duration - bolus vs a timed infusion. Different durations at the
    same dose give very different peak shapes, so a 3 h infusion must not be scored
    against a bolus simulation. Returns 'bolus' or minutes (as a string)."""
    t = (s or "").lower()
    if "bolus" in t:
        return "bolus"
    m = re.search(r"(\d+\.?\d*)\s*(h|hr|hour|min)\b", t)
    if not m:
        return None
    val = float(m.group(1))
    mins = val * 60 if m.group(2).startswith("h") else val
    return f"{round(mins)}min"


def _phenotype(s: str | None) -> str | None:
    """Metabolizer phenotype (extensive/poor/intermediate/ultra-rapid) - an EM and a
    PM arm of the same study/dose have very different clearance, so they must not be
    cross-matched."""
    t = (s or "").lower()
    if "poor metaboli" in t:
        return "pm"
    if "extensive metaboli" in t:
        return "em"
    # 'pm' / 'em' as a standalone token (not inside a word like 'system'/'problem')
    m = re.search(r"(?<![a-z])(pm|em)(?![a-z])", t)
    return m.group(1) if m else None


def _match_score(obs: dict, pred) -> int | None:
    """Score a candidate (observed, predicted) pairing; None = incompatible.

    Route and dose are HARD constraints (a mismatch rules the pairing out).
    Study is a soft PREFERENCE: an exact-study match scores higher so it wins the
    best-match search, but a simulation named only by route+dose (e.g.
    'PO SD 50 mg', no author) still matches the right observed arm by route+dose.
    Requiring at least route or dose to positively match prevents a study-less
    or route-less simulation from matching everything."""
    o_study, o_route, o_dose = _obs_key(obs)
    o_name, p_name = obs.get("dataset", ""), getattr(pred, "simulation", "")
    score, hard = 0, False
    # study: HARD when BOTH sides carry a real author+year id (from the field or the
    # raw name). Same study is a strong positive; DIFFERENT studies rule the pairing
    # out - so an observation cannot pile onto another study's simulation just because
    # its route matched. Simulations named only by route/dose yield no study -> skip.
    o_st = _study_token(o_study) or _study_token(o_name)
    p_st = _study_token(pred.study) or _study_token(p_name)
    if o_st and p_st:
        if o_st != p_st:
            return None                              # different study -> out
        score += 3; hard = True
    o_r, p_r = _norm_route(o_route), _norm_route(pred.route)
    if o_r and p_r:
        if o_r != p_r:
            return None                              # route mismatch -> out
        score += 2; hard = True
    o_d = _dose_canon(o_dose) or _dose_from_name(o_name)
    p_d = _dose_canon(pred.dose) or _dose_from_name(p_name)
    if o_d and p_d:
        if o_d[1] != p_d[1] or abs(o_d[0] - p_d[0]) > 1e-6 + 0.01 * o_d[0]:
            return None                              # dose mismatch -> out
        score += 2; hard = True
    # formulation + food state: HARD when BOTH sides name them (like route/dose),
    # so a fed arm never pairs with a fasted simulation of the same study/dose.
    o_form, p_form = _formulation(o_name), _formulation(p_name)
    if o_form and p_form:
        if o_form != p_form:
            return None                              # tablet vs capsule -> out
        score += 1
    o_food, p_food = _food(o_name), _food(p_name)
    if o_food and p_food:
        if o_food != p_food:
            return None                              # fed vs fasted -> out
        score += 1
    o_inf, p_inf = _infusion(o_name), _infusion(p_name)
    if o_inf and p_inf:
        if o_inf != p_inf:
            return None                              # bolus vs 3 h infusion -> out
        score += 1
    o_ph, p_ph = _phenotype(o_name), _phenotype(p_name)
    if o_ph and p_ph:
        if o_ph != p_ph:
            return None                              # EM vs PM -> out
        score += 1
    o_sch, p_sch = _schedule(o_name), _schedule(p_name)
    if o_sch and p_sch:
        if o_sch != p_sch:
            return None                              # single-dose vs multiple-dose -> out
        score += 1
    # need a positive study/route/dose match; a bare formulation/food is too weak.
    return score if hard else None


def linkage_from_snapshot(snap: dict) -> dict[str, str]:
    """The EXPLICIT observed-dataset -> simulation pairing OSP records in each
    simulation's OutputMappings/ObservedData. This is the authoritative linkage the
    modeller set up (e.g. 'Shah2006...Tab_Fed' -> 'Tizanidine 8mg po tablet fed'), so
    scoring never has to re-guess it from names. A dataset the reference maps to NO
    simulation (an extra arm the model was not built/validated on - e.g. Digoxin's DDI
    control arms) is simply absent, and is then excluded from scoring rather than
    force-matched. Returns {observed_dataset_name: simulation_name}."""
    link: dict[str, str] = {}
    for s in snap.get("Simulations") or []:
        nm = s.get("Name")
        for om in s.get("OutputMappings") or []:
            od = om.get("ObservedData")
            if od and nm:
                link.setdefault(od, nm)
        for od in s.get("ObservedData") or []:
            if isinstance(od, str) and nm:
                link.setdefault(od, nm)
    return link


def map_predictions(profiles: list, observed: list[dict], linkage: dict | None = None):
    """-> (predicted_profiles[list of dicts], unmatched[list of dataset names]).

    When ``linkage`` (from ``linkage_from_snapshot``) is given, each observation is
    paired with its simulation EXACTLY as OSP recorded it, and observations the model
    maps to no simulation are left unmatched (excluded), not heuristically force-fit.
    Without a linkage it falls back to name-based ``_match_score`` matching."""
    predicted_profiles, unmatched = [], []
    prof_by_sim = {p.simulation: p for p in profiles} if linkage else {}
    for o in observed:
        ds = o["dataset"]
        best = None
        if linkage is not None:
            sim = linkage.get(ds)
            best = prof_by_sim.get(sim) if sim else None    # unmapped or no-profile -> None
        else:
            best_score = 0
            for p in profiles:
                s = _match_score(o, p)
                if s and s > best_score:
                    best, best_score = p, s
        if best:
            predicted_profiles.append({
                "dataset": ds, "time_h": best.time_h,
                "pred_conc_mg_L": best.conc_mg_L,
                "_from_simulation": best.simulation})
        else:
            unmatched.append(ds)
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
        "pk_parameters": pk_parameters(observed, predicted_profiles),
    }


# --------------------------------------------------------------------------- #
# PK-parameter goodness-of-fit (Cmax / tmax / AUC / t1/2), the classic PBPK
# evaluation table (Kuepfer 2016): GMFE captures the whole curve, but a modeler
# also checks whether the PEAK, its timing, the EXPOSURE and the TERMINAL slope
# each land. Observed and predicted are evaluated on the SAME (observed) time
# grid so the comparison is apples-to-apples and independent of grid density.
# --------------------------------------------------------------------------- #

def _pk_params(times: list[float], concs: list[float]) -> "dict | None":
    pts = [(t, c) for t, c in zip(times, concs)
           if _finite(t) is not None and _finite(c) is not None and c > 0]
    pts.sort()
    if len(pts) < 2:
        return None
    cmax = max(c for _, c in pts)
    tmax = next(t for t, c in pts if c == cmax)
    auc = sum((pts[i][0] - pts[i - 1][0]) * (pts[i][1] + pts[i - 1][1]) / 2.0
              for i in range(1, len(pts)))                 # linear trapezoid
    # terminal half-life: log-linear slope over the points from tmax onward
    tail = [(t, c) for t, c in pts if t >= tmax and c > 0]
    thalf = None
    if len(tail) >= 3:
        xs = [t for t, _ in tail]
        ys = [math.log(c) for _, c in tail]
        n = len(xs)
        mx, my = sum(xs) / n, sum(ys) / n
        den = sum((x - mx) ** 2 for x in xs)
        if den > 0:
            slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / den
            if slope < 0:
                thalf = math.log(2) / (-slope)
    return {"cmax": cmax, "tmax": tmax, "auc": auc, "thalf": thalf}


def _fold(a: "float | None", r: "float | None") -> "float | None":
    if a is None or r is None or r == 0:
        return None
    return round(a / r, 3)


def pk_parameters(observed: list[dict], predicted_profiles: list[dict]) -> dict[str, Any]:
    """Per-dataset observed vs predicted Cmax/tmax/AUC/t1/2 with fold-errors, plus
    the median |log| fold-error per metric across datasets (GMFE-style)."""
    preds = {p["dataset"]: p for p in predicted_profiles}
    rows, agg = [], {"cmax": [], "auc": [], "tmax": [], "thalf": []}
    for o in observed:
        pr = preds.get(o["dataset"])
        if not pr:
            continue
        ot = list(o["time_h"])
        obs = _pk_params(ot, list(o["conc_mg_L"]))
        # predicted sampled at the SAME observed times
        psamp = [_interp(pr["time_h"], pr["pred_conc_mg_L"], t) for t in ot]
        sim = _pk_params(ot, [p if p is not None else 0.0 for p in psamp])
        if not obs or not sim:
            continue
        row = {"dataset": o["dataset"], "study": o.get("study"),
               "route": (_norm_route(_obs_key(o)[1]) or "NA")}
        for k in ("cmax", "tmax", "auc", "thalf"):
            f = _fold(sim.get(k), obs.get(k))
            row[k + "_obs"] = (round(obs.get(k), 4) if obs.get(k) is not None else None)
            row[k + "_pred"] = (round(sim.get(k), 4) if sim.get(k) is not None else None)
            row[k + "_fold"] = f
            if f and f > 0:
                agg[k].append(abs(math.log(f)))
        rows.append(row)
    summary = {k: (round(math.exp(sum(v) / len(v)), 3) if v else None)
               for k, v in agg.items()}   # geometric-mean fold-error per metric
    return {"per_dataset": rows, "gmfe_by_metric": summary}


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


def plausibility(params: list[dict], hard_only: bool = False) -> list[dict]:
    """Flag implausible parameter values. With ``hard_only`` (minimal/purist mode)
    only PHYSICAL IMPOSSIBILITY is flagged — fraction unbound outside (0, 1], a
    non-positive clearance/permeability, a GFR fraction outside [0, 1]. The softer
    biological-plausibility opinions (a lipophilicity range, clearance vs hepatic
    blood flow) are dropped: judging those is the modeling skill under test, so the
    model makes that call itself."""
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
            if not hard_only and not -2 <= v <= 7:      # soft opinion
                bad = "lipophilicity expected within [-2, 7]"
        elif "clearance" in name or re.search(r"\bcl\b", name):
            if v <= 0:
                bad = "clearance must be > 0"
            elif not hard_only and "l/min" in unit and v > 1.5:   # soft opinion
                bad = "clearance > ~1.5 L/min exceeds adult hepatic blood flow"
        elif "permeab" in name and v <= 0:
            bad = "permeability must be > 0"
        elif "gfr" in name:
            # GFR fraction is a physiological fraction (unbound drug filtered at GFR).
            # Outside [0, 1] is genuinely unphysical -> flag. A value below 1 is NOT
            # flagged: a GFR fraction < 1 is common and legitimate (renal clearance
            # below the passive filtration ceiling), so the harness offers no opinion
            # on where within [0, 1] it should sit - that is a modeling judgement.
            if not 0 <= v <= 1:
                bad = "GFR fraction is a physiological fraction, must be in [0, 1]"
        if bad:
            flags.append({"parameter": p.get("parameter"), "value": v, "message": bad})
    return flags


def _downsample_curve(pts: list, cap: int = 80) -> list:
    """Thin a simulated curve to ~cap evenly-spaced points for a compact overlay
    (the raw PK-Sim grid is far denser than a plot needs)."""
    n = len(pts)
    if n <= cap:
        return pts
    step = n / float(cap)
    out = [pts[int(i * step)] for i in range(cap)]
    if out[-1] is not pts[-1]:
        out.append(pts[-1])                       # always keep the terminal point
    return out


def overlay_series(observed: list[dict], predicted: list[dict],
                   cap: int = 80) -> list[dict]:
    """Per-dataset observed vs simulated curves for a fit overlay (fit-vision).

    ``observed`` entries carry ``time_h``/``conc_mg_L``; ``predicted`` entries
    (from :func:`map_predictions`) carry ``time_h``/``pred_conc_mg_L`` keyed by
    ``dataset``. Returns ``[{study, route, dose, observed:[[t,c]],
    simulated:[[t,c]]}]`` - the shape :mod:`fit_plot` renders."""
    pred_by = {p.get("dataset"): p for p in (predicted or [])}
    out = []
    for o in observed or []:
        p = pred_by.get(o.get("dataset"))
        obs_pts = [[t, c] for t, c in zip(o.get("time_h", []), o.get("conc_mg_L", []))
                   if c is not None]
        sim_pts = ([[t, c] for t, c in zip(p.get("time_h", []), p.get("pred_conc_mg_L", []))
                    if c is not None] if p else [])
        out.append({"study": o.get("study") or o.get("dataset"),
                    "route": o.get("route", ""), "dose": o.get("dose", ""),
                    "observed": obs_pts, "simulated": _downsample_curve(sim_pts, cap)})
    return out
