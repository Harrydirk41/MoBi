"""Model-agnostic engine for the downstream QSP tasks (trial design, calibration,
Vpop generation, drug design, held-out validation).

This is the 2-7 analogue of qsp_core.py: pure-Python primitives with NO model
vocabulary. A run's CSV columns, the drug/parameter catalogs and the clinical
reference numbers all arrive as arguments (the caller reads them off a
QSPTaskConfig, which is loaded from projects/<name>/tasks.json). Nothing here
names a drug, a cytokine or a trial.

Two families of function:
  * run summarizers -- turn a run_vpop result into population response RATES, given
    a semantic column map (which CSV column plays which role).
  * scorers / spec builders / numerical routines -- score a predicted response,
    assemble dose / override / sample strings, weight a Vpop to target moments, run
    a bounded 1-D fit. These never touch column names.
"""

from __future__ import annotations

from typing import Any, Optional

_INF = (float("inf"), float("-inf"))


def _finite(xs) -> list[float]:
    return [x for x in xs
            if isinstance(x, (int, float)) and x == x and x not in _INF]


def _rate(flags) -> Optional[float]:
    """Percent of finite flags set to 1 (a model-set responder flag)."""
    xs = _finite(flags)
    if not xs:
        return None
    return round(100.0 * sum(1 for x in xs if x >= 0.5) / len(xs), 1)


def _mean(xs) -> Optional[float]:
    xs = _finite(xs)
    return round(sum(xs) / len(xs), 3) if xs else None


# ----------------------------------------------------------- run summarizers -- #
def summarize_run(res: dict, columns: dict) -> dict[str, Any]:
    """A run_vpop result -> {first_line, second_line, severity} population rates.

    ``columns`` is the semantic column map (from the task config):
      {patient, first_line{role->col}, subgroup_flag, second_line{role->col},
       severity{baseline, readout}}
    first_line rates are over all patients; second_line rates are restricted to the
    subgroup where ``subgroup_flag`` == 1 (e.g. the inadequate responders), so the
    caller never hardcodes which flag or which columns a given model emits.
    """
    cols = res.get("columns") or {}
    pcol = columns.get("patient", "patient")
    n = len(_finite(cols.get(pcol, [])))

    first = {"n": n}
    for role, colname in (columns.get("first_line") or {}).items():
        first[role] = _rate(cols.get(colname, []))

    flag = columns.get("subgroup_flag")
    nonresp = cols.get(flag, []) if flag else []

    def _isset(v):
        return isinstance(v, (int, float)) and v == v and v >= 0.5

    n_ir = sum(1 for v in nonresp if _isset(v))

    def _among_subgroup(colname):
        col = cols.get(colname, [])
        vals = [v for v, nr in zip(col, nonresp)
                if isinstance(v, (int, float)) and v == v and _isset(nr)]
        if not vals:
            return None
        return round(100.0 * sum(1 for v in vals if v >= 0.5) / len(vals), 1)

    second = {"n_subgroup": n_ir}
    for role, colname in (columns.get("second_line") or {}).items():
        second[role] = _among_subgroup(colname)

    dcols = columns.get("severity") or {}
    sev = {"baseline_mean": _mean(cols.get(dcols.get("baseline", ""), [])),
           "readout_mean": _mean(cols.get(dcols.get("readout", ""), []))}
    return {"first_line": first, "second_line": second, "severity": sev}


def ir_mask(run: dict, columns: dict, acr_key: str = None,
            threshold: float = 3.2) -> dict:
    """Classify inadequate responders from a run: a patient is an IR if they did NOT
    reach the response flag (== 0) AND still have active disease (severity readout >
    threshold). Returns {patient_id: bool}, keyed by patient so arms align.

    ``columns`` supplies patient / first_line role->col / severity.readout; ``acr_key``
    is a first_line ROLE (e.g. 'ACR50') resolved through the column map, defaulting
    to the second first_line role if omitted.
    """
    cols = run.get("columns") or {}
    fl = columns.get("first_line") or {}
    roles = list(fl.keys())
    role = acr_key if acr_key in fl else (roles[1] if len(roles) > 1 else roles[0])
    pcol = columns.get("patient", "patient")
    acol = fl[role]
    dcol = (columns.get("severity") or {}).get("readout", "")
    pats = cols.get(pcol, [])
    acr = cols.get(acol, [])
    das = cols.get(dcol, [])
    m: dict[int, bool] = {}
    for p, a, d in zip(pats, acr, das):
        if not isinstance(p, (int, float)):
            continue
        if not (isinstance(a, (int, float)) and a == a):
            continue
        if not (isinstance(d, (int, float)) and d == d):
            continue
        m[int(p)] = (a < 0.5) and (d > threshold)
    return m


def response_in_subgroup(run: dict, ids: set, columns: dict,
                         roles=None) -> dict:
    """Percent responders on each first_line role, restricted to the patient ids in
    the subgroup. ``roles`` selects which first_line roles to report (default all)."""
    cols = run.get("columns") or {}
    fl = columns.get("first_line") or {}
    pcol = columns.get("patient", "patient")
    pats = cols.get(pcol, [])
    out: dict[str, Any] = {"n": len(ids)}
    for role in (roles or list(fl.keys())):
        colname = fl.get(role, role)
        col = cols.get(colname, [])
        vals = [v for p, v in zip(pats, col)
                if isinstance(p, (int, float)) and int(p) in ids
                and isinstance(v, (int, float)) and v == v]
        out[role] = round(100.0 * sum(1 for v in vals if v >= 0.5) / len(vals), 1) \
            if vals else None
    return out


# ------------------------------------------------------------- spec builders -- #
def build_dose_spec(first_line: list[str], second_line: Optional[list[str]] = None,
                    switch_day: Optional[float] = None,
                    dose_scale: Optional[float] = None) -> str:
    """Assemble the ';'-joined --dose string. A second-line dose gets a '*scale'
    suffix (multiply the dose amount) and/or a '@switch_day' suffix (start after the
    first-line readout). Token form is NAME[*SCALE][@START]."""
    parts = list(first_line or [])
    for nm in (second_line or []):
        tok = nm
        if dose_scale is not None and abs(dose_scale - 1.0) > 1e-9:
            tok = f"{tok}*{dose_scale:g}"
        if switch_day is not None:
            tok = f"{tok}@{switch_day:g}"
        parts.append(tok)
    return ";".join(p for p in parts if p)


def build_sample_spec(bounds: dict) -> str:
    """{'F_TNFa': (lo, hi, 'log'), ...} -> 'F_TNFa,lo,hi,log;...' for the sampler.
    A 2-tuple defaults to linear scale."""
    parts = []
    for name, b in (bounds or {}).items():
        lo, hi = b[0], b[1]
        scale = b[2] if len(b) > 2 else "lin"
        parts.append(f"{name},{lo:g},{hi:g},{scale}")
    return ";".join(parts)


def build_override_spec(overrides: dict) -> str:
    """{'KD_TCZ': 2.5e-12} -> 'KD_TCZ=2.5e-12' for the harness param_overrides arg."""
    return ";".join(f"{k}={float(v):.6g}" for k, v in (overrides or {}).items())


# ------------------------------------------------------------------ scorers -- #
def score_flagship(predicted: dict, target: dict) -> dict[str, Any]:
    """Score a predicted response against a target (both as percents, matching
    endpoint keys). Metric = mean absolute error in percentage points across the
    endpoints present in both."""
    keys = list(dict.fromkeys(list(predicted.keys()) + list(target.keys())))
    per = {}
    errs = []
    for k in keys:
        p, t = predicted.get(k), target.get(k)
        if isinstance(p, (int, float)) and isinstance(t, (int, float)):
            e = round(abs(p - t), 1)
            per[k] = {"predicted": p, "target": t, "abs_error_pp": e}
            errs.append(e)
    mae = round(sum(errs) / len(errs), 2) if errs else None
    return {"mae_pp": mae, "n_endpoints": len(errs), "per_endpoint": per}


def score_min_dose(second_line: dict, dose_scale: float, target_pct: float,
                   endpoint: str = None) -> dict[str, Any]:
    """Score the 'minimum effective dose' objective: the regimen must reach
    ``endpoint`` >= ``target_pct`` (%) at the LOWEST dose. A regimen that meets the
    target at a smaller scale scores better; one that misses fails regardless of
    dose. ``endpoint`` defaults to the first response role in ``second_line``."""
    if endpoint is None:
        endpoint = next((k for k in second_line if k != "n_subgroup"), None)
    achieved = second_line.get(endpoint)
    met = isinstance(achieved, (int, float)) and achieved >= target_pct
    return {
        "endpoint": endpoint,
        "target_pct": target_pct,
        "achieved_pct": achieved,
        "dose_scale": dose_scale,
        "target_met": bool(met),
        "response_per_dose": (round(achieved / dose_scale, 1)
                              if met and dose_scale else None),
        "verdict": ("meets target at dose x%g" % dose_scale if met
                    else "MISSES target (%s %s%% < %s%%)"
                    % (endpoint, achieved, target_pct)),
    }


def score_vpop(das_values, target: dict) -> dict[str, Any]:
    """Score a generated population's baseline severity against the clinical target:
    the fraction inside the active band (yield), the accepted cohort's mean/sd, and a
    distribution distance = |mean-target_mean| + |sd-target_sd| (on accepted
    patients). Lower distance + higher yield = a better virtual population."""
    lo, hi = target["band"]
    xs = _finite(das_values)
    n = len(xs)
    accepted = [x for x in xs if lo <= x <= hi]
    yield_pct = round(100.0 * len(accepted) / n, 1) if n else None
    if accepted:
        m = sum(accepted) / len(accepted)
        sd = (sum((x - m) ** 2 for x in accepted) / len(accepted)) ** 0.5
        dist = round(abs(m - target["mean"]) + abs(sd - target["sd"]), 3)
    else:
        m = sd = dist = None
    return {
        "n": n, "n_accepted": len(accepted), "yield_pct": yield_pct,
        "accepted_mean": round(m, 3) if m is not None else None,
        "accepted_sd": round(sd, 3) if sd is not None else None,
        "target_mean": target["mean"], "target_sd": target["sd"],
        "band": target["band"], "distribution_distance": dist,
    }


def select_to_moments(das_values, target: dict) -> dict[str, Any]:
    """Numerically SELECT a virtual population from a sampled pool to match a target
    severity distribution -- the standard QSP prevalence-weighting method (which the
    paper's genetic algorithm implements). Each in-band candidate gets an importance
    weight = target_density / pool_density, so the reweighted population's moments
    match the target. Returns the weighted mean/sd, the distance to target, and the
    effective sample size (ESS); a low ESS means the pool did not cover the target
    range.

    The agent's job is to sample a pool WIDE enough to span the target; this routine
    does the selection."""
    import math
    lo, hi = target["band"]
    m, s = target["mean"], target["sd"]
    pool = [x for x in _finite(das_values) if lo <= x <= hi]
    n_in = len(pool)
    if n_in < 5:
        return {"n_pool": len(_finite(das_values)), "n_inband": n_in,
                "ok": False, "reason": "too few in-band candidates to select from"}

    nbins = max(6, int(round(n_in ** 0.5)))
    binw = (hi - lo) / nbins

    def _bin(x):
        return min(max(int((x - lo) / (hi - lo) * nbins), 0), nbins - 1)

    counts = [0] * nbins
    for x in pool:
        counts[_bin(x)] += 1

    def pool_density(x):
        return max(counts[_bin(x)], 1) / (n_in * binw)

    def target_density(x):
        return math.exp(-0.5 * ((x - m) / s) ** 2) / (s * math.sqrt(2 * math.pi))

    w = [target_density(x) / pool_density(x) for x in pool]
    sw = sum(w)
    if sw <= 0:
        return {"n_pool": len(_finite(das_values)), "n_inband": n_in,
                "ok": False, "reason": "degenerate weights"}
    w = [wi / sw for wi in w]
    wmean = sum(wi * x for wi, x in zip(w, pool))
    wsd = (sum(wi * (x - wmean) ** 2 for wi, x in zip(w, pool))) ** 0.5
    ess = 1.0 / sum(wi * wi for wi in w)
    return {
        "ok": True, "n_pool": len(_finite(das_values)), "n_inband": n_in,
        "weighted_mean": round(wmean, 3), "weighted_sd": round(wsd, 3),
        "target_mean": m, "target_sd": s,
        "distribution_distance": round(abs(wmean - m) + abs(wsd - s), 3),
        "effective_sample_size": round(ess, 1),
        "ess_fraction": round(ess / n_in, 3),
    }


def numeric_fit_1d(evaluate, lo: float, hi: float, log: bool = True,
                   max_evals: int = 20) -> dict[str, Any]:
    """Minimize a scalar objective over ONE parameter with a bounded numerical
    optimizer (scipy Brent). ``evaluate(value) -> error`` runs the forward model; the
    optimizer chooses the values. Returns the fitted value, the error there, the
    number of evaluations, and the full (value, error) trace so the caller can judge
    identifiability."""
    import math
    from scipy.optimize import minimize_scalar
    trace: list[dict] = []

    def obj(u):
        val = 10.0 ** u if log else u
        err = evaluate(val)
        e = err if isinstance(err, (int, float)) else 1e9
        trace.append({"value": val, "error": err})
        return e

    a, b = (math.log10(lo), math.log10(hi)) if log else (lo, hi)
    minimize_scalar(obj, bounds=(a, b), method="bounded",
                    options={"maxiter": max_evals, "xatol": (b - a) / 500.0})
    best = min(trace, key=lambda t: (t["error"] if isinstance(t["error"], (int, float))
                                     else 1e9))
    return {"fitted": best["value"], "error": best["error"], "n_evals": len(trace),
            "trace": sorted(trace, key=lambda t: t["value"])}


def score_fit(predicted: dict, target: dict, fitted: dict,
              references: dict) -> dict[str, Any]:
    """Score a calibration: PRIMARY is how well the fitted parameter(s) reproduce the
    observed response (MAE vs target); SECONDARY is how far each fitted value sits
    from its literature reference (log10-fold)."""
    import math
    fit = score_flagship(predicted, target)
    params = {}
    for name, val in (fitted or {}).items():
        ref = (references or {}).get(name)
        lf = None
        if isinstance(val, (int, float)) and isinstance(ref, (int, float)) \
                and val > 0 and ref > 0:
            lf = round(math.log10(val / ref), 2)
        params[name] = {"fitted": val, "reference": ref, "log10_fold_from_ref": lf}
    return {"acr_mae_pp": fit.get("mae_pp"), "per_endpoint": fit.get("per_endpoint"),
            "parameters": params}


# --------------------------------------------------- clinical reference data -- #
def trial_target(trials: dict, drug: str, week: int,
                 correction: str = "raw") -> Optional[dict]:
    """Return {endpoint: pct} for a drug/week from a clinical-trials table.
    ``correction`` is 'raw' (drug arm) or 'pcorr' (drug minus placebo, floored at 0).
    ``trials`` is the config's clinical_trials dict."""
    arm = (trials.get(drug) or {}).get("weeks", {}).get(week) \
        or (trials.get(drug) or {}).get("weeks", {}).get(str(week))
    if not arm:
        return None
    drug_arm = arm["drug"]
    keys = [k for k in drug_arm if k.upper().startswith("ACR")]
    if correction == "raw":
        return {k: drug_arm[k] for k in keys}
    pl = arm.get("placebo", {})
    return {k: round(max(0.0, drug_arm[k] - pl.get(k, 0.0)), 1) for k in keys}
