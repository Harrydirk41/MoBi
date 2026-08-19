"""Summarize and score an RA QSP virtual-trial run (the SimBiology analogue of
osp_ddi's scoring layer).

The Vantage RA model computes the whole clinical trial as timed events, so a Vpop
run returns each patient's model-set response flags (see engines/simbiology.py and
RA_QSP_AGENT_TASK.md). This module turns that per-patient CSV into population
response RATES for the two arms - first-line (day 284) and the held-out
second-line TCZ-in-MTX-inadequate-responders (day 600) - and scores a predicted
flagship response against a target (the paper's held-out validation, or the
model's reference-protocol output).

Kept pure-Python and engine-free so the tools and the example scripts share one
definition of "the arms" and "the score".
"""

from __future__ import annotations

from typing import Any, Optional

# drug -> (long name, mechanism, the dose names the model ships)
DRUG_CATALOG: dict[str, dict[str, Any]] = {
    "MTX": {
        "drug": "methotrexate",
        "modality": "conventional synthetic DMARD (small molecule)",
        "mechanism": "antifolate; broad anti-inflammatory - the first-line anchor",
        "doses": ["MTX_15mg_Q1W_SC_t200", "MTX_9mg_Q1W_SC_t200",
                  "MTX_2p5mgx3_Q1W_Oral_t200"],
    },
    "ADA": {
        "drug": "adalimumab",
        "modality": "biologic (anti-TNF-alpha monoclonal antibody)",
        "mechanism": "neutralizes TNF-alpha",
        "doses": ["ADA40mg_Q2W_SC_t200"],
    },
    "TCZ": {
        "drug": "tocilizumab",
        "modality": "biologic (anti-IL-6-receptor monoclonal antibody)",
        "mechanism": "blocks IL-6 signaling",
        "doses": ["TCZ8mgkg_Q4W_IV_t200", "TCZ4mgkg_Q4W_IV_t200"],
    },
    "SEC": {
        "drug": "secukinumab",
        "modality": "biologic (anti-IL-17A monoclonal antibody)",
        "mechanism": "neutralizes IL-17A",
        "doses": ["SEC_150mg_Q4W_t200", "SEC_75mg_Q4W_t200"],
    },
    "ana": {
        "drug": "anakinra",
        "modality": "biologic (IL-1 receptor antagonist)",
        "mechanism": "blocks IL-1 signaling",
        "doses": ["ana_100mg"],
    },
}

# PD parameters an agent can be asked to CALIBRATE (Stage-4 fit task). Reference
# values transcribed from the paper's ESM2 (Model_parameters, PD section). The
# "true" value is the one that reproduces the observed clinical response; the
# literature reference is a secondary sanity check, not ground truth.
FIT_PARAMS: dict[str, dict[str, Any]] = {
    "KD_TCZ": {
        "unit": "M (molar)",
        "reference": 2.5e-12,
        "meaning": "tocilizumab-IL-6R dissociation constant (binding affinity). "
                   "LOWER KD = tighter binding = more potent IL-6R blockade = "
                   "higher ACR response.",
        "search_range": (1e-13, 1e-9),
        "log_scale": True,
        "source": "ESM2 PD parameters (reference 2.5e-12 M)",
    },
}

# Disease-driver parameters an agent can SAMPLE to build a virtual population
# (Stage-3). Nominal = the Vpop1 mean; (lo, hi) = the Vpop1 observed span. These
# are the pro-inflammatory amplification factors (F_*) and cell-baseline growth
# rates (kg_*) that set disease severity, so their spread drives the baseline
# DAS28-CRP distribution.
VPOP_PARAMS: dict[str, dict[str, Any]] = {
    "F_TNFa":  {"nominal": 10.0,   "span": (0.034, 55.9), "meaning": "TNF-alpha amplification factor"},
    "F_IL6":   {"nominal": 8.39,   "span": (0.035, 44.4), "meaning": "IL-6 amplification factor"},
    "F_IL17":  {"nominal": 13.1,   "span": (0.071, 99.8), "meaning": "IL-17 amplification factor"},
    "F_RANTES":{"nominal": 3.47,   "span": (0.024, 28.3), "meaning": "RANTES/CCL5 amplification factor"},
    "F_GMCSF": {"nominal": 5.54,   "span": (0.027, 37.1), "meaning": "GM-CSF amplification factor"},
    "kg_FLS_Baseline":       {"nominal": 5.18e5, "span": (2.1e4, 4.5e6), "meaning": "fibroblast-like synoviocyte baseline growth"},
    "kg_Macrophage_Baseline":{"nominal": 3.42e5, "span": (8.0e3, 7.5e6), "meaning": "macrophage baseline growth"},
}

# clinical target for the baseline DAS28-CRP distribution (the Vpop1 reference,
# an active-RA population): mean/sd to match, and the active-disease band.
VPOP_TARGET = {"mean": 5.12, "sd": 1.24, "band": (3.2, 8.0)}

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


def summarize_run(res: dict) -> dict[str, Any]:
    """A run_vpop result -> {first_line, second_line, das28} population rates.

    first_line  = the model's ACR20/50/70/remission flags at day 284 (all patients).
    second_line = the MTX_NonResp_TCZ_* flags at day 600, denominator = the
                  MTX-inadequate responders (MTX_NonResp==1) only - the flagship.
    """
    cols = res.get("columns") or {}
    n = len(_finite(cols.get("patient", [])))

    first = {
        "n": n,
        "ACR20": _rate(cols.get("ACR20", [])),
        "ACR50": _rate(cols.get("ACR50", [])),
        "ACR70": _rate(cols.get("ACR70", [])),
        "remission": _rate(cols.get("Rem", [])),
    }

    nonresp = cols.get("MTX_NonResp", [])
    def _isnr(v):
        return isinstance(v, (int, float)) and v == v and v >= 0.5
    n_ir = sum(1 for v in nonresp if _isnr(v))

    def _among_ir(key):
        col = cols.get(key, [])
        vals = [v for v, nr in zip(col, nonresp)
                if isinstance(v, (int, float)) and v == v and _isnr(nr)]
        if not vals:
            return None
        return round(100.0 * sum(1 for v in vals if v >= 0.5) / len(vals), 1)

    second = {
        "n_MTX_IR": n_ir,
        "ACR20": _among_ir("TCZ_ACR20"),
        "ACR50": _among_ir("TCZ_ACR50"),
        "ACR70": _among_ir("TCZ_ACR70"),
        "remission": _among_ir("TCZ_Rem"),
    }

    das = {
        "baseline_mean": _mean(cols.get("DAS28_base", [])),
        "readout_mean": _mean(cols.get("DAS28_read", [])),
    }
    return {"first_line": first, "second_line": second, "das28": das}


def build_dose_spec(first_line: list[str], second_line: Optional[list[str]] = None,
                    switch_day: Optional[float] = None,
                    dose_scale: Optional[float] = None) -> str:
    """Assemble the ';'-joined --dose string. A second-line dose gets a '*scale'
    suffix (multiply the dose amount, a continuous dose level) and/or a
    '@switch_day' suffix (start after the first-line readout). Token form is
    NAME[*SCALE][@START]."""
    parts = list(first_line or [])
    for nm in (second_line or []):
        tok = nm
        if dose_scale is not None and abs(dose_scale - 1.0) > 1e-9:
            tok = f"{tok}*{dose_scale:g}"
        if switch_day is not None:
            tok = f"{tok}@{switch_day:g}"
        parts.append(tok)
    return ";".join(p for p in parts if p)


def score_min_dose(second_line: dict, dose_scale: float, acr20_target: float,
                   endpoint: str = "ACR20") -> dict[str, Any]:
    """Score the 'minimum effective dose' objective: the second-line regimen must
    reach ``endpoint`` >= ``acr20_target`` (%) in the MTX-IR subgroup at the LOWEST
    dose. ``dose_scale`` is the committed dose multiple (1.0 = the labeled dose).
    A regimen that meets the target at a smaller scale scores better; one that
    misses the target fails regardless of dose. 'Just use the max dose' meets the
    target but wastes drug, so it is not optimal."""
    achieved = second_line.get(endpoint)
    met = isinstance(achieved, (int, float)) and achieved >= acr20_target
    return {
        "endpoint": endpoint,
        "target_pct": acr20_target,
        "achieved_pct": achieved,
        "dose_scale": dose_scale,
        "target_met": bool(met),
        # efficiency: response delivered per unit dose (higher = less wasteful);
        # only meaningful when the target is met
        "response_per_dose": (round(achieved / dose_scale, 1)
                              if met and dose_scale else None),
        "verdict": ("meets target at dose x%g" % dose_scale if met
                    else "MISSES target (%s %s%% < %s%%)"
                    % (endpoint, achieved, acr20_target)),
    }


def build_sample_spec(bounds: dict) -> str:
    """{'F_TNFa': (lo, hi, 'log'), ...} -> 'F_TNFa,lo,hi,log;...' for sb_sample_vpop.
    A 2-tuple defaults to linear scale."""
    parts = []
    for name, b in (bounds or {}).items():
        lo, hi = b[0], b[1]
        scale = b[2] if len(b) > 2 else "lin"
        parts.append(f"{name},{lo:g},{hi:g},{scale}")
    return ";".join(parts)


def score_vpop(das_values, target: dict = None) -> dict[str, Any]:
    """Score a generated population's baseline DAS28-CRP against the clinical target:
    the fraction inside the active-RA band (yield), the accepted cohort's mean/sd,
    and a distribution distance = |mean-target_mean| + |sd-target_sd| (on accepted
    patients). Lower distance + higher yield = a better virtual population."""
    tgt = target or VPOP_TARGET
    lo, hi = tgt["band"]
    xs = _finite(das_values)
    n = len(xs)
    accepted = [x for x in xs if lo <= x <= hi]
    yield_pct = round(100.0 * len(accepted) / n, 1) if n else None
    if accepted:
        m = sum(accepted) / len(accepted)
        sd = (sum((x - m) ** 2 for x in accepted) / len(accepted)) ** 0.5
        dist = round(abs(m - tgt["mean"]) + abs(sd - tgt["sd"]), 3)
    else:
        m = sd = dist = None
    return {
        "n": n, "n_accepted": len(accepted), "yield_pct": yield_pct,
        "accepted_mean": round(m, 3) if m is not None else None,
        "accepted_sd": round(sd, 3) if sd is not None else None,
        "target_mean": tgt["mean"], "target_sd": tgt["sd"], "band": tgt["band"],
        "distribution_distance": dist,
    }


def build_override_spec(overrides: dict) -> str:
    """{'KD_TCZ': 2.5e-12} -> 'KD_TCZ=2.5e-12' for the harness param_overrides arg."""
    return ";".join(f"{k}={float(v):.6g}" for k, v in (overrides or {}).items())


def score_fit(predicted: dict, target: dict, fitted: dict,
              references: dict) -> dict[str, Any]:
    """Score a Stage-4 calibration: PRIMARY is how well the fitted parameter(s)
    reproduce the observed response (ACR MAE vs target); SECONDARY is how far each
    fitted value sits from its literature reference (log10-fold)."""
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


def score_flagship(predicted: dict, target: dict) -> dict[str, Any]:
    """Score a predicted second-line response against a target (both as percents,
    keys ACR20/ACR50/ACR70[/remission]). Metric = mean absolute error in
    percentage points across the endpoints present in both."""
    keys = ["ACR20", "ACR50", "ACR70", "remission"]
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
