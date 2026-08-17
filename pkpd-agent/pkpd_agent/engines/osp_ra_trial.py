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
                    switch_day: Optional[float] = None) -> str:
    """Assemble the ';'-joined --dose string. Second-line doses get a '@switch_day'
    suffix so they start after the first-line readout (the sequential switch)."""
    parts = list(first_line or [])
    for nm in (second_line or []):
        parts.append(f"{nm}@{switch_day:g}" if switch_day is not None else nm)
    return ";".join(p for p in parts if p)


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
