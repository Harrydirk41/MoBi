"""Stage-1 (Layer 5): can an LLM supply physiologically-grounded PARAMETER priors?

The layer below network topology is the ~500 constants that make the model quantitative
(rates, EC50s, max-fold effects, clearances). In practice a modeler sets literature-based
priors / bounds for these before fitting - that is a reasoning task (units + physiology ->
order of magnitude), and the Vantage RA paper documents 130 of them with real values
(ESM2). We give the LLM each parameter's NAME + UNITS + cell context and score its value
guess by order-of-magnitude error against the model's value.

The honesty check is the NAIVE BASELINE: predict, for each parameter, the geometric mean
of all model values sharing its units. That "knows the units and the dataset's scale" but
zero parameter-specific biology. If the LLM only matches this baseline, it is doing units
bookkeeping, not biology; to matter it must BEAT it.
"""

from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass

_DATA = os.path.join(os.path.dirname(__file__), "data", "ra_params_esm2.json")


@dataclass(frozen=True)
class Param:
    name: str
    units: str
    section: str
    value: float

    def dimensionless(self) -> bool:
        return self.units.lower() in ("dimensionless", "", "none", "fraction")


def load_truth(path: str = _DATA) -> list[Param]:
    with open(path, encoding="utf-8") as fh:
        return [Param(d["name"], d["units"], d["section"], float(d["value"]))
                for d in json.load(fh)]


def prompt_view(truth: list[Param]) -> list[dict]:
    """What the agent is shown: name, units, cell context - NOT the value."""
    return [{"name": p.name, "units": p.units, "cell_context": p.section} for p in truth]


def _geomean(vals: list[float]) -> float:
    vals = [abs(v) for v in vals if v]
    if not vals:
        return 1.0
    return math.exp(sum(math.log(v) for v in vals) / len(vals))


def unit_geomean_baseline(truth: list[Param]) -> dict[str, float]:
    """Naive predictor: for each parameter, the geometric mean of all model values that
    share its units. Knows units + empirical scale, no per-parameter biology."""
    by_unit: dict[str, list[float]] = {}
    for p in truth:
        by_unit.setdefault(p.units, []).append(p.value)
    gm = {u: _geomean(vs) for u, vs in by_unit.items()}
    return {p.name: gm[p.units] for p in truth}


def _log10_errs(pred: dict[str, float], truth: list[Param]) -> list[tuple[Param, float]]:
    out = []
    for p in truth:
        v = pred.get(p.name)
        if v is None or v == 0 or p.value == 0:
            continue
        out.append((p, abs(math.log10(abs(v) / abs(p.value)))))
    return out


def _summ(errs: list[float]) -> dict:
    if not errs:
        return {"n": 0}
    s = sorted(errs)
    return {
        "n": len(errs),
        "median_log10_err": round(s[len(s) // 2], 2),
        "within_3x": round(sum(e <= math.log10(3) for e in errs) / len(errs), 3),
        "within_10x": round(sum(e <= 1.0 for e in errs) / len(errs), 3),
        "within_100x": round(sum(e <= 2.0 for e in errs) / len(errs), 3),
    }


def score_params(pred: dict[str, float], truth: list[Param]) -> dict:
    """Order-of-magnitude scoring of predicted parameter values, split by dimensionless
    vs dimensional (dimensionless fold-effects cluster near 1 and are easy), plus the
    naive unit-geomean baseline on the SAME covered set for an honest comparison."""
    pairs = _log10_errs(pred, truth)
    covered = [p for p, _ in pairs]
    errs = [e for _, e in pairs]
    dimless = [e for (p, e) in pairs if p.dimensionless()]
    dimens = [e for (p, e) in pairs if not p.dimensionless()]

    base_pred = unit_geomean_baseline(truth)
    base_errs = [e for _, e in _log10_errs(base_pred, covered)]

    worst = sorted(pairs, key=lambda pe: -pe[1])[:12]
    return {
        "n_predicted": len(pred), "n_scored": len(errs), "n_truth": len(truth),
        "overall": _summ(errs),
        "dimensionless": _summ(dimless),
        "dimensional": _summ(dimens),
        "naive_unit_geomean_baseline": _summ(base_errs),
        "beats_baseline": (bool(errs) and bool(base_errs)
                           and _summ(errs)["median_log10_err"]
                           < _summ(base_errs)["median_log10_err"]),
        "worst_misses": [{"name": p.name, "units": p.units, "true": p.value,
                          "pred": pred.get(p.name), "log10_err": round(e, 2)}
                         for p, e in worst],
    }


def clean_predictions(items) -> dict[str, float]:
    """Agent predictions -> {name: value}. Accepts a list of {name, value} or a dict."""
    out: dict[str, float] = {}
    if isinstance(items, dict):
        items = [{"name": k, "value": v} for k, v in items.items()]
    for it in items or []:
        try:
            nm = str(it["name"]).strip()
            out[nm] = float(it["value"])
        except (KeyError, TypeError, ValueError):
            continue
    return out
