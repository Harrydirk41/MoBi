"""Model-agnostic primitives for the Stage-1 benchmarks: Edge, param scoring, edge scoring.

Pure, vocabulary-free helpers shared by the general pipeline (QSPModel, the loop tools, the
LLM extractor). Nothing here knows about any particular disease or model - node names are
just strings. This is the general core that replaced the RA-specific ra_network / ra_params
modules.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

_SIGN = {"Pro": 1, "Anti": -1, "Hill": 1}


# --------------------------------------------------------------- edges ------ #
@dataclass(frozen=True)
class Edge:
    source: str
    sign: int          # +1 promote, -1 inhibit
    target: str
    process: str = ""  # Sec / Prolif / Influx / ...

    def pair(self) -> tuple[str, str]:
        return (self.source, self.target)

    def signed(self) -> tuple[str, int, str]:
        return (self.source, self.sign, self.target)


def score_network(proposed: list[Edge], truth: list[Edge], sign_aware: bool = True) -> dict:
    """Precision / recall / F1 of a proposed edge set vs the truth, with hit/missed/extra."""
    def key(e: Edge):
        return e.signed() if sign_aware else e.pair()

    P = {key(e) for e in proposed}
    T = {key(e) for e in truth}
    hit, missed, extra = P & T, T - P, P - T
    prec = len(hit) / len(P) if P else 0.0
    rec = len(hit) / len(T) if T else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
    return {
        "sign_aware": sign_aware, "n_proposed": len(P), "n_truth": len(T),
        "hit": len(hit), "missed": len(missed), "extra": len(extra),
        "precision": round(prec, 3), "recall": round(rec, 3), "f1": round(f1, 3),
        "missed_edges": sorted(str(x) for x in missed),
        "extra_edges": sorted(str(x) for x in extra),
    }


def score_signs(pred: dict, truth: list[Edge]) -> dict:
    """Isolated sign accuracy vs the majority-class baseline. pred: (source,target)->sign."""
    total = len(truth)
    if total == 0:
        return {"n": 0}
    correct = sum(1 for e in truth
                  if pred.get((e.source, e.target)) is not None
                  and (1 if pred[(e.source, e.target)] >= 0 else -1) == e.sign)
    n_pos = sum(1 for e in truth if e.sign > 0)
    maj = max(n_pos, total - n_pos) / total
    acc = correct / total
    return {"n": total, "correct": correct, "accuracy": round(acc, 3),
            "majority_baseline": round(maj, 3), "beats_majority": acc > maj,
            "frac_positive": round(n_pos / total, 3)}


# ---------------------------------------------------------- parameters ------ #
# Dimensional units that are MODEL-SCALING (per-molecule / per-mL normalization), not
# physiological - unknowable from physiology, so excluded from the fair param test.
_MODEL_SCALING_UNITS = {"nanogram/(molecule*day)", "molecule/(milliliter*day)",
                        "molecule/(ml*day)", "ng/(molecule*day)"}


@dataclass(frozen=True)
class Param:
    name: str
    units: str
    section: str
    value: float

    def dimensionless(self) -> bool:
        return self.units.lower() in ("dimensionless", "", "none", "fraction")

    def model_scaling(self) -> bool:
        return self.units.lower() in _MODEL_SCALING_UNITS

    def physiological(self) -> bool:
        return (not self.dimensionless()) and (not self.model_scaling())


def _geomean(vals: list[float]) -> float:
    vals = [abs(v) for v in vals if v]
    return math.exp(sum(math.log(v) for v in vals) / len(vals)) if vals else 1.0


def unit_geomean_baseline(truth: list) -> dict[str, float]:
    """Naive predictor: each param -> geomean of all model values sharing its units."""
    by_unit: dict[str, list[float]] = {}
    for p in truth:
        by_unit.setdefault(p.units, []).append(p.value)
    gm = {u: _geomean(vs) for u, vs in by_unit.items()}
    return {p.name: gm[p.units] for p in truth}


def _log10_errs(pred: dict[str, float], truth: list) -> list[tuple]:
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
    return {"n": len(errs), "median_log10_err": round(s[len(s) // 2], 2),
            "within_3x": round(sum(e <= math.log10(3) for e in errs) / len(errs), 3),
            "within_10x": round(sum(e <= 1.0 for e in errs) / len(errs), 3),
            "within_100x": round(sum(e <= 2.0 for e in errs) / len(errs), 3)}


def score_params(pred: dict[str, float], truth: list) -> dict:
    """Order-of-magnitude scoring, split dimensionless / physiological / model-scaling,
    against the naive unit-geomean baseline (fair subset = physiological)."""
    pairs = _log10_errs(pred, truth)
    covered = [p for p, _ in pairs]
    errs = [e for _, e in pairs]
    phys = [e for (p, e) in pairs if p.physiological()]
    base_pred = unit_geomean_baseline(truth)
    base_pairs = _log10_errs(base_pred, covered)
    base_errs = [e for _, e in base_pairs]
    base_phys = [e for (p, e) in base_pairs if p.physiological()]
    ll = _summ(phys).get("median_log10_err")
    bb = _summ(base_phys).get("median_log10_err")
    worst = sorted(pairs, key=lambda pe: -pe[1])[:12]
    return {
        "n_predicted": len(pred), "n_scored": len(errs), "n_truth": len(truth),
        "overall": _summ(errs),
        "dimensionless": _summ([e for (p, e) in pairs if p.dimensionless()]),
        "dimensional": _summ([e for (p, e) in pairs if not p.dimensionless()]),
        "physiological": _summ(phys),
        "physiological_baseline": _summ(base_phys),
        "model_scaling": _summ([e for (p, e) in pairs if p.model_scaling()]),
        "beats_physiological_baseline": (ll is not None and bb is not None and ll < bb),
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
            out[str(it["name"]).strip()] = float(it["value"])
        except (KeyError, TypeError, ValueError):
            continue
    return out
