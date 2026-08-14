"""Numerical parameter identification for the OSP PBPK engine.

The LLM decides the MODEL - structure (methods/processes) and which parameters
to estimate vs fix. This module does the mechanical part the LLM shouldn't: it
numerically fits the chosen parameters to the observed data with a derivative-
free optimizer (Nelder-Mead in log space), exactly like PK-Sim's Parameter
Identification module, but headless.

Each objective evaluation runs PK-Sim, so it is expensive; the optimizer fits
against a representative SUBSET of simulations (chosen to span routes and doses)
and evaluates the final fit on the full set. It reports parameters that hit
their bounds - the signal of an unidentifiable parameter or a wrong structure,
which is what the LLM reasons about between rounds.
"""

from __future__ import annotations

import math
import re
from typing import Any

from . import osp_score
from .osp_cli import OSPCli


# --------------------------------------------------------------------------- #
# choose a representative subset of simulations to fit against (for speed)
# --------------------------------------------------------------------------- #

def pick_subset(sim_names: list[str], k: int = 4) -> list[str]:
    """One simulation per (route, dose) bucket, capped at k, spanning routes."""
    buckets: dict[tuple, str] = {}
    for n in sim_names:
        _, route, dose = OSPCli._parse_sim_name(n)
        key = (route, dose)
        buckets.setdefault(key, n)
    picked = list(buckets.values())
    # ensure both routes represented if present
    iv = [n for n in picked if OSPCli._parse_sim_name(n)[1] == "IV"]
    po = [n for n in picked if OSPCli._parse_sim_name(n)[1] == "PO"]
    out = []
    for pool in (iv, po):
        out.extend(pool[:max(1, k // 2)])
    for n in picked:
        if n not in out and len(out) < k:
            out.append(n)
    return out[:k] or sim_names[:k]


# --------------------------------------------------------------------------- #
# objective: mean squared log fold error over matched points
# --------------------------------------------------------------------------- #

def _log_sse(observed: list[dict], predicted: list[dict]) -> tuple[float, int]:
    preds = {p["dataset"]: p for p in predicted}
    total, n = 0.0, 0
    for o in observed:
        pr = preds.get(o["dataset"])
        if not pr:
            continue
        for t, c in zip(o["time_h"], o["conc_mg_L"]):
            of, tf = osp_score._finite(c), osp_score._finite(t)
            if of is None or of <= 0 or tf is None:
                continue
            p = osp_score._interp(pr["time_h"], pr["pred_conc_mg_L"], tf)
            if p is None or p <= 0:
                continue
            total += (math.log(p) - math.log(of)) ** 2
            n += 1
    return (total / n if n else 1e9), n


# --------------------------------------------------------------------------- #
# optimize
# --------------------------------------------------------------------------- #

def run_optimization(cli: OSPCli, snapshot_path: str, observed: list[dict],
                     estimate: dict[str, list], fix: dict | None = None,
                     structure: dict | None = None,
                     fit_simulations: list[str] | None = None,
                     max_evals: int = 30, on_eval=None) -> dict[str, Any]:
    """Fit ``estimate`` (name -> [lo, hi]) to the observed data.

    ``on_eval(i, values, log_sse)`` (optional) is called after each objective
    evaluation for live progress.

    Returns optimized values, the full-set fit, evals used, and any parameters
    that hit their bounds (identifiability / structure warning).
    """
    import numpy as np
    from scipy.optimize import minimize

    names = list(estimate.keys())
    los = np.array([float(estimate[n][0]) for n in names], float)
    his = np.array([float(estimate[n][1]) for n in names], float)
    if np.any(los <= 0):
        return {"ok": False, "message": "bounds must be positive (log-space fit)"}

    # subset of simulations to fit against
    all_sims = cli.simulation_names(snapshot_path)
    subset = fit_simulations or pick_subset(all_sims, k=4)
    subset_studies = {osp_score._norm_study(OSPCli._parse_sim_name(s)[0]) for s in subset}
    observed_sub = [o for o in observed
                    if osp_score._norm_study(osp_score._obs_key(o)[0]) in subset_studies] \
        or observed

    base_edits: dict[str, Any] = dict(structure or {})
    base_edits.setdefault("parameters", {})
    base_edits["parameters"].update(fix or {})

    lx, hx = np.log10(los), np.log10(his)
    history: list[dict] = []

    def eval_at(values: dict[str, float], sims):
        edits = {**base_edits,
                 "parameters": {**base_edits["parameters"], **values}}
        # during the fit, prune the snapshot to the subset so snap builds ONLY
        # those simulations (snap otherwise rebuilds all of them every eval).
        res = cli.build_and_run(snapshot_path, edits=edits, simulations=sims,
                                prune_simulations=sims is subset)
        if not res["ok"]:
            return None, res
        predicted, _ = osp_score.map_predictions(
            res["profiles"],
            observed_sub if sims is subset else observed)
        return predicted, res

    def objective(x):
        xc = np.clip(x, lx, hx)
        penalty = float(np.sum((x - xc) ** 2)) * 10.0
        values = {n: float(10 ** xc[i]) for i, n in enumerate(names)}
        predicted, res = eval_at(values, subset)
        if predicted is None:
            history.append({"values": values, "obj": None, "error": res.get("message")})
            if on_eval:
                on_eval(len(history), values, None)
            return 1e6 + penalty
        sse, npts = _log_sse(observed_sub, predicted)
        history.append({"values": values, "log_sse": round(sse, 5), "n": npts})
        if on_eval:
            on_eval(len(history), values, round(sse, 5))
        return sse + penalty

    x0 = (lx + hx) / 2.0
    res = minimize(objective, x0, method="Nelder-Mead",
                   options={"maxfev": max_evals, "xatol": 0.02, "fatol": 1e-3})
    best = np.clip(res.x, lx, hx)
    optimized = {n: float(10 ** best[i]) for i, n in enumerate(names)}

    # parameters that hit a bound (within ~3% in log space) -> identifiability flag
    at_bound = []
    for i, n in enumerate(names):
        span = hx[i] - lx[i]
        if span > 0 and (best[i] - lx[i] < 0.03 * span):
            at_bound.append({"parameter": n, "value": optimized[n], "bound": "lower"})
        elif span > 0 and (hx[i] - best[i] < 0.03 * span):
            at_bound.append({"parameter": n, "value": optimized[n], "bound": "upper"})

    # LOCAL SENSITIVITY (identifiability EVIDENCE, not a hard-coded opinion):
    # perturb each fitted parameter around the optimum and measure how much the
    # objective actually moves. A parameter the data barely responds to is weakly
    # identified - but we report the NUMBER and let the agent/report reason about
    # it, rather than asserting which parameters are unidentifiable.
    sensitivity = _local_sensitivity(eval_at, _log_sse, observed_sub, subset,
                                     optimized, names, best, lx, hx)

    # final fit on the FULL observed set
    predicted_full, res_full = eval_at(optimized, None)
    if predicted_full is None:
        return {"ok": False, "message": f"final run failed: {res_full.get('message')}",
                "optimized": optimized}
    score = osp_score.score_fit(observed, predicted_full)
    return {
        "ok": True,
        "optimized": optimized,
        "fit": score["overall"],
        "by_route": score["by_route"],
        "worst_datasets": [{"dataset": d["dataset"], "route": d["route"],
                            "gmfe": d["gmfe"], "bias": d["bias"]}
                           for d in score["per_dataset"][:3]],
        "params_at_bound": at_bound,
        "sensitivity": sensitivity,
        "n_evals": len(history),
        "fit_simulations": subset,
    }


def _local_sensitivity(eval_at, log_sse, observed_sub, subset, optimized, names,
                       best, lx, hx, step: float = 0.15) -> dict:
    """One-at-a-time local sensitivity around the fitted optimum.

    For each estimated parameter, multiply/divide it by 10**step (others held at
    the optimum) and measure |Δ log-objective|. Normalised to the most sensitive
    parameter so the result is a data-driven relative ranking; a ``relative``
    near 0 means the data hardly constrain that parameter (weakly identifiable).
    Returns {param: {obj_change, relative}} - the CONSUMER decides what it means."""
    import numpy as np
    base_pred, _ = eval_at(optimized, subset)
    if base_pred is None:
        return {}
    base = log_sse(observed_sub, base_pred)[0]
    raw: dict[str, float] = {}
    for i, n in enumerate(names):
        deltas = []
        for sign in (+1.0, -1.0):
            xp = float(np.clip(best[i] + sign * step, lx[i], hx[i]))
            if abs(xp - best[i]) < 1e-9:
                continue
            pert = dict(optimized)
            pert[n] = float(10 ** xp)
            pred, _ = eval_at(pert, subset)
            if pred is None:
                continue
            deltas.append(abs(log_sse(observed_sub, pred)[0] - base))
        raw[n] = (sum(deltas) / len(deltas)) if deltas else 0.0
    top = max(raw.values()) if raw else 0.0
    return {n: {"obj_change": round(v, 5),
                "relative": round(v / top, 3) if top > 0 else 0.0}
            for n, v in raw.items()}
