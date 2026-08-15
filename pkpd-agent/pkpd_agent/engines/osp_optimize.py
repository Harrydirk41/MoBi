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

# abort an optimization this many failed evals in if the model has never once
# built/run - the failure is structural, so more evals cannot help.
_ABORT_AFTER = 3


class _ModelNeverRan(Exception):
    """Raised inside the objective when the model fails to build/run repeatedly
    before any successful evaluation - a structural error, not a fitting one."""

    def __init__(self, detail):
        self.detail = detail or "no error detail from PK-Sim"
        super().__init__(self.detail)


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
    # the fit runs in log10 space, so bounds must be strictly positive with hi>lo.
    # A non-positive UPPER bound (or hi<=lo) is genuinely broken -> fail. But a
    # lower bound of 0 is a common, benign ">= 0" intent; clamp it to a small
    # fraction of the upper bound rather than failing the whole call and forcing
    # the agent to burn a retry.
    if np.any(his <= 0) or np.any(his <= los):
        return {"ok": False, "message": "each parameter needs bounds [lo, hi] with "
                "0 <= lo < hi and hi > 0 (the fit is in log space)"}
    clamped = [names[i] for i in range(len(los)) if los[i] <= 0]
    los = np.where(los <= 0, his * 1e-6, los)

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

    state = {"n_ok": 0}                       # did the model ever run at all?

    def objective(x):
        xc = np.clip(x, lx, hx)
        penalty = float(np.sum((x - xc) ** 2)) * 10.0
        values = {n: float(10 ** xc[i]) for i, n in enumerate(names)}
        predicted, res = eval_at(values, subset)
        if predicted is None:
            err = res.get("message")
            history.append({"values": values, "obj": None, "error": err})
            if on_eval:
                on_eval(len(history), values, None, err)
            # if the model has NEVER run and keeps failing to build/run, the
            # problem is the model configuration (structure/methods/process), not
            # the parameter values - every future eval will fail identically.
            # Abort early with the real error instead of burning the whole budget.
            if state["n_ok"] == 0 and len(history) >= _ABORT_AFTER:
                raise _ModelNeverRan(err)
            return 1e6 + penalty
        state["n_ok"] += 1
        sse, npts = _log_sse(observed_sub, predicted)
        history.append({"values": values, "log_sse": round(sse, 5), "n": npts})
        if on_eval:
            on_eval(len(history), values, round(sse, 5))
        return sse + penalty

    x0 = (lx + hx) / 2.0
    try:
        res = minimize(objective, x0, method="Nelder-Mead",
                       options={"maxfev": max_evals, "xatol": 0.02, "fatol": 1e-3})
    except _ModelNeverRan as exc:
        return {"ok": False, "run_never_succeeded": True,
                "message": (f"the model failed to build/run on all {len(history)} "
                            f"attempts before any parameter could be fitted: "
                            f"{exc.detail}. This is a STRUCTURE/configuration "
                            "problem (calculation methods, processes, or a "
                            "parameter key that does not exist), NOT a parameter-"
                            "value or bounds problem - fix the model setup, not the "
                            "bounds."),
                "n_evals": len(history)}
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
        "bounds_clamped": clamped or None,
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
    coll = _collinearity(eval_at, log_sse, observed_sub, subset, optimized,
                         names, best, lx, hx, base, step)
    out = {}
    for n, v in raw.items():
        partner, c = coll.get(n, (None, 0.0))
        out[n] = {"obj_change": round(v, 5),
                  "relative": round(v / top, 3) if top > 0 else 0.0,
                  "collinearity": c, "collinear_with": partner}
    return out


def _collinearity(eval_at, log_sse, observed_sub, subset, optimized, names,
                  best, lx, hx, base, step: float = 0.15) -> dict:
    """Pairwise collinearity (trade-off) around the fitted optimum.

    One-at-a-time sensitivity is blind to *correlated* parameters: two clearances
    acting on the same parent can each look influential while only their sum is
    identifiable (their split trades off freely). At the optimum the gradient is
    ~0, so a perturbation's effect is curvature-driven. For a pair (i, j) we
    compare the two joint diagonal moves - both-up (i+, j+) and opposed (i+, j-).
    If one direction is much flatter than the other, there is a near-flat
    trade-off direction: the pair is collinear and their individual values are
    not separately identifiable (only their combination is constrained).

    collinearity = 1 - min(d_up, d_opp) / max(d_up, d_opp), in [0, 1]; ~1 means a
    flat trade-off exists. Returns {param: (strongest_partner, collinearity)}.
    Capped at 8 fitted parameters to bound the O(n^2) extra evaluations."""
    import numpy as np
    n = len(names)
    best_pair = {nm: (None, 0.0) for nm in names}
    if n < 2 or n > 8:
        return best_pair

    def joint_obj(i, dj_i, j, dj_j):
        pert = dict(optimized)
        for k, dk in ((i, dj_i), (j, dj_j)):
            xp = float(np.clip(best[k] + dk, lx[k], hx[k]))
            pert[names[k]] = float(10 ** xp)
        pred, _ = eval_at(pert, subset)
        if pred is None:
            return None
        return abs(log_sse(observed_sub, pred)[0] - base)

    for i in range(n):
        for j in range(i + 1, n):
            d_up = joint_obj(i, +step, j, +step)
            d_opp = joint_obj(i, +step, j, -step)
            if d_up is None or d_opp is None:
                continue
            hi, lo = max(d_up, d_opp), min(d_up, d_opp)
            if hi <= 1e-9:
                continue
            c = round(1.0 - lo / hi, 3)
            if c > best_pair[names[i]][1]:
                best_pair[names[i]] = (names[j], c)
            if c > best_pair[names[j]][1]:
                best_pair[names[j]] = (names[i], c)
    return best_pair
