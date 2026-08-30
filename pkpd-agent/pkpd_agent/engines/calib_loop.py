"""Closed-loop calibration: let the LLM DIAGNOSE non-identifiability and ACQUIRE the
information that fixes it, one round at a time - instead of judging bounds once, up front.

The one-shot experiment showed that bounds cannot separate coupled regulators fit to a
single target: a static prior clips the box but adds no equation. A modeller instead runs a
LOOP - fit, see the scatter, and RESPOND by acquiring an isolating experiment for each
parameter (a single-cytokine perturbation that pins one regulator). This module is the
loop's brain: parse which cytokine isolates which regulator, a 1-D bracketed solve to pin a
parameter from its isolating experiment, and a pluggable LLM action-picker. The heavy part
(running perturbation simulations) lives in the runner; everything here is pure and testable.
"""

from __future__ import annotations

import re

from .llm_structure import _parse_json


def regulator_cytokine(param_name: str) -> str:
    """The cytokine that isolates a regulator: 'FLSProlif_MaxbyIL6' -> 'IL6'. The isolating
    experiment for this parameter is 'elevate only <cytokine>, measure the cell's response'."""
    m = re.search(r"[Mm]axby([A-Za-z0-9]+)", param_name)
    return m.group(1) if m else ""


def solve_1d(evaluator, target: float, lo: float, hi: float,
             tol: float = 1e-3, iters: int = 40) -> float:
    """Bisection solve: find the parameter value in [lo, hi] whose ``evaluator(value)`` (a
    monotonic response) equals ``target``. Isolating experiments make each regulator's response
    monotonic in that one parameter, so a 1-D solve pins it. Returns the best value found;
    brackets are widened once if the target is not enclosed. ``evaluator`` is injected (the
    runner's version runs a perturbation simulation; tests pass a synthetic function)."""
    flo, fhi = evaluator(lo) - target, evaluator(hi) - target
    if flo * fhi > 0:                                  # target not bracketed: widen once
        lo, hi = lo / 10, hi * 10
        flo, fhi = evaluator(lo) - target, evaluator(hi) - target
        if flo * fhi > 0:                              # still not bracketed: return nearer end
            return lo if abs(flo) < abs(fhi) else hi
    a, b = lo, hi
    for _ in range(iters):
        mid = (a * b) ** 0.5 if a > 0 and b > 0 else (a + b) / 2   # geometric (params are log-scale)
        fm = evaluator(mid) - target
        if abs(fm) <= tol * max(abs(target), 1e-9):
            return mid
        if (evaluator(a) - target) * fm <= 0:
            b = mid
        else:
            a = mid
    return (a * b) ** 0.5 if a > 0 and b > 0 else (a + b) / 2


_SYS = ("You are calibrating a QSP model. Several coupled parameters cannot be separated by "
        "fitting them to one aggregate steady-state value - the fit is under-determined. Your "
        "job is to decide, each round, what to do next: acquire an ISOLATING experiment for one "
        "parameter (a single-cytokine perturbation that pins it), or stop. Prefer acquiring an "
        "experiment while any coupled parameter is still unpinned. Output JSON only.")


def choose_experiment(state: dict, call) -> dict:
    """Ask the LLM for the next action given the loop state {pinned, remaining, errors}. Returns
    {"action": "experiment", "cytokine": X} or {"action": "stop"}. Pluggable ``call`` for tests."""
    remaining = state.get("remaining", [])
    user = ("State:\n"
            f"  parameters still un-pinned (coupled, non-identifiable): {remaining}\n"
            f"  already pinned this way: {state.get('pinned', [])}\n"
            f"  current mean recovery error: {state.get('mean_error', '?')}\n\n"
            "Each un-pinned parameter has an isolating single-cytokine perturbation experiment "
            "available (elevate that one cytokine, measure the cell's response). Choose the next "
            'action. Return JSON {"action": "experiment"|"stop", "cytokine": <name or null>, '
            '"reason": "one phrase"}. Pick a cytokine from the un-pinned list.')
    d = _parse_json(call(_SYS, user))
    if not isinstance(d, dict) or d.get("action") == "stop":
        return {"action": "stop"}
    cyt = d.get("cytokine")
    if cyt in remaining:
        return {"action": "experiment", "cytokine": cyt, "reason": d.get("reason")}
    return {"action": "stop"}                          # invalid pick -> stop rather than loop
