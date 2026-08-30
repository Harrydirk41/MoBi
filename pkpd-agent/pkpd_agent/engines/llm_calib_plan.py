"""L1, done by the AGENT: given the real-data inventory (which parameters have literature
values, which are missing, and what targets/experiments are available), the LLM must plan the
calibration AND judge identifiability - decide, for each missing parameter, whether the
available data can pin it, and if not, what data would. This is the no-oracle judgment: knowing
what you can and cannot determine, without an answer key. We then grade the agent's judgment
against the known identifiability truth.

Pure: the LLM call is injected; grading compares the agent's determinable/not verdicts to truth.
"""

from __future__ import annotations

from .llm_structure import _parse_json

_SYS = (
    "You are calibrating a QSP model subsystem from LIMITED real data. Some parameters have "
    "literature values; others must be inferred from the available targets/experiments. For "
    "EACH missing parameter decide whether the available data can identify it, and how. Reason "
    "about identifiability: a single steady-state target constrains only ONE overall scale (it "
    "pins a baseline rate); it cannot separate several shape parameters that trade off against "
    "each other - those need a per-perturbation experiment (e.g. a single-cytokine dose-"
    "response). Be honest when a parameter is NOT identifiable from what is available. JSON only.")


def propose_plan(known: list, missing: list, available_data: list, call) -> dict:
    """The LLM plans the calibration and judges identifiability. ``known`` / ``missing`` are
    parameter-name lists; ``available_data`` describes the targets/experiments on hand. Returns
    {plan: [{param, determinable: bool, source, method|reason}]}."""
    user = ("Parameters WITH a literature value (fixed): " + ", ".join(known) +
            "\nParameters MISSING (must be inferred): " + ", ".join(missing) +
            "\nAvailable data to calibrate against:\n  - " + "\n  - ".join(available_data) +
            '\n\nFor each MISSING parameter, decide if the available data can identify it. '
            'Return JSON {"plan": [{"param": name, "determinable": true|false, '
            '"source": "which data / method", "reason": "one phrase"}]}.')
    d = _parse_json(call(_SYS, user))
    out = {}
    for e in (d.get("plan") or []):
        if isinstance(e, dict) and e.get("param"):
            out[e["param"]] = {"determinable": bool(e.get("determinable")),
                               "source": e.get("source"), "reason": e.get("reason")}
    return {"plan": out}


def grade_plan(plan: dict, truth: dict) -> dict:
    """Compare the agent's determinable/not verdicts to the known truth ``{param: bool}``.
    Returns accuracy plus the parameters it got wrong (a false 'determinable' is the dangerous
    error - claiming it can pin a parameter the data cannot identify)."""
    p = plan.get("plan", plan)
    correct, wrong, overclaim = 0, [], []
    for param, is_det in truth.items():
        got = p.get(param, {}).get("determinable")
        if got is None:
            wrong.append((param, "no verdict"))
        elif got == is_det:
            correct += 1
        else:
            wrong.append((param, f"said {got}, truth {is_det}"))
            if got and not is_det:                     # claimed identifiable when it is not
                overclaim.append(param)
    n = len(truth) or 1
    return {"accuracy": round(correct / n, 3), "correct": correct, "n": len(truth),
            "wrong": wrong, "overclaimed": overclaim}
