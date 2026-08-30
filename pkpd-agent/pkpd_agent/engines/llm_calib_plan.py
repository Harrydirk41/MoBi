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


_SYS_OPEN = (
    "You are planning to build and calibrate a QSP model subsystem from the data you actually "
    "have. You are given the full parameter list (some parameters have a literature value, some "
    "do not) and a list of the datasets available. Plan the calibration: for each parameter "
    "WITHOUT a literature value, decide whether the AVAILABLE data can identify it; if it "
    "cannot, state what experiment WOULD identify it and whether that experiment is in your "
    "available data. Be rigorous about identifiability and about what you are missing - do NOT "
    "assume a dataset exists unless it is listed. JSON only.")


def propose_plan_open(params: list, available_data: list, call) -> dict:
    """Harder, un-scaffolded plan: ``params`` is [{name, has_literature_value}]; the agent must
    itself work out which unknowns the data can pin, what experiment each needs, and whether that
    experiment is present in ``available_data`` (which may contain distractors and does NOT
    announce what is missing). Returns {plan: {param: {determinable, needs, needs_available}}}."""
    plist = "\n".join(f"  {p['name']} "
                      f"({'has literature value' if p.get('has_literature_value') else 'NO value'})"
                      for p in params)
    user = ("Full parameter list:\n" + plist +
            "\n\nDatasets available to you:\n  - " + "\n  - ".join(available_data) +
            '\n\nFor each parameter with NO literature value, return JSON {"plan": [{"param": '
            'name, "determinable": true|false, "needs": "the experiment that would identify it", '
            '"needs_available": true|false (is that experiment in your available datasets?)}]}.')
    d = _parse_json(call(_SYS_OPEN, user))
    out = {}
    for e in (d.get("plan") or []):
        if isinstance(e, dict) and e.get("param"):
            out[e["param"]] = {"determinable": bool(e.get("determinable")),
                               "needs": e.get("needs"),
                               "needs_available": e.get("needs_available")}
    return {"plan": out}


def grade_open(plan: dict, truth: dict, needs_kw: list) -> dict:
    """Grade the un-scaffolded plan on three axes: identifiability accuracy, and - for the
    NOT-identifiable params - whether the agent (a) named the right missing experiment (a keyword
    in ``needs_kw``, e.g. 'dose', 'perturb', 'titration') and (b) correctly flagged it as NOT in
    the available data. The un-scaffolded win is realising, unprompted, that the identifying data
    is absent."""
    p = plan.get("plan", plan)
    id_correct, named_need, flagged_absent, overclaim = 0, 0, 0, []
    non_ident = [k for k, v in truth.items() if not v]
    for param, is_det in truth.items():
        v = p.get(param, {})
        got = v.get("determinable")
        if got == is_det:
            id_correct += 1
        if got and not is_det:
            overclaim.append(param)
        if not is_det:                                 # should need an absent perturbation exp
            needs = (v.get("needs") or "").lower()
            if any(kw in needs for kw in needs_kw):
                named_need += 1
            if v.get("needs_available") is False:
                flagged_absent += 1
    n, m = len(truth) or 1, len(non_ident) or 1
    return {"id_accuracy": round(id_correct / n, 3),
            "named_missing_experiment": f"{named_need}/{len(non_ident)}",
            "flagged_data_absent": f"{flagged_absent}/{len(non_ident)}",
            "overclaimed": overclaim}


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
