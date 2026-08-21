"""LLM tasks.json drafter: propose the downstream-task role assignments from a model.

The 2-7 tasks need a projects/<name>/tasks.json. Most of its fields are DERIVABLE from
the model itself - which parameters are disease drivers, which are calibratable PD
constants, which are druggable pathway knobs, and which states are the clinical
readout. This extractor reads network.json and drafts those, so the author reviews a
draft instead of writing from scratch. It is the tasks-layer analogue of
llm_structure.py (structure), and shares its safety discipline:

  1. The drafter LLM is SEPARATE from any benchmarked LLM.
  2. ``compare_tasks`` regresses the draft against a KNOWN-GOOD tasks.json (the
     hand-written RA one) - the role assignments must reproduce before the draft is
     trusted on a new model. Trust the LLM to READ naming; trust the regression to
     prove it read right.

What the drafter CANNOT invent (left as explicit TODO stubs for the author, because it
is external data or lives in the .sbproj, not network.json):
  * clinical_trials / refractory_target / vpop_target - real trial numbers (author's
    validation data),
  * drugs (dose NAMES) - SimBiology dose objects, not in network.json,
  * timeline - the model's event days (not dumped by sb_network_json.m).

The drafter is pluggable via ``call_fn(system, user) -> str`` (tests use a stub).
"""

from __future__ import annotations

import json

from .llm_structure import _parse_json, default_call  # reuse tolerant parser + Claude call

_SYS = ("You are a QSP-model analyst. You read a model's parameters, species and rule "
        "expressions and assign each to its ROLE in a virtual-trial workflow, using the "
        "EXACT names given. You report only what the model encodes; you do not invent "
        "clinical numbers. Output JSON only, no prose.")


def classify_parameters(params: list[dict], rules: list[str], call) -> dict:
    """Assign parameters to task roles from name/units/context. Returns
    {disease_drivers, calibratable, druggable} lists of parameter names (a parameter
    may appear in more than one - a pathway amplification factor is both a Vpop driver
    and a design target)."""
    compact = [{"name": p.get("name"), "units": p.get("units"), "value": p.get("value")}
               for p in params]
    user = (
        "Assign each parameter to any of these roles (a parameter may take more than "
        "one, or none):\n"
        "- 'disease_drivers': pro-inflammatory amplification factors / baseline growth "
        "rates that set disease SEVERITY - the knobs you would vary to build a virtual "
        "population (e.g. pathway amplification factors, cell baseline growth rates).\n"
        "- 'druggable': the subset of disease drivers a drug could plausibly SUPPRESS "
        "(a pathway knob a therapeutic targets).\n"
        "- 'calibratable': PD / binding constants you would FIT to reproduce an observed "
        "drug response (e.g. a dissociation constant KD, an EC50, a kon/koff).\n"
        'Return JSON {"disease_drivers": [...], "druggable": [...], "calibratable": '
        "[...]} using the exact parameter names.\n\n"
        f"PARAMETERS:\n{json.dumps(compact)}\n\nRULES (context):\n" + "\n".join(rules[:60]))
    r = _parse_json(call(_SYS, user))
    return {k: [n for n in r.get(k, []) if isinstance(n, str)]
            for k in ("disease_drivers", "druggable", "calibratable")}


def classify_readout_states(species: list[str], rules: list[str], call) -> dict:
    """Identify the clinical-readout states and their trial roles. Returns
    {first_line_flags, subgroup_flag, second_line_flags, severity_states} using exact
    species names - the pieces that assemble into readout_states + run_columns."""
    user = (
        "These are a QSP model's states and the rules that set them. Identify the "
        "CLINICAL-READOUT states and their role in a two-line trial:\n"
        "- 'first_line_flags': response flags set at the first-line readout (e.g. "
        "graded response levels + a remission flag).\n"
        "- 'subgroup_flag': the single flag marking inadequate responders escalated to "
        "second line (empty string if the model has none).\n"
        "- 'second_line_flags': the response flags for that subgroup at the second "
        "readout (parallel to first_line_flags).\n"
        "- 'severity_states': the continuous disease-severity state and its captured "
        "baseline (e.g. a composite activity score and its baseline copy).\n"
        'Return JSON {"first_line_flags": [...], "subgroup_flag": "...", '
        '"second_line_flags": [...], "severity_states": [...]} using exact names.\n\n'
        f"STATES:\n{json.dumps(species)}\n\nRULES:\n" + "\n".join(rules[:80]))
    r = _parse_json(call(_SYS, user))
    return {
        "first_line_flags": [s for s in r.get("first_line_flags", []) if isinstance(s, str)],
        "subgroup_flag": r.get("subgroup_flag") or "",
        "second_line_flags": [s for s in r.get("second_line_flags", []) if isinstance(s, str)],
        "severity_states": [s for s in r.get("severity_states", []) if isinstance(s, str)],
    }


def _readout_states(ro: dict) -> list[str]:
    """Assemble the role-ordered readout_states list from the classified pieces."""
    out = list(ro["first_line_flags"])
    if ro["subgroup_flag"]:
        out.append(ro["subgroup_flag"])
    out += list(ro["second_line_flags"])
    out += list(ro["severity_states"])
    return out


_TODO = "TODO: author must supply (external clinical data or .sbproj dose names)"


def draft_tasks(data: dict, call, name: str = "QSP model") -> dict:
    """Draft a tasks.json dict from a network.json dict. Model-derivable fields are
    filled; genuinely-external fields carry a TODO stub for the author. Never trusted
    until compare_tasks regresses it against a known-good reference."""
    species = [s["name"] for s in data.get("species", [])]
    params = data.get("parameters", [])
    rules = [r.get("rule", "") if isinstance(r, dict) else str(r)
             for r in data.get("rules", [])]

    roles = classify_parameters(params, rules, call)
    ro = classify_readout_states(species, rules, call)
    unit_of = {p.get("name"): p.get("units", "") for p in params}

    def _pmap(names):
        return {n: {"meaning": f"TODO: describe {n}", "units": unit_of.get(n, "")}
                for n in names}

    return {
        "name": name,
        "_derived_from": "network.json (llm_tasks.draft_tasks) - REVIEW before use",
        "readout_states": _readout_states(ro),
        "readout_roles": ro,
        "vpop_drivers": {n: {"nominal": None, "span": None,
                             "meaning": f"TODO: describe {n}"}
                         for n in roles["disease_drivers"]},
        "design_targets": {n: {"pathway": f"TODO: pathway for {n}", "analogue": "TODO",
                               "note": ""} for n in roles["druggable"]},
        "fit_params": {n: {"unit": unit_of.get(n, ""), "reference": None,
                           "meaning": f"TODO: describe {n}", "search_range": None,
                           "log_scale": True} for n in roles["calibratable"]},
        # external / .sbproj-only -> explicit stubs, not guesses
        "drugs": _TODO,
        "timeline": _TODO,
        "vpop_target": _TODO,
        "clinical_trials": _TODO,
        "refractory_target": _TODO,
    }


def compare_tasks(draft: dict, reference: dict) -> dict:
    """Regression: do the draft's role assignments reproduce a known-good tasks.json?
    Reports precision/recall/F1 for each derivable field (as name sets)."""
    def prf(pred: set, truth: set) -> dict:
        hit = len(pred & truth)
        p = hit / len(pred) if pred else 0.0
        r = hit / len(truth) if truth else 0.0
        f = 2 * p * r / (p + r) if (p + r) else 0.0
        return {"precision": round(p, 3), "recall": round(r, 3), "f1": round(f, 3),
                "hit": hit, "n_pred": len(pred), "n_truth": len(truth),
                "missed": sorted(truth - pred)[:15], "extra": sorted(pred - truth)[:15]}

    def keys(d, field):
        v = d.get(field)
        return set(v.keys()) if isinstance(v, dict) else set()

    return {
        "readout_states": prf(set(draft.get("readout_states", [])),
                              set(reference.get("readout_states", []))),
        "vpop_drivers": prf(keys(draft, "vpop_drivers"), keys(reference, "vpop_drivers")),
        "design_targets": prf(keys(draft, "design_targets"),
                              keys(reference, "design_targets")),
        "fit_params": prf(keys(draft, "fit_params"), keys(reference, "fit_params")),
    }
