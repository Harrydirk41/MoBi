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

from .llm_structure import _parse_json  # reuse the tolerant JSON extractor

_SYS = ("You are a QSP-model analyst. You read a model's parameters, species and rule "
        "expressions and assign each to its ROLE in a virtual-trial workflow, using the "
        "EXACT names given. You report only what the model encodes; you do not invent "
        "clinical numbers. Output JSON only, no prose.")


def default_call(config):
    """A real Claude call: (system, user) -> text. Larger max_tokens than the structure
    extractor (role lists over hundreds of parameters), and a clear error on an empty
    reply so a truncated/blocked response is diagnosable instead of a cryptic JSON error."""
    import anthropic
    client = anthropic.Anthropic()

    def call(system: str, user: str) -> str:
        resp = client.messages.create(
            model=config.model, max_tokens=16000, system=system,
            messages=[{"role": "user", "content": user}])
        text = "".join(getattr(b, "text", "") for b in resp.content
                       if getattr(b, "type", None) == "text")
        if not text.strip():
            kinds = [getattr(b, "type", "?") for b in resp.content]
            raise RuntimeError(
                f"LLM returned no text (stop_reason={resp.stop_reason}, "
                f"content blocks={kinds}). If stop_reason is 'max_tokens', the reply "
                "was truncated - the drafter batches parameters to avoid this, so re-run.")
        return text
    return call


def default_web_call(config, max_searches: int = 8, tool_version: str = "web_search_20260209"):
    """Like ``default_call`` but with the server-side web-search tool enabled, so the model
    reads the literature ITSELF during the call (Route B) instead of being handed pre-fetched
    text. The search runs on Anthropic's infrastructure and loops server-side; we just read the
    final text blocks. Non-deterministic by nature (live search) - that is the trade vs the
    cached Route-A material. Falls back to the older tool id if the newer one is rejected."""
    import anthropic
    client = anthropic.Anthropic()

    def call(system: str, user: str) -> str:
        for ver in (tool_version, "web_search_20250305"):
            try:
                resp = client.messages.create(
                    model=config.model, max_tokens=16000, system=system,
                    tools=[{"type": ver, "name": "web_search", "max_uses": max_searches}],
                    messages=[{"role": "user", "content": user}])
                break
            except anthropic.BadRequestError:
                if ver == "web_search_20250305":
                    raise
        text = "".join(getattr(b, "text", "") for b in resp.content
                       if getattr(b, "type", None) == "text")
        if not text.strip():
            kinds = [getattr(b, "type", "?") for b in resp.content]
            raise RuntimeError(
                f"web-enabled LLM returned no text (stop_reason={resp.stop_reason}, "
                f"content blocks={kinds}).")
        return text
    return call


_BOUNDS_SYS = (
    "You are a systems-biology modeller setting PHYSIOLOGICALLY PLAUSIBLE bounds for model "
    "parameters before calibration. For each parameter, give a lower and upper bound that a "
    "reasonable modeller would allow it to vary within - tight enough to reflect biology "
    "(e.g. a cytokine's fold-change effect on a cell rate is modest, not 100x), wide enough "
    "not to exclude the true value. Reason ONLY from the parameter's name/meaning and general "
    "biology, NOT from any specific number you might recall. Output JSON only.")


def propose_bounds(params: list, call) -> dict:
    """Ask the LLM for plausible [lo, hi] bounds per parameter - the 'biological judgment' that
    regularizes an ill-posed calibration. ``params`` is [{name, units, meaning?}]; returns
    {name: (lo, hi)} for the entries the LLM bounded sensibly (lo>0, hi>lo). Pluggable ``call``
    so tests use a stub."""
    lines = "\n".join(
        f"- {p['name']} ({p.get('units', '?')}): {p.get('meaning', '')}".rstrip()
        for p in params)
    user = ("Give plausible calibration bounds for these parameters:\n" + lines +
            '\n\nReturn JSON {"bounds": [{"name": ..., "lo": number, "hi": number, '
            '"basis": "one phrase"}]}. Bounds must be positive with hi > lo.')
    d = _parse_json(call(_BOUNDS_SYS, user))
    out: dict = {}
    for b in (d.get("bounds") or []):
        if not isinstance(b, dict):
            continue
        name, lo, hi = b.get("name"), b.get("lo"), b.get("hi")
        try:
            lo, hi = float(lo), float(hi)
        except (TypeError, ValueError):
            continue
        if name and 0 < lo < hi:
            out[name] = (lo, hi)
    return out


def _merge_roles(dst: dict, src: dict) -> None:
    for k in ("disease_drivers", "druggable", "calibratable"):
        for n in src.get(k, []):
            if isinstance(n, str) and n not in dst[k]:
                dst[k].append(n)


def classify_parameters(params: list[dict], rules: list[str], call,
                        batch: int = 120) -> dict:
    """Assign parameters to task roles from name/units/context. Returns
    {disease_drivers, calibratable, druggable} lists of parameter names (a parameter
    may appear in more than one - a pathway amplification factor is both a Vpop driver
    and a design target). Batched over the parameter list so a many-hundred-parameter
    model does not truncate the reply."""
    out = {"disease_drivers": [], "druggable": [], "calibratable": []}
    rule_ctx = "\n".join(rules[:60])
    for k in range(0, len(params), batch):
        chunk = params[k:k + batch]
        compact = [{"name": p.get("name"), "units": p.get("units"),
                    "value": p.get("value")} for p in chunk]
        user = (
            "Assign each parameter to any of these roles. Be STRICT: the vast majority of "
            "parameters take NONE. Return a parameter ONLY if it is a top-level knob for "
            "that role, never a per-reaction mechanistic coefficient.\n"
            "- 'disease_drivers': GLOBAL, model-wide amplification factors or baseline "
            "growth/production rates that set overall disease SEVERITY - the handful of "
            "knobs you would vary to build a virtual population. EXCLUDE per-interaction "
            "coefficients (a parameter naming two specific species, a max-effect strength, "
            "a Hill/EC50/half-effect constant, a rule-intermediate 'effect' term).\n"
            "- 'druggable': the subset of those GLOBAL disease drivers a therapeutic could "
            "suppress. Same exclusions.\n"
            "- 'calibratable': target-binding or potency constants of a DRUG you would FIT "
            "to reproduce an observed response (a dissociation constant KD, an EC50, a "
            "kon/koff). EXCLUDE PK disposition parameters (clearance CL, volume, "
            "bioavailability F, absorption) and endogenous-biology coefficients.\n"
            "When unsure, leave it out - a short high-precision list is the goal, and the "
            "author will add anything you missed.\n"
            'Return JSON {"disease_drivers": [...], "druggable": [...], "calibratable": '
            "[...]} using the exact parameter names from THIS batch.\n\n"
            f"PARAMETERS:\n{json.dumps(compact)}\n\nRULES (context):\n" + rule_ctx)
        _merge_roles(out, _parse_json(call(_SYS, user)))
    return out


_VPOP_SYS = (
    "You are a QSP-model analyst. From a model's full parameter list you choose the "
    "subset to VARY across a virtual population: the parameters that plausibly differ "
    "between patients and shape disease severity - typically cell growth / turnover / "
    "migration rates, mediator secretion rates, and pathway up/down-regulation factors. "
    "You EXCLUDE parameters that are structurally fixed across patients - clearance, "
    "decay, and binding/potency constants fixed from data, and pure PK disposition. You "
    "infer each parameter's mechanistic category from its NAME alone, for THIS model - do "
    "not assume any particular disease's vocabulary. Use the EXACT names given. Output "
    "JSON only, no prose.")


def propose_vpop_set(params: list[dict], call, batch: int = 120) -> dict:
    """Select the virtual-population VARIED set from the full model parameter list - the
    paper's category-based approach (a broad set grouped by mechanistic class, not a
    strict handful like ``classify_parameters``). ``params`` is [{name, value, units?},
    ...] (e.g. from ``SimBiologyEngine.list_parameters``). Returns {selected:[names],
    categories:{name: inferred category}, rationale, n_candidates, n_selected}. General:
    categories are inferred per-model from the names, with no fixed disease vocabulary."""
    named = [p for p in params if p.get("name")]
    selected: list[str] = []
    categories: dict[str, str] = {}
    rationale = ""
    for k in range(0, len(named), batch):
        chunk = named[k:k + batch]
        compact = [{"name": p.get("name"), "units": p.get("units"),
                    "value": p.get("value")} for p in chunk]
        batch_names = {p["name"] for p in compact}
        user = (
            "From these model parameters, select the ones to VARY across a virtual "
            "population (patient-to-patient variability that drives severity), grouped by "
            "the mechanistic category you infer from each NAME. Exclude fixed constants "
            "(clearance / decay / binding / PK disposition).\n"
            'Return JSON {"selected": [names], "categories": {name: short category}, '
            '"rationale": "one sentence on your grouping"} using EXACT names from THIS '
            f"batch only.\n\nPARAMETERS:\n{json.dumps(compact)}")
        try:
            d = _parse_json(call(_VPOP_SYS, user))
        except Exception:
            continue
        for n in (d.get("selected") or []):
            if n in batch_names and n not in selected:
                selected.append(n)
        for n, c in (d.get("categories") or {}).items():
            if n in batch_names:
                categories[n] = c
        if not rationale and d.get("rationale"):
            rationale = str(d["rationale"])
    return {"selected": selected, "categories": categories, "rationale": rationale,
            "n_candidates": len(named), "n_selected": len(selected)}


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
        "_review": ("readout_states is a structural draft; vpop_drivers / design_targets "
                    "/ fit_params are high-recall CANDIDATE POOLS - prune to the knobs "
                    "your task actually exposes (which params to expose is a modeler's "
                    "choice, not derivable from the model)."),
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
