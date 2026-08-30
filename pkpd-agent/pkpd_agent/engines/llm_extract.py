"""Data-acquisition benchmark (A): can an LLM READ a cited paper and EXTRACT the value the
modellers used? Given a parameter's provenance (reference, figure, experiment description),
a web-enabled LLM searches for the paper and returns the number; we grade it against MOESM2's
"value from reference". Separates the two failure modes the modelling bottleneck hinges on:
retrieval (did it find/read the paper) and extraction (did it get the right number - often
buried in a figure image, not the text).

Pure here (prompt building + grading); the web-enabled call is injected by the runner.
"""

from __future__ import annotations

from .llm_structure import _parse_json

_SYS = ("You are extracting a specific quantitative value from a cited biology paper, to "
        "reproduce a QSP model parameter. Search for the exact paper, read it, and report the "
        "number the modellers would have used. Be honest: if the value is only in a figure "
        "image you cannot read, or you cannot find the paper, say so - a wrong guess is worse "
        "than admitting you could not get it. Output JSON only.")


def build_prompt(prov: dict) -> str:
    """A retrieval+extraction task from a parameter's provenance record."""
    return (
        f"Parameter: {prov.get('name')} ({prov.get('units') or '?'})\n"
        f"Cited reference: {prov.get('reference') or '(none)'}\n"
        f"Figure/Table: {prov.get('figure') or '(unspecified)'}\n"
        f"Experiment: {prov.get('experiment') or '(not described)'}\n"
        f"Meaning: {prov.get('description') or ''}\n\n"
        "Find this paper, read the cited figure/table, and report the numeric value the "
        'modellers extracted (a fold-change, rate, or concentration). Return JSON '
        '{"found_paper": true|false, "value": number|null, "in_figure_only": true|false, '
        '"note": "one phrase on what you could/could not read"}.')


def extract_value(prov: dict, call) -> dict:
    """Run the retrieval+extraction with a web-enabled ``call(system, user)->str``. Returns
    {found_paper, value, in_figure_only, note} (value may be None if not extractable)."""
    d = _parse_json(call(_SYS, build_prompt(prov)))
    if not isinstance(d, dict):
        return {"found_paper": False, "value": None, "in_figure_only": None, "note": "parse fail"}
    v = d.get("value")
    try:
        v = float(v) if v is not None else None
    except (TypeError, ValueError):
        v = None
    return {"found_paper": bool(d.get("found_paper")), "value": v,
            "in_figure_only": d.get("in_figure_only"), "note": d.get("note")}


def extract_value_vision(prov: dict, image_paths: list, call) -> dict:
    """Extraction given the actual FIGURE IMAGE(S): pass the images to a multimodal
    ``call(system, user, image_paths)->str`` and read the value off the chart. Tests the true
    capability once the modality gap (figure vs text) is removed."""
    user = build_prompt(prov) + ("\n\nThe cited figure image(s) are attached above. Read the "
                                 "value directly from the chart/table.")
    d = _parse_json(call(_SYS, user, image_paths))
    if not isinstance(d, dict):
        return {"found_paper": True, "value": None, "in_figure_only": True, "note": "parse fail"}
    v = d.get("value")
    try:
        v = float(v) if v is not None else None
    except (TypeError, ValueError):
        v = None
    return {"found_paper": True, "value": v, "in_figure_only": d.get("in_figure_only"),
            "note": d.get("note")}


def grade(pred_value, truth_value, tol: float = 0.25) -> dict:
    """Compare an extracted value to the reference value: a hit if within ``tol`` relative
    error. Returns {extracted, hit, rel_err}."""
    if pred_value is None or truth_value in (None, 0):
        return {"extracted": pred_value is not None, "hit": False, "rel_err": None}
    rel = abs(pred_value - truth_value) / abs(truth_value)
    return {"extracted": True, "hit": rel <= tol, "rel_err": round(rel, 3)}
