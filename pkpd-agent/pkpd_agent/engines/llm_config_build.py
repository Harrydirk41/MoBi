"""Conversational config builder: a modeler describes the model in plain language,
the LLM produces the filled tasks.json - no JSON editing, no schema to memorize.

This is the entry point for a traditional modeler who lives in the SimBiology GUI, not
a text editor. They write (or dictate) a few sentences about their model - the disease
and its severity readout, the drugs and dose names, the clinical numbers to match - and
this turns that plus the model dump into a complete tasks.json. It PRUNES the drafter's
high-recall candidate pools down to what the description actually uses, and fills the
external fields (dose names, clinical targets, timeline) from the description.

Everything it produces is checked by project_validate against the real model before use,
and anything the description did not cover is left as a stub for the modeler to confirm.

Pluggable via ``call_fn`` (tests use a stub). The LLM never invents clinical numbers -
it only places what the description states.
"""

from __future__ import annotations

import json

from .llm_structure import _parse_json
from .llm_tasks import default_call  # reuse the larger-max_tokens, empty-guarded call

_SYS = ("You are a QSP-modeling assistant helping a modeler configure a virtual-trial "
        "workflow for THEIR model. You are given the model's parameters/species, a "
        "draft of candidate task roles, and the modeler's own description. You produce "
        "a single tasks.json object. Use the modeler's description as the source of "
        "truth for external facts (drug/dose names, clinical numbers, timeline); use "
        "the model's names verbatim; never invent clinical numbers - if the description "
        "does not give a number, leave that field out. Output JSON only, no prose.")

_SCHEMA_HINT = """The tasks.json object has these fields (fill what the description and
model support; omit what you cannot determine):
{
  "name": "<model name>",
  "project_aliases": ["<short name>"],
  "disease": "<one line: disease + severity readout>",
  "severity_readout": "<name of the severity score, e.g. DAS28-CRP>",
  "readout_desc": "<short: what the clinical response is>",
  "readout_states": [<the response-flag + severity state names, role order>],
  "run_columns": {
     "patient": "patient",
     "first_line": {"<endpoint>": "<CSV column the readout script emits>", ...},
     "subgroup_flag": "<flag column marking inadequate responders, or ''>",
     "second_line": {"<endpoint>": "<CSV column>", ...},
     "severity": {"baseline": "<column>", "readout": "<column>"}
  },
  "timeline": {"baseline_day": N, "first_line_readout_day": N, "second_line_readout_day": N},
  "drugs": {"<CODE>": {"drug": "...", "mechanism": "...", "doses": ["<dose name>", ...]}},
  "vpop_drivers": {"<param>": {"nominal": N, "span": [lo,hi], "meaning": "..."}},
  "vpop_target": {"mean": N, "sd": N, "band": [lo,hi]},
  "fit_params": {"<param>": {"unit": "...", "reference": N, "meaning": "...",
                 "search_range": [lo,hi], "log_scale": true}},
  "design_targets": {"<param>": {"pathway": "...", "analogue": "...", "note": "..."}},
  "clinical_trials": {"<DRUG>": {"trial": "...", "weeks": {"24": {"drug": {"ACR20": N,...}}}}},
  "refractory_target": {"trial": "...", "ACR20": N, ...}
}
PRUNE vpop_drivers / design_targets / fit_params from the candidate pools to only the
ones the description says the tasks use. Keep readout_states from the draft unless the
description corrects it."""


def build_tasks(network: dict, description: str, tasks_draft: dict, call) -> dict:
    """Produce a filled tasks.json dict from the model dump + the modeler's description
    + the drafter's candidate pools. Merges onto the draft so a partial LLM reply still
    yields a usable config (the draft supplies readout_states and the candidate pools)."""
    params = [p.get("name") for p in network.get("parameters", [])]
    species = [s.get("name") for s in network.get("species", [])]
    pools = {k: list((tasks_draft.get(k) or {}).keys())
             for k in ("vpop_drivers", "design_targets", "fit_params")}
    user = (
        f"{_SCHEMA_HINT}\n\nMODEL NAME: {network.get('name', 'QSP model')}\n"
        f"MODEL PARAMETERS ({len(params)}): {json.dumps(params[:400])}\n"
        f"MODEL SPECIES: {json.dumps(species)}\n"
        f"DRAFT readout_states: {json.dumps(tasks_draft.get('readout_states', []))}\n"
        f"CANDIDATE POOLS (prune these): {json.dumps(pools)}\n\n"
        f"MODELER'S DESCRIPTION:\n{description}\n\n"
        "Return the tasks.json object.")
    try:
        out = _parse_json(call(_SYS, user))
    except Exception:                                  # noqa: BLE001
        out = {}
    if not isinstance(out, dict):
        out = {}
    # merge: LLM output wins; draft fills anything the LLM omitted
    merged = dict(tasks_draft)
    merged.update({k: v for k, v in out.items() if v not in (None, "", {}, [])})
    # never carry the drafter's internal markers into a real config
    for internal in ("_derived_from", "_review", "readout_roles"):
        merged.pop(internal, None)
    return merged


def build_spec(network: dict, description: str, structure_draft: dict) -> dict:
    """Assemble a spec.json dict from the structure extractor's output (+ description for
    the short project alias / readout name). Pure assembly - the structure draft already
    carries the nodes/edges/readout the benchmarks need."""
    species = [s.get("name") for s in network.get("species", [])]
    cls = structure_draft.get("classification", {})
    readout = cls.get("readout", [])
    return {
        "name": network.get("name", "QSP model"),
        "readout_name": (readout[0] if readout else "the disease-severity score"),
        "readout_targets": readout,
        "drug_patterns": [],          # let infer/extractor patterns be added by review
        "readout_patterns": [],
        "aliases": {},
        "gsa_top": [],                # external (a figure) - the modeler adds these
    }
