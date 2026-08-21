"""Conversational onboarding tools: a modeler describes the model, the agent builds the
project config - no JSON editing, no schema, no hand-run CLI steps.

The agent gets four tools and drives the whole build itself: inspect the model, build a
config from the modeler's plain-language description, correct any specific field the
validator flags (or the modeler mentions), and save the project once the validator is
clean. It is the GUI-free entry point - the modeler supplies words, the agent supplies
the config and reports back in plain English.

  * ``onboard_inspect`` (observe) - the model's species / parameter counts and a few
    names, so the agent (and modeler) can talk about it.
  * ``onboard_build``   (act)     - extract structure + draft roles + fill from the
    description; validate against the model; return the plain-English report.
  * ``onboard_set``     (act)     - set ONE field (dotted path) the validator flagged
    or the modeler corrected; re-validate.
  * ``onboard_save``    (evaluate)- write projects/<name>/{spec,tasks}.json (only when
    there are no validation ERRORS).

ctx: {network (dict), description (str), name (str), call (LLM fn), out_dir?}.
"""

from __future__ import annotations

import json
import os

from ..engines import llm_structure as LS
from ..engines import llm_tasks as LT
from ..engines import llm_config_build as LC
from ..engines import project_validate as PV
from .registry import Tool, ToolRegistry, ToolResult


def _set_path(d: dict, path: str, value):
    """Set a dotted path (e.g. 'vpop_target.band') in a nested dict, creating dicts."""
    keys = path.split(".")
    cur = d
    for k in keys[:-1]:
        cur = cur.setdefault(k, {})
        if not isinstance(cur, dict):
            raise ValueError(f"'{path}' passes through a non-object at '{k}'")
    cur[keys[-1]] = value


def _summary(tasks: dict) -> dict:
    def n(f):
        v = tasks.get(f)
        return len(v) if isinstance(v, (dict, list)) else 0
    stub = [f for f in ("drugs", "vpop_target", "clinical_trials", "refractory_target")
            if not tasks.get(f) or (isinstance(tasks.get(f), str)
                                    and str(tasks.get(f)).startswith("TODO"))]
    return {"readout_states": n("readout_states"), "vpop_drivers": n("vpop_drivers"),
            "design_targets": n("design_targets"), "fit_params": n("fit_params"),
            "still_stubbed": stub}


def register_onboarding_tools(registry: ToolRegistry, config, ctx: dict) -> None:
    network: dict = ctx["network"]
    description: str = ctx.get("description", "")
    name: str = ctx["name"]
    call = ctx["call"]
    out_dir = ctx.get("out_dir") or os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..", "..", "projects"))

    def inspect(a, s):
        sp = [x.get("name") for x in network.get("species", [])]
        pr = [x.get("name") for x in network.get("parameters", [])]
        return ToolResult.success(
            f"Model '{network.get('name', name)}': {len(sp)} species, {len(pr)} "
            "parameters. Build a config from the modeler's description, then validate.",
            n_species=len(sp), n_parameters=len(pr),
            example_species=sp[:12], example_parameters=pr[:12],
            have_description=bool(description.strip()))

    def build(a, s):
        desc = a.get("description") or description
        structure = LS.extract_structure(network, call)
        draft = LT.draft_tasks(network, call, name=network.get("name", name))
        tasks = LC.build_tasks(network, desc, draft, call)
        spec = LC.build_spec(network, desc, structure)
        report = PV.validate_project(tasks, spec, network)
        s.put("onboard_tasks", tasks)
        s.put("onboard_spec", spec)
        s.put("onboard_report", report)
        return ToolResult.success(
            "built config from the description. " + PV.format_report(report),
            summary=_summary(tasks), errors=report["errors"],
            warnings=report["warnings"])

    def set_field(a, s):
        tasks = s.get("onboard_tasks")
        if tasks is None:
            return ToolResult.error("build the config first (onboard_build).")
        path, value = a.get("path"), a.get("value")
        if not path:
            return ToolResult.error("give {path, value}, e.g. path 'vpop_target.band' "
                                    "value [6, 20].")
        target = s.get("onboard_spec") if path.startswith("spec.") else tasks
        p = path[5:] if path.startswith("spec.") else path
        try:
            _set_path(target, p, value)
        except ValueError as e:
            return ToolResult.error(str(e))
        report = PV.validate_project(tasks, s.get("onboard_spec"), network)
        s.put("onboard_report", report)
        return ToolResult.success(f"set {path}. " + PV.format_report(report),
                                  errors=report["errors"], warnings=report["warnings"])

    def save(a, s):
        tasks, spec = s.get("onboard_tasks"), s.get("onboard_spec")
        report = s.get("onboard_report") or {}
        if tasks is None:
            return ToolResult.error("nothing to save - build the config first.")
        if report.get("errors"):
            return ToolResult.error(
                "config still has ERRORS - fix them with onboard_set before saving:\n"
                + PV.format_report(report))
        folder = os.path.join(out_dir, name)
        os.makedirs(folder, exist_ok=True)
        with open(os.path.join(folder, "tasks.json"), "w", encoding="utf-8") as fh:
            json.dump(tasks, fh, indent=2)
        with open(os.path.join(folder, "spec.json"), "w", encoding="utf-8") as fh:
            json.dump(spec, fh, indent=2)
        s.put("onboard_saved", folder)
        return ToolResult.success(
            f"saved projects/{name}/tasks.json + spec.json. Remaining warnings (fill "
            f"when you have the data): {report.get('warnings')}", folder=folder)

    registry.register(Tool(
        name="onboard_inspect",
        description="OBSERVE the model: species / parameter counts and a few names. Call first.",
        input_schema={"type": "object", "properties": {}}, handler=inspect, phase="observe"))
    registry.register(Tool(
        name="onboard_build",
        description=("ACT: build the project config (spec + tasks) from the modeler's "
                     "plain-language description of the disease, drugs, clinical numbers "
                     "and parameters. Returns a plain-English validation report. Pass an "
                     "updated 'description' to rebuild with more detail."),
        input_schema={"type": "object", "properties": {
            "description": {"type": "string"}}}, handler=build, phase="act"))
    registry.register(Tool(
        name="onboard_set",
        description=("ACT: set ONE config field the validator flagged or the modeler "
                     "corrected, by dotted path. e.g. path 'vpop_target.band' value "
                     "[6,20]; path 'spec.gsa_top' value [...]. Re-validates."),
        input_schema={"type": "object", "properties": {
            "path": {"type": "string"}, "value": {}}}, handler=set_field, phase="act"))
    registry.register(Tool(
        name="onboard_save",
        description=("COMMIT: write projects/<name>/{spec,tasks}.json. Refuses while "
                     "there are validation ERRORS (fix with onboard_set first)."),
        input_schema={"type": "object", "properties": {}}, handler=save, phase="evaluate"))
