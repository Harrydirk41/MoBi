"""LLM-loop tools for the OSP PBPK benchmark: inspect + try-a-model.

Two tools give Claude the closed loop:

  * ``osp_inspect``   (observe) - the current model (parameters, distribution
    method, processes), the literature priors, and the observed-data overview.
  * ``osp_try_model`` (act)     - apply an edit spec, run PK-Sim headless, and
    return the fit (GMFE, %-within-2-fold, and per-route BIAS so the agent knows
    which way to move each parameter), plus plausibility flags.

The task context (the base snapshot to edit, the observed data, the PK-Sim CLI)
is captured in the handlers, so these are registered per-task, not globally.
"""

from __future__ import annotations

import json
from typing import Any

from ..engines.osp_cli import OSPCli
from ..engines import osp_score
from ..engines.snapshot_edit import PARTITION_METHODS, PERMEABILITY_METHODS
from .registry import Tool, ToolRegistry, ToolResult


def _current_model(snapshot_path: str) -> dict[str, Any]:
    with open(snapshot_path, encoding="utf-8") as fh:
        comp = (json.load(fh).get("Compounds") or [{}])[0]
    params = []
    seen = set()

    def walk(o):
        if isinstance(o, dict):
            nm = o.get("Name")
            if isinstance(nm, str) and isinstance(o.get("Value"), (int, float)) \
                    and nm not in seen:
                seen.add(nm)
                params.append({"name": nm, "value": o["Value"],
                               "unit": o.get("Unit", "")})
            for v in o.values():
                walk(v)
        elif isinstance(o, list):
            for v in o:
                walk(v)

    walk(comp)
    return {
        "parameters": params,
        "calculation_methods": comp.get("CalculationMethods") or [],
        "processes": [{"molecule": p.get("Molecule"),
                       "internal": p.get("InternalName")}
                      for p in comp.get("Processes") or []],
    }


def _observed_overview(observed: list[dict]) -> dict[str, Any]:
    routes: dict[str, dict] = {}
    for o in observed:
        r = (o.get("route") or "NA")
        routes.setdefault(r, {"n_datasets": 0, "studies": set(), "doses": set()})
        routes[r]["n_datasets"] += 1
        if o.get("study"):
            routes[r]["studies"].add(o["study"])
        if o.get("dose"):
            routes[r]["doses"].add(str(o["dose"]))
    for r in routes.values():
        r["studies"] = sorted(r["studies"])
        r["doses"] = sorted(r["doses"])
    return {"n_datasets": len(observed), "by_route": routes}


def register_osp_loop_tools(registry: ToolRegistry, config, ctx: dict) -> None:
    """ctx: {cli: OSPCli, snapshot_path, observed: [...], input: {...}}."""
    cli: OSPCli = ctx["cli"]
    snapshot_path: str = ctx["snapshot_path"]
    observed: list[dict] = ctx["observed"]
    inp: dict = ctx.get("input") or {}

    # -- observe -------------------------------------------------------- #
    def inspect(args: dict, session) -> ToolResult:
        gd = inp.get("given_data", {}) or {}
        bg = inp.get("background") or {}
        return ToolResult.success(
            "task description: objective, known biology, literature priors, "
            "current model, and the observed clinical data",
            objective=inp.get("objective"),
            background=bg.get("description"),
            known_biology=bg.get("literature_facts"),
            compound_identity=gd.get("compound_identity"),
            literature_physicochemical=gd.get("literature_physicochemical")
            or inp.get("literature_physicochemical"),
            unknowns_guidance=inp.get("unknowns_guidance"),
            evaluation_rubric=inp.get("evaluation_rubric"),
            current_model=_current_model(snapshot_path),
            observed_overview=_observed_overview(observed),
            valid_partition_methods=PARTITION_METHODS,
            valid_permeability_methods=PERMEABILITY_METHODS,
        )

    # -- act (run + score) ---------------------------------------------- #
    def try_model(args: dict, session) -> ToolResult:
        edits = args.get("edits") or {}
        res = cli.build_and_run(snapshot_path, edits=edits)
        if not res["ok"]:
            return ToolResult.error(
                f"PK-Sim run failed: {res['message']}",
                edits_applied=res.get("edits_applied"))

        predicted, unmatched = osp_score.map_predictions(res["profiles"], observed)
        score = osp_score.score_fit(observed, predicted)
        applied = res.get("edits_applied") or {}
        param_list = [{"parameter": k, "value": v, "unit": ""}
                      for k, v in applied.get("parameters", {}).items()]
        flags = osp_score.plausibility(param_list)

        overall = score["overall"]
        # track history + best
        hist = session.get("osp_history") or []
        hist.append({"edits": edits, "gmfe": overall.get("gmfe")})
        session.put("osp_history", hist)
        best = session.get("osp_best_gmfe")
        if overall.get("gmfe") is not None and (best is None or overall["gmfe"] < best):
            session.put("osp_best_gmfe", overall["gmfe"])
            session.put("osp_best_edits", edits)

        worst = [{"dataset": d["dataset"], "route": d["route"],
                  "gmfe": d["gmfe"], "bias": d["bias"]}
                 for d in score["per_dataset"][:3]]
        return ToolResult.success(
            f"GMFE {overall.get('gmfe')} overall "
            f"(within2fold {overall.get('within_2fold_pct')}%); "
            f"best so far {session.get('osp_best_gmfe')}",
            gmfe_overall=overall.get("gmfe"),
            within_2fold_pct=overall.get("within_2fold_pct"),
            bias_overall=overall.get("bias"),
            by_route=score["by_route"],
            worst_datasets=worst,
            parameter_flags=flags,
            edits_applied=applied,
            not_found=applied.get("not_found"),
            n_matched=len(predicted), n_total=len(observed),
            best_gmfe_so_far=session.get("osp_best_gmfe"),
            iteration=len(hist),
        )

    registry.register(Tool(
        name="osp_inspect",
        description=(
            "OBSERVE the PBPK task: returns the current compound model "
            "(parameters with values/units, the distribution & permeability "
            "calculation methods, and the metabolizing/clearance processes), the "
            "literature physicochemical priors, and an overview of the observed "
            "clinical datasets (routes, studies, doses). Call this first."),
        input_schema={"type": "object", "properties": {}},
        handler=inspect, phase="observe"))

    registry.register(Tool(
        name="osp_try_model",
        description=(
            "ACT: apply an edit spec to the model, run PK-Sim headless, and score "
            "the fit against the observed data. edits = {parameters:{name:value}, "
            "calculation_methods:{partition:..,permeability:..}, "
            "processes:{Molecule:true/false}}. Returns GMFE and %-within-2-fold "
            "overall and PER ROUTE, plus the per-route geometric BIAS (>1 = the "
            "model over-predicts -> lower exposure/raise clearance; <1 = "
            "under-predicts). Use the bias and worst datasets to decide the next "
            "edit. Iterate to minimize GMFE."),
        input_schema={
            "type": "object",
            "properties": {
                "edits": {
                    "type": "object",
                    "description": "parameters / calculation_methods / processes",
                    "properties": {
                        "parameters": {"type": "object"},
                        "calculation_methods": {"type": "object"},
                        "processes": {"type": "object"},
                    },
                },
            },
            "required": ["edits"],
        },
        handler=try_model, phase="act"))
