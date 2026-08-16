"""LLM-loop tools for the DDI benchmark: inspect the interaction + try-a-model.

The DDI analogue of osp_loop_tools. Where the single-compound loop tunes drug
disposition, the DDI loop tunes the perpetrator's INTERACTION parameters (Ki,
kinact/K_kinact_half, EC50/Emax) so the model reproduces the observed change in
the victim's exposure (AUCR = AUC_treatment / AUC_control):

  * ``ddi_inspect``   (observe) - the interaction structure (perpetrator, victim,
    mechanisms with their current parameters), the control/treatment pairs, the
    observed interaction ratios, and the a-priori identifiability guidance.
  * ``ddi_try_model`` (act)     - apply interaction_parameters, run every
    control/treatment arm headless, and return the predicted vs observed AUCR
    (GMFE on the ratio, fraction within 2-fold).

The victim's own disposition is FIXED (it is a validated model); the unknown is
the interaction. Registered per-task, so the handlers close over the DDI context.
"""

from __future__ import annotations

from typing import Any

from ..engines.osp_cli import OSPCli
from ..engines import osp_ddi
from ..engines import osp_catalog
from .registry import Tool, ToolRegistry, ToolResult


def register_osp_ddi_loop_tools(registry: ToolRegistry, config, ctx: dict) -> None:
    """ctx: {cli, snapshot_path, ddi, victim, observed_ratios, input}."""
    cli: OSPCli = ctx["cli"]
    snapshot_path: str = ctx["snapshot_path"]
    ddi: dict = ctx["ddi"]
    victim: str = ctx["victim"]
    observed_ratios: list = ctx.get("observed_ratios") or []
    inp: dict = ctx.get("input") or {}

    def _mechanism_catalog() -> list[dict]:
        out = []
        for p in ddi.get("perpetrators") or []:
            for m in p.get("mechanisms") or []:
                entry = osp_catalog.interaction_by_internal_name(m.get("internal_name")) or {}
                out.append({
                    "perpetrator": p.get("name"),
                    "internal_name": m.get("internal_name"),
                    "target": m.get("target"),
                    "current_parameters": m.get("parameters"),
                    "parameters": [q.get("name") for q in entry.get("parameters") or []],
                    "identifiability": entry.get("identifiability"),
                    "validated": entry.get("validated"),
                })
        return out

    # -- observe -------------------------------------------------------- #
    def inspect(args: dict, session) -> ToolResult:
        recs = osp_ddi.ddi_identifiability(
            ddi, n_observed_ratios=len(observed_ratios) or None)
        return ToolResult.success(
            "DDI task: perpetrator -> victim interaction, the mechanisms and their "
            "current parameters, the control/treatment pairs, the observed "
            "interaction ratios, and the identifiability guidance",
            objective=inp.get("objective"),
            perpetrators=[p.get("name") for p in ddi.get("perpetrators") or []],
            victim=victim,
            mechanisms=_mechanism_catalog(),
            control_treatment_pairs=ddi.get("pairs"),
            observed_interaction_ratios=observed_ratios,
            unknowns_guidance=inp.get("unknowns_guidance"),
            how_scored=inp.get("how_scored"),
            identifiability_recommendations=recs,
            edit_spec_help={
                "interaction_parameters": "[{perpetrator, internal_name, target, "
                "parameters:{name: value}}] - set/estimate the perpetrator's "
                "interaction parameters; the victim's disposition is fixed"},
        )

    # -- act (run all arms + score the ratio) --------------------------- #
    def try_model(args: dict, session) -> ToolResult:
        ip = args.get("interaction_parameters")
        edits = {"interaction_parameters": ip} if ip else None
        out = osp_ddi.run_ddi_prediction(
            cli, snapshot_path, ddi, victim, edits=edits,
            observed_ratios=observed_ratios)
        if not out.get("ok"):
            return ToolResult.error(f"DDI run failed: {out.get('message')}",
                                    edits_applied=out.get("edits_applied"))
        score = out.get("score") or {}
        gmfe = score.get("gmfe_aucr")
        hist = session.get("ddi_history") or []
        hist.append({"interaction_parameters": ip, "gmfe_aucr": gmfe})
        session.put("ddi_history", hist)
        best = session.get("ddi_best_gmfe")
        if gmfe is not None and (best is None or gmfe < best):
            session.put("ddi_best_gmfe", gmfe)
            session.put("ddi_best_edits", edits)
        return ToolResult.success(
            f"AUCR GMFE {gmfe} (within2fold "
            f"{score.get('within_2fold_pct')}%); best so far "
            f"{session.get('ddi_best_gmfe')}",
            gmfe_aucr=gmfe,
            within_2fold_pct=score.get("within_2fold_pct"),
            per_arm=score.get("per_arm"),
            predicted_ratios=out.get("predicted_ratios"),
            recommendations=out.get("recommendations"),
            best_gmfe_so_far=session.get("ddi_best_gmfe"),
            iteration=len(hist),
        )

    registry.register(Tool(
        name="ddi_inspect",
        description=(
            "OBSERVE the DDI task: the perpetrator(s), the victim, each interaction "
            "mechanism with its parameters (Ki, kinact/K_kinact_half, EC50/Emax) and "
            "current values, the control/treatment simulation pairs, the observed "
            "interaction ratios (AUCR), and a-priori identifiability guidance. Call "
            "this first."),
        input_schema={"type": "object", "properties": {}},
        handler=inspect, phase="observe"))

    registry.register(Tool(
        name="ddi_try_model",
        description=(
            "ACT: set the perpetrator's interaction parameters and run every "
            "control/treatment arm headless, returning the predicted vs observed "
            "interaction ratio (AUCR = AUC_treatment/AUC_control): GMFE on the "
            "ratio, fraction within 2-fold, and per-arm fold error. Do not change "
            "the victim's disposition - the interaction is the unknown."),
        input_schema={"type": "object", "properties": {
            "interaction_parameters": {
                "type": "array",
                "description": "[{perpetrator, internal_name, target, parameters}]",
                "items": {"type": "object"}}}},
        handler=try_model, phase="act"))
