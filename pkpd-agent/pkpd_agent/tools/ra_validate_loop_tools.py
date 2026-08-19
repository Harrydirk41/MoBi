"""LLM-loop tools for the RA Stage-6 task: reproduce the paper's HELD-OUT validation.

The paper validated its model by predicting tocilizumab's response in a REFRACTORY
subpopulation that inadequately responded to prior therapy, and comparing to a real
trial. That population is not a shipped model flag - it must be constructed by
running each prior therapy and intersecting the non-responders. The agent designs
that two-stage selection; the harness runs the arms and classifies.

  * ``validate_inspect`` (observe) - the goal, the available therapy arms, the
    readout, and the standard inadequate-responder convention.
  * ``validate_run``     (act)     - given the arms and the IR criteria, run the
    prior-therapy arms, intersect the inadequate responders, run TCZ, and return the
    TCZ response IN that refractory subgroup vs the real comparator.
  * ``validate_finalize``(evaluate)- commit the validation design to be scored.
"""

from __future__ import annotations

from typing import Any

from ..engines.simbiology import SimBiologyEngine
from ..engines import osp_ra_trial
from ..engines import ra_clinical_reference as ra_clin
from .registry import Tool, ToolRegistry, ToolResult


def register_ra_validate_loop_tools(registry: ToolRegistry, config, ctx: dict) -> None:
    """ctx: {sb, vpop, limit, comparator (dict ACR20/50/70)}."""
    sb: SimBiologyEngine = ctx["sb"]
    vpop: str = ctx["vpop"]
    limit: int = int(ctx.get("limit") or 60)
    comparator: dict = ctx.get("comparator") or ra_clin.REFRACTORY_TCZ

    def inspect(args: dict, session) -> ToolResult:
        return ToolResult.success(
            "RA held-out validation: predict tocilizumab's response in the "
            "REFRACTORY subpopulation that inadequately responds to prior therapies, "
            "and compare to a real trial.",
            objective="build the inadequate-responder population by running the prior "
                      "therapies and intersecting their non-responders, then give TCZ "
                      "and read its response in that subgroup",
            available_arms={
                "prior_therapies": ["MTX_15mg_Q1W_SC_t200", "ADA40mg_Q2W_SC_t200"],
                "tcz_arm": "MTX_15mg_Q1W_SC_t200;TCZ8mgkg_Q4W_IV_t200@285"},
            inadequate_responder_convention=(
                "clinically, an inadequate responder has ACR<50 (did not reach ACR50) "
                "AND active disease DAS28-CRP>3.2 after therapy; you may set the "
                "criteria"),
            comparator=f"real trial {comparator.get('trial')}",
            edit_spec_help={
                "prior_arms": "[dose, ...] - each therapy whose non-responders define "
                              "the refractory population (intersection of failures), "
                              "e.g. ['MTX_15mg_Q1W_SC_t200','ADA40mg_Q2W_SC_t200']",
                "tcz_arm": "the TCZ arm to read the refractory response from",
                "acr_key": "ACR level for non-response, 'ACR50' (default) or 'ACR20'",
                "das_threshold": "DAS28-CRP above which disease is active (default 3.2)"},
        )

    def run(args: dict, session) -> ToolResult:
        prior = args.get("prior_arms") or []
        tcz_arm = args.get("tcz_arm") or "MTX_15mg_Q1W_SC_t200;TCZ8mgkg_Q4W_IV_t200@285"
        acr_key = args.get("acr_key") or "ACR50"
        das_thr = float(args.get("das_threshold") or 3.2)
        if len(prior) < 1:
            return ToolResult.error(
                "give prior_arms - the therapies whose non-responders define the "
                "refractory population (e.g. MTX then ADA).")

        # run each prior therapy, classify inadequate responders, intersect
        masks = []
        counts = {}
        for dose in prior:
            r = sb.run_vpop(vpop, dose=dose, stop_time=400.0, readout_day=284.0,
                            limit=limit)
            m = osp_ra_trial.ir_mask(r, acr_key=acr_key, das_threshold=das_thr)
            masks.append(m)
            counts[dose] = sum(1 for v in m.values() if v)
        # intersection of failures across all prior therapies
        common = set.intersection(*[set(m) for m in masks]) if masks else set()
        refractory = {p for p in common if all(m.get(p) for m in masks)}

        # TCZ arm, ACR read at the day-600 second-line readout
        tcz = sb.run_vpop(vpop, dose=tcz_arm, stop_time=700.0, readout_day=600.0,
                          limit=limit)
        resp = osp_ra_trial.response_in_subgroup(tcz, refractory)
        pred = {k: resp.get(k) for k in ("ACR20", "ACR50", "ACR70")}
        score = osp_ra_trial.score_flagship(pred, comparator)
        mae = score.get("mae_pp")

        hist = session.get("val_history") or []
        entry = {"prior_arms": prior, "tcz_arm": tcz_arm, "acr_key": acr_key,
                 "das_threshold": das_thr, "n_refractory": len(refractory),
                 "predicted": pred, "mae": mae}
        hist.append(entry)
        session.put("val_history", hist)

        return ToolResult.success(
            f"IR per arm {counts}; refractory (failed all) n={len(refractory)}; "
            f"TCZ response there ACR20 {pred.get('ACR20')}% ACR50 {pred.get('ACR50')}% "
            f"ACR70 {pred.get('ACR70')}% vs real {comparator.get('trial','').split('(')[0]}"
            f"(ACR20 {comparator.get('ACR20')}) -> MAE {mae} pp",
            per_arm_IR=counts, n_refractory=len(refractory), predicted=pred,
            comparator=comparator, mae_pp=mae, per_endpoint=score.get("per_endpoint"),
            iteration=len(hist))

    def finalize(args: dict, session) -> ToolResult:
        hist = session.get("val_history") or []
        if not hist:
            return ToolResult.error("nothing to finalize - run validate_run first.")
        # match the design if given, else the last run
        prior = args.get("prior_arms")
        match = None
        if prior is not None:
            match = next((h for h in reversed(hist) if h.get("prior_arms") == prior), None)
        match = match or hist[-1]
        session.put("val_final", match)
        return ToolResult.success(
            f"committed validation: refractory n={match['n_refractory']}, TCZ "
            f"ACR20/50/70 {match['predicted']}, MAE {match['mae']} pp vs "
            f"{comparator.get('trial')}",
            **match)

    registry.register(Tool(
        name="validate_inspect",
        description=(
            "OBSERVE the held-out validation task: the goal (predict TCZ in the "
            "refractory population and compare to a real trial), the available "
            "prior-therapy and TCZ arms, and the inadequate-responder convention. "
            "Call this first."),
        input_schema={"type": "object", "properties": {}},
        handler=inspect, phase="observe"))

    registry.register(Tool(
        name="validate_run",
        description=(
            "ACT: run the prior therapies, classify each arm's inadequate responders "
            "(ACR below the level AND DAS28-CRP still active), INTERSECT them to build "
            "the refractory population, then run TCZ and return its response in that "
            "subgroup vs the real comparator trial (MAE). Design the selection: which "
            "prior therapies define refractoriness, and the IR criteria."),
        input_schema={"type": "object", "properties": {
            "prior_arms": {"type": "array", "items": {"type": "string"},
                           "description": "prior-therapy doses; refractory = failed all"},
            "tcz_arm": {"type": "string", "description": "the TCZ arm to read"},
            "acr_key": {"type": "string", "description": "'ACR50' (default) or 'ACR20'"},
            "das_threshold": {"type": "number", "description": "active-disease DAS28 (3.2)"}}},
        handler=run, phase="act"))

    registry.register(Tool(
        name="validate_finalize",
        description=(
            "COMMIT the validation design you recommend (already run). Scores the "
            "refractory TCZ response against the real comparator trial. Call once "
            "before finishing."),
        input_schema={"type": "object", "properties": {
            "prior_arms": {"type": "array", "items": {"type": "string"}}}},
        handler=finalize, phase="evaluate"))
