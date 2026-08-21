"""Model-agnostic loop tools for the held-out VALIDATION task.

The paper validated its model by predicting a therapy's response in a REFRACTORY
subpopulation that inadequately responded to prior therapy, and comparing to a real
trial. That population is not a shipped model flag - it must be constructed by running
each prior therapy and intersecting the non-responders. The agent designs that
two-stage selection; the harness runs the arms and classifies. The available arms and
the real comparator come off the QSPTaskConfig.

  * ``validate_inspect`` (observe) - the goal, the available arms, the IR convention.
  * ``validate_run``     (act)     - run prior arms, intersect IRs, run the test arm,
    return its response in the refractory subgroup vs the real comparator.
  * ``validate_finalize``(evaluate)- commit the validation design to be scored.
"""

from __future__ import annotations

from ..engines.simbiology import SimBiologyEngine
from ..engines import qsp_tasks
from ..engines.qsp_config import QSPTaskConfig
from .registry import Tool, ToolRegistry, ToolResult

_ENDPOINTS = ("ACR20", "ACR50", "ACR70")


def register_qsp_validate_loop_tools(registry: ToolRegistry, config, ctx: dict) -> None:
    """ctx: {cfg, sb, vpop, limit, comparator?}."""
    cfg: QSPTaskConfig = ctx["cfg"]
    sb: SimBiologyEngine = ctx["sb"]
    vpop: str = ctx["vpop"]
    limit: int = int(ctx.get("limit") or 60)
    comparator: dict = ctx.get("comparator") or cfg.refractory_target
    arms = cfg.validate_arms
    default_test_arm = arms.get("tcz_arm", "")
    second_readout = cfg.timeline.get("second_line_readout_day", 600.0)
    first_readout = cfg.timeline.get("first_line_readout_day", 284.0)
    comp_endpoints = [k for k in _ENDPOINTS if k in comparator] \
        or [k for k in comparator if str(k).upper().startswith("ACR")]

    def inspect(args: dict, session) -> ToolResult:
        return ToolResult.success(
            f"{cfg.name} held-out validation: predict the test therapy's response in "
            "the REFRACTORY subpopulation that inadequately responds to prior "
            "therapies, and compare to a real trial.",
            objective="build the inadequate-responder population by running the prior "
                      "therapies and intersecting their non-responders, then give the "
                      "test therapy and read its response in that subgroup",
            available_arms=arms,
            inadequate_responder_convention=(
                "clinically, an inadequate responder did not reach the ACR level AND "
                "still has active disease after therapy; you may set the criteria"),
            comparator=f"real trial {comparator.get('trial')}",
            edit_spec_help={
                "prior_arms": "[dose, ...] - each therapy whose non-responders define "
                              "the refractory population (intersection of failures)",
                "test_arm": "the arm to read the refractory response from",
                "acr_key": "first-line role for non-response, e.g. 'ACR50' or 'ACR20'",
                "das_threshold": "severity above which disease is active (default 3.2)"})

    def run(args: dict, session) -> ToolResult:
        prior = args.get("prior_arms") or arms.get("prior_therapies") or []
        test_arm = args.get("test_arm") or args.get("tcz_arm") or default_test_arm
        acr_key = args.get("acr_key") or "ACR50"
        das_thr = float(args.get("das_threshold") or 3.2)
        if len(prior) < 1:
            return ToolResult.error("give prior_arms - the therapies whose non-responders "
                                    "define the refractory population.")

        masks = []
        counts = {}
        for dose in prior:
            r = sb.run_vpop(vpop, dose=dose, stop_time=400.0, readout_day=first_readout,
                            limit=limit)
            m = cfg.ir_mask(r, acr_key=acr_key, threshold=das_thr)
            masks.append(m)
            counts[dose] = sum(1 for v in m.values() if v)
        common = set.intersection(*[set(m) for m in masks]) if masks else set()
        refractory = {p for p in common if all(m.get(p) for m in masks)}

        test = sb.run_vpop(vpop, dose=test_arm, stop_time=second_readout + 100.0,
                           readout_day=second_readout, limit=limit)
        resp = cfg.response_in_subgroup(test, refractory)
        pred = {k: resp.get(k) for k in comp_endpoints}
        score = qsp_tasks.score_flagship(pred, comparator)
        mae = score.get("mae_pp")

        hist = session.get("val_history") or []
        hist.append({"prior_arms": prior, "test_arm": test_arm, "acr_key": acr_key,
                     "das_threshold": das_thr, "n_refractory": len(refractory),
                     "predicted": pred, "mae": mae})
        session.put("val_history", hist)

        return ToolResult.success(
            f"IR per arm {counts}; refractory (failed all) n={len(refractory)}; test "
            f"response there {pred} vs real "
            f"{str(comparator.get('trial', '')).split('(')[0]}-> MAE {mae} pp",
            per_arm_IR=counts, n_refractory=len(refractory), predicted=pred,
            comparator=comparator, mae_pp=mae, per_endpoint=score.get("per_endpoint"),
            iteration=len(hist))

    def finalize(args: dict, session) -> ToolResult:
        hist = session.get("val_history") or []
        if not hist:
            return ToolResult.error("nothing to finalize - run validate_run first.")
        prior = args.get("prior_arms")
        match = None
        if prior is not None:
            match = next((h for h in reversed(hist) if h.get("prior_arms") == prior), None)
        match = match or hist[-1]
        session.put("val_final", match)
        return ToolResult.success(
            f"committed validation: refractory n={match['n_refractory']}, test response "
            f"{match['predicted']}, MAE {match['mae']} pp vs {comparator.get('trial')}",
            **match)

    registry.register(Tool(
        name="validate_inspect",
        description=("OBSERVE the held-out validation task: the goal (predict the test "
                     "therapy in the refractory population and compare to a real trial), "
                     "the available prior-therapy and test arms, and the IR convention. "
                     "Call this first."),
        input_schema={"type": "object", "properties": {}},
        handler=inspect, phase="observe"))
    registry.register(Tool(
        name="validate_run",
        description=("ACT: run the prior therapies, classify each arm's inadequate "
                     "responders (below the ACR level AND still active), INTERSECT them "
                     "to build the refractory population, then run the test arm and "
                     "return its response in that subgroup vs the real comparator (MAE). "
                     "Design the selection: which prior therapies, and the IR criteria."),
        input_schema={"type": "object", "properties": {
            "prior_arms": {"type": "array", "items": {"type": "string"}},
            "test_arm": {"type": "string"}, "acr_key": {"type": "string"},
            "das_threshold": {"type": "number"}}},
        handler=run, phase="act"))
    registry.register(Tool(
        name="validate_finalize",
        description=("COMMIT the validation design you recommend (already run). Scores "
                     "the refractory test response against the real comparator trial. "
                     "Call once."),
        input_schema={"type": "object", "properties": {
            "prior_arms": {"type": "array", "items": {"type": "string"}}}},
        handler=finalize, phase="evaluate"))
