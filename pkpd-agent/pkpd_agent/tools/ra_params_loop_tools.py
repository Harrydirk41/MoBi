"""LLM-loop tools for Stage-1 Layer 5: estimate the model's PARAMETER values.

The agent is shown each parameter's name, units, and cell context and must predict its
value from physiology (units + role -> order of magnitude). Scored by order-of-magnitude
error against the model, and - the honest part - against a naive unit-geomean baseline
that knows the scale of each unit but no biology. Values are revealed only at finalize.

  * ``param_inspect``  (observe) - the list of parameters to estimate (name/units/cell),
    the objective, and how it is scored. No values, no baseline.
  * ``param_estimate`` (act)     - submit {name, value} guesses; accumulate across calls.
  * ``param_finalize`` (evaluate)- score all guesses (order-of-magnitude, split by
    dimensionless vs dimensional, vs the naive baseline). Terminal.
"""

from __future__ import annotations

from ..engines import ra_params as P
from .registry import Tool, ToolRegistry, ToolResult


def register_ra_params_loop_tools(registry: ToolRegistry, config, ctx: dict) -> None:
    """ctx: {truth: list[Param]} - the answer key, never shown."""
    truth: list = ctx["truth"]

    def inspect(args: dict, session) -> ToolResult:
        view = P.prompt_view(truth)
        n_dim = sum(1 for p in truth if not p.dimensionless())
        return ToolResult.success(
            f"RA parameter-estimation task: predict the value of each of {len(truth)} "
            "model parameters from its name, units, and cell context. This is the "
            "quantitative layer - set a physiologically grounded prior for each.",
            objective="predict a numeric value for every parameter; you are scored on "
                      "ORDER-OF-MAGNITUDE accuracy (log10 error) vs the real model",
            scoring=("median log10 error and fraction within 3x / 10x, reported "
                     "separately for dimensionless fold-effects (cluster near 1, easy) "
                     f"and the {n_dim} DIMENSIONAL parameters (rates 1/day, secretion "
                     "ng/molecule/day, concentrations M) where physiology actually helps"),
            guidance=[
                "For '1/day' rates: think biological turnover - cell death/proliferation "
                "~0.01-1 /day, cytokine clearance faster.",
                "For 'dimensionless' Max-fold effects: the model keeps them modest, "
                "typically 0.3-3 (a fold CHANGE, so ~1).",
                "Use the cell context and the source/target in the name (e.g. "
                "'kcl_VEGF' = VEGF clearance rate; 'HalfEffectConc_..._byIL6' = an IL-6 "
                "EC50 concentration)."],
            parameters=view)

    def estimate(args: dict, session) -> ToolResult:
        preds = P.clean_predictions(args.get("predictions") or args.get("values") or [])
        acc = dict(session.get("param_preds") or {})
        acc.update(preds)
        session.put("param_preds", acc)
        remaining = [p.name for p in truth if p.name not in acc]
        return ToolResult.success(
            f"accepted {len(preds)} value(s); {len(acc)}/{len(truth)} parameters "
            f"estimated ({len(remaining)} remaining).",
            n_estimated=len(acc), n_remaining=len(remaining),
            still_missing=remaining[:20])

    def finalize(args: dict, session) -> ToolResult:
        preds = dict(session.get("param_preds") or {})
        if not preds:
            return ToolResult.error("estimate parameter values before finalizing.")
        score = P.score_params(preds, truth)
        session.put("param_final", score)
        ov, base = score["overall"], score["naive_unit_geomean_baseline"]
        dim = score["dimensional"]
        return ToolResult.success(
            f"scored {score['n_scored']}/{score['n_truth']} parameters. "
            f"Overall median log10 err {ov['median_log10_err']} (within 10x "
            f"{ov['within_10x']}). DIMENSIONAL median {dim['median_log10_err']} "
            f"(within 3x {dim['within_3x']}). Naive baseline median "
            f"{base['median_log10_err']}. Beats baseline: {score['beats_baseline']}.",
            overall=ov, dimensional=dim, dimensionless=score["dimensionless"],
            naive_baseline=base, beats_baseline=score["beats_baseline"])

    registry.register(Tool(
        name="param_inspect",
        description=("OBSERVE the parameter-estimation task: the list of parameters "
                     "(name, units, cell context), the objective, and the order-of-"
                     "magnitude scoring. No values are given. Call first."),
        input_schema={"type": "object", "properties": {}},
        handler=inspect, phase="observe"))

    registry.register(Tool(
        name="param_estimate",
        description=("ACT: submit value guesses as a list of {name, value}. Accumulates "
                     "across calls; reports how many parameters remain. Cover them all "
                     "before finalizing."),
        input_schema={"type": "object", "properties": {
            "predictions": {"type": "array", "items": {"type": "object", "properties": {
                "name": {"type": "string"},
                "value": {"type": "number"}}}}}},
        handler=estimate, phase="act"))

    registry.register(Tool(
        name="param_finalize",
        description=("COMMIT all estimates; scored by order-of-magnitude error vs the "
                     "model, split dimensionless/dimensional, against a naive baseline. "
                     "Terminal - call once."),
        input_schema={"type": "object", "properties": {}},
        handler=finalize, phase="evaluate"))
