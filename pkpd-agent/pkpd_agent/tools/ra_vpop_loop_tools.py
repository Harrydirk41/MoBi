"""LLM-loop tools for the RA Stage-3 task: GENERATE a virtual population.

The agent chooses which disease-driver parameters to vary and over what bounds,
samples a candidate population, simulates each to its untreated disease baseline,
and matches the resulting baseline DAS28-CRP distribution to the clinical target.
Too-wide bounds push patients out of the active-disease band (low yield, too much
spread); too-narrow bounds kill phenotypic diversity. Finding the bounds that
reproduce a realistic, diverse active-RA population is the task.

  * ``vpop_inspect`` (observe) - the disease-driver parameters (meaning, nominal,
    observed span), the clinical target distribution, and the active band.
  * ``vpop_sample``  (act)     - set per-parameter bounds, generate + simulate a
    candidate population, and return the baseline DAS28-CRP distribution + yield.
  * ``vpop_finalize``(evaluate)- commit the sampling design to be scored.
"""

from __future__ import annotations

from typing import Any

from ..engines.simbiology import SimBiologyEngine
from ..engines import osp_ra_trial
from .registry import Tool, ToolRegistry, ToolResult


def register_ra_vpop_loop_tools(registry: ToolRegistry, config, ctx: dict) -> None:
    """ctx: {sb, n_samples, baseline_day, seed, target}."""
    sb: SimBiologyEngine = ctx["sb"]
    n_samples: int = int(ctx.get("n_samples") or 60)
    baseline_day: float = float(ctx.get("baseline_day") or 200.0)
    seed: int = int(ctx.get("seed") or 1)
    target: dict = ctx.get("target") or osp_ra_trial.VPOP_TARGET

    def _catalog() -> list[dict]:
        out = []
        for name, p in osp_ra_trial.VPOP_PARAMS.items():
            out.append({"name": name, "meaning": p["meaning"],
                        "nominal": p["nominal"], "observed_span": p["span"]})
        return out

    # -- observe -------------------------------------------------------- #
    def inspect(args: dict, session) -> ToolResult:
        return ToolResult.success(
            "RA virtual-population generation: sample disease-driver parameters so "
            "the untreated baseline DAS28-CRP distribution matches an active-RA "
            "clinical target.",
            objective="build a virtual population whose baseline DAS28-CRP matches "
                      f"mean {target['mean']}, sd {target['sd']}, mostly inside the "
                      f"active band {target['band']}",
            disease_driver_parameters=_catalog(),
            clinical_target=target,
            n_candidates_per_sample=n_samples,
            edit_spec_help={
                "bounds": "{param: [lo, hi, scale]} - the sampling range per "
                          "parameter; scale is 'lin' or 'log' (use 'log' for the "
                          "wide-ranging factors). Each candidate draws every listed "
                          "parameter independently. Vary enough drivers, over wide "
                          "enough bounds, to get phenotypic diversity - but not so "
                          "wide that most patients fall outside the active band."},
        )

    # -- act ------------------------------------------------------------ #
    def sample(args: dict, session) -> ToolResult:
        bounds = args.get("bounds") or {}
        if not bounds:
            return ToolResult.error(
                "no bounds given - pass {bounds:{param:[lo,hi,scale]}}. "
                "Call vpop_inspect for the parameter names and spans.")
        spec = osp_ra_trial.build_sample_spec(bounds)
        r = sb.sample_vpop(spec, n_samples=n_samples, baseline_day=baseline_day,
                           seed=seed)
        das = (r.get("columns") or {}).get("DAS28_base", [])
        score = osp_ra_trial.score_vpop(das, target)

        hist = session.get("vpop_history") or []
        hist.append({"bounds": bounds, "score": score})
        session.put("vpop_history", hist)
        best = session.get("vpop_best_dist")
        d = score.get("distribution_distance")
        if d is not None and (best is None or d < best):
            session.put("vpop_best_dist", d)
            session.put("vpop_best_bounds", bounds)

        return ToolResult.success(
            f"sampled {score.get('n')} candidates: yield {score.get('yield_pct')}% "
            f"in band {target['band']}; accepted DAS28 mean "
            f"{score.get('accepted_mean')} sd {score.get('accepted_sd')} "
            f"(target {target['mean']}/{target['sd']}); distance "
            f"{score.get('distribution_distance')} (best {session.get('vpop_best_dist')})",
            bounds=bounds, **score,
            best_distance_so_far=session.get("vpop_best_dist"), iteration=len(hist))

    # -- evaluate ------------------------------------------------------- #
    def finalize(args: dict, session) -> ToolResult:
        bounds = args.get("bounds") or session.get("vpop_best_bounds")
        if not bounds:
            return ToolResult.error("nothing to finalize - run vpop_sample first.")
        hist = session.get("vpop_history") or []
        match = next((h for h in reversed(hist) if h.get("bounds") == bounds), None)
        if match is None:
            return ToolResult.error("that exact bounds set was not sampled - "
                                    "vpop_sample it first, then finalize.")
        session.put("vpop_final", match)
        s = match["score"]
        return ToolResult.success(
            f"committed vpop design: yield {s.get('yield_pct')}%, DAS28 "
            f"{s.get('accepted_mean')}/{s.get('accepted_sd')}, distance "
            f"{s.get('distribution_distance')}",
            bounds=bounds, score=s)

    registry.register(Tool(
        name="vpop_inspect",
        description=(
            "OBSERVE the virtual-population task: the disease-driver parameters "
            "(pro-inflammatory amplification factors and cell-baseline growth rates, "
            "each with its nominal value and observed span), the clinical target "
            "baseline DAS28-CRP distribution, and the active-disease band. Call first."),
        input_schema={"type": "object", "properties": {}},
        handler=inspect, phase="observe"))

    registry.register(Tool(
        name="vpop_sample",
        description=(
            "ACT: set the sampling bounds per disease-driver parameter and generate a "
            "candidate virtual population, simulating each patient to its untreated "
            "baseline. Returns the baseline DAS28-CRP yield (fraction in the active "
            "band), the accepted cohort's mean/sd, and the distance to the target "
            "distribution. Iterate the bounds to raise yield and match mean/sd with "
            "realistic spread."),
        input_schema={"type": "object", "properties": {
            "bounds": {"type": "object",
                       "description": "{param: [lo, hi, scale]}, scale 'lin' or 'log'"}}},
        handler=sample, phase="act"))

    registry.register(Tool(
        name="vpop_finalize",
        description=(
            "COMMIT the sampling design you recommend (already run with vpop_sample). "
            "Scores the generated population against the clinical target distribution. "
            "Call once before finishing."),
        input_schema={"type": "object", "properties": {
            "bounds": {"type": "object"}}},
        handler=finalize, phase="evaluate"))
