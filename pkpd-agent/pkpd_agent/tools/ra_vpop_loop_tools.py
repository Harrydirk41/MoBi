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
    enable_select: bool = bool(ctx.get("enable_select", True))
    n_pool: int = int(ctx.get("n_pool") or 80)

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

    # -- act via NUMERICAL selection (agent designs, routine selects) --- #
    def select(args: dict, session) -> ToolResult:
        bounds = args.get("bounds") or {}
        if not bounds:
            return ToolResult.error(
                "no bounds given - pass {bounds:{param:[lo,hi,scale]}}. Sample a pool "
                "WIDE enough to span the target range; the routine does the selection.")
        pool_n = int(args.get("n_pool") or n_pool)
        spec = osp_ra_trial.build_sample_spec(bounds)
        r = sb.sample_vpop(spec, n_samples=pool_n, baseline_day=baseline_day, seed=seed)
        das = (r.get("columns") or {}).get("DAS28_base", [])
        sel = osp_ra_trial.select_to_moments(das, target)
        if not sel.get("ok"):
            return ToolResult.error(
                f"selection failed: {sel.get('reason')} (pool {sel.get('n_pool')}, "
                f"in-band {sel.get('n_inband')}). Widen the bounds so the pool spans "
                f"the active band {target['band']} and the target mean {target['mean']}.",
                **sel)

        hist = session.get("vpop_history") or []
        hist.append({"bounds": bounds, "score": sel, "via": "numeric"})
        session.put("vpop_history", hist)
        best = session.get("vpop_best_dist")
        d = sel.get("distribution_distance")
        if d is not None and (best is None or d < best):
            session.put("vpop_best_dist", d)
            session.put("vpop_best_bounds", bounds)

        return ToolResult.success(
            f"pool {sel['n_pool']}, {sel['n_inband']} in band -> reweighted to target: "
            f"DAS28 mean {sel['weighted_mean']} sd {sel['weighted_sd']} "
            f"(target {sel['target_mean']}/{sel['target_sd']}); distance "
            f"{sel['distribution_distance']}; effective sample size "
            f"{sel['effective_sample_size']} ({int(sel['ess_fraction']*100)}% of "
            f"in-band - low means the pool barely covers the target)",
            bounds=bounds, **sel,
            best_distance_so_far=session.get("vpop_best_dist"))

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
        # numeric selection reports weighted_mean/sd + ESS; raw sampling reports
        # accepted_mean/sd + yield
        mean = s.get("weighted_mean", s.get("accepted_mean"))
        sd = s.get("weighted_sd", s.get("accepted_sd"))
        extra = (f"ESS {s.get('effective_sample_size')}" if match.get("via") == "numeric"
                 else f"yield {s.get('yield_pct')}%")
        return ToolResult.success(
            f"committed vpop design ({match.get('via', 'raw')}): DAS28 {mean}/{sd}, "
            f"distance {s.get('distribution_distance')}, {extra}",
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
            "ACT: sample a candidate population at the given bounds and report the "
            "RAW baseline DAS28-CRP distribution (yield in the active band, mean/sd, "
            "distance). Use this to PROBE - see where a set of bounds lands the raw "
            "distribution - not to hand-tune the population to target. For the actual "
            "population, prefer vpop_select, which reweights a pool to the target."),
        input_schema={"type": "object", "properties": {
            "bounds": {"type": "object",
                       "description": "{param: [lo, hi, scale]}, scale 'lin' or 'log'"}}},
        handler=sample, phase="act"))

    if enable_select:
        registry.register(Tool(
            name="vpop_select",
            description=(
                "ACT: build the population NUMERICALLY. You set up the design - which "
                "disease drivers to vary and WIDE bounds that span the target range - "
                "and a prevalence-weighting routine (the standard QSP method, what the "
                "paper's genetic algorithm implements) reweights the sampled pool so "
                "its baseline DAS28-CRP matches the target moments. Returns the "
                "reweighted mean/sd, the distance, and the EFFECTIVE SAMPLE SIZE. Your "
                "job is choosing drivers + bounds wide enough to cover the target, and "
                "judging the result: a low effective sample size means the pool barely "
                "covers the target (widen the bounds); a healthy one means a diverse, "
                "on-target population."),
            input_schema={"type": "object", "properties": {
                "bounds": {"type": "object",
                           "description": "{param: [lo, hi, scale]}; sample WIDE"},
                "n_pool": {"type": "number",
                           "description": "candidate pool size (default 80); each is a "
                                          "baseline sim, so keep modest"}}},
            handler=select, phase="act"))

    registry.register(Tool(
        name="vpop_finalize",
        description=(
            "COMMIT the sampling design you recommend (already run with vpop_sample). "
            "Scores the generated population against the clinical target distribution. "
            "Call once before finishing."),
        input_schema={"type": "object", "properties": {
            "bounds": {"type": "object"}}},
        handler=finalize, phase="evaluate"))
