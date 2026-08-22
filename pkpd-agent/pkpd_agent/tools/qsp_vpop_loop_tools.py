"""Model-agnostic loop tools for the virtual-population GENERATION task.

The agent chooses which disease-driver parameters to vary and over what bounds,
samples a candidate population, simulates each to its untreated baseline, and matches
the resulting baseline severity distribution to the clinical target. The driver
catalog and target come off the QSPTaskConfig.

  * ``vpop_inspect`` (observe) - the disease-driver parameters, the target, the band.
  * ``vpop_sample``  (act)     - set bounds, sample + simulate, return the RAW distribution.
  * ``vpop_select``  (act)     - reweight a wide pool to the target moments (numerical).
  * ``vpop_finalize``(evaluate)- commit the sampling design to be scored.
"""

from __future__ import annotations

from ..engines.simbiology import SimBiologyEngine
from ..engines import qsp_tasks
from ..engines.qsp_config import QSPTaskConfig
from .registry import Tool, ToolRegistry, ToolResult


def register_qsp_vpop_loop_tools(registry: ToolRegistry, config, ctx: dict) -> None:
    """ctx: {cfg, sb, n_samples, baseline_day, seed, target?, enable_select?, n_pool?}."""
    cfg: QSPTaskConfig = ctx["cfg"]
    sb: SimBiologyEngine = ctx["sb"]
    n_samples: int = int(ctx.get("n_samples") or 60)
    baseline_day: float = float(ctx.get("baseline_day")
                                or cfg.timeline.get("baseline_day", 200.0))
    seed: int = int(ctx.get("seed") or 1)
    target: dict = ctx.get("target") or cfg.vpop_target
    sev = cfg.severity_readout
    base_col = (cfg.run_columns.get("severity") or {}).get("baseline", "")
    enable_select: bool = bool(ctx.get("enable_select", True))
    n_pool: int = int(ctx.get("n_pool") or 80)

    def _catalog() -> list[dict]:
        return [{"name": name, "meaning": p.get("meaning"), "nominal": p.get("nominal"),
                 "observed_span": p.get("span")} for name, p in cfg.vpop_drivers.items()]

    def inspect(args: dict, session) -> ToolResult:
        return ToolResult.success(
            f"{cfg.name} virtual-population generation: sample disease-driver "
            f"parameters so the untreated baseline {sev} distribution matches a "
            "clinical target.",
            objective=f"build a virtual population whose baseline {sev} matches mean "
                      f"{target['mean']}, sd {target['sd']}, mostly inside the active "
                      f"band {target['band']}",
            disease_driver_parameters=_catalog(),
            clinical_target=target,
            n_candidates_per_sample=n_samples,
            edit_spec_help={"bounds": "{param: [lo, hi, scale]} - the sampling range "
                            "per parameter; scale is 'lin' or 'log'. Vary enough "
                            "drivers, over wide enough bounds, for phenotypic diversity "
                            "- but not so wide that most patients leave the active band."})

    def sample(args: dict, session) -> ToolResult:
        bounds = args.get("bounds") or {}
        if not bounds:
            return ToolResult.error("no bounds given - pass {bounds:{param:[lo,hi,scale]}}. "
                                    "Call vpop_inspect for the parameter names and spans.")
        spec = qsp_tasks.build_sample_spec(bounds)
        r = sb.sample_vpop(spec, n_samples=n_samples, baseline_day=baseline_day, seed=seed)
        das = (r.get("columns") or {}).get(base_col, [])
        score = qsp_tasks.score_vpop(das, target)

        hist = session.get("vpop_history") or []
        hist.append({"bounds": bounds, "score": score})
        session.put("vpop_history", hist)
        d = score.get("distribution_distance")
        best = session.get("vpop_best_dist")
        if d is not None and (best is None or d < best):
            session.put("vpop_best_dist", d)
            session.put("vpop_best_bounds", bounds)

        return ToolResult.success(
            f"sampled {score.get('n')} candidates: yield {score.get('yield_pct')}% in "
            f"band {target['band']}; accepted {sev} mean {score.get('accepted_mean')} "
            f"sd {score.get('accepted_sd')} (target {target['mean']}/{target['sd']}); "
            f"distance {d} (best {session.get('vpop_best_dist')})",
            bounds=bounds, **score,
            best_distance_so_far=session.get("vpop_best_dist"), iteration=len(hist))

    def select(args: dict, session) -> ToolResult:
        bounds = args.get("bounds") or {}
        if not bounds:
            return ToolResult.error("no bounds given - pass {bounds:{param:[lo,hi,scale]}}. "
                                    "Sample a pool WIDE enough to span the target range.")
        pool_n = int(args.get("n_pool") or n_pool)
        spec = qsp_tasks.build_sample_spec(bounds)
        r = sb.sample_vpop(spec, n_samples=pool_n, baseline_day=baseline_day, seed=seed)
        das = (r.get("columns") or {}).get(base_col, [])
        sel = qsp_tasks.select_to_moments(das, target)
        if not sel.get("ok"):
            return ToolResult.error(
                f"selection failed: {sel.get('reason')} (pool {sel.get('n_pool')}, "
                f"in-band {sel.get('n_inband')}). Widen the bounds so the pool spans "
                f"the active band {target['band']} and target mean {target['mean']}.", **sel)

        hist = session.get("vpop_history") or []
        hist.append({"bounds": bounds, "score": sel, "via": "numeric"})
        session.put("vpop_history", hist)
        d = sel.get("distribution_distance")
        best = session.get("vpop_best_dist")
        if d is not None and (best is None or d < best):
            session.put("vpop_best_dist", d)
            session.put("vpop_best_bounds", bounds)

        return ToolResult.success(
            f"pool {sel['n_pool']}, {sel['n_inband']} in band -> reweighted to target: "
            f"{sev} mean {sel['weighted_mean']} sd {sel['weighted_sd']} (target "
            f"{sel['target_mean']}/{sel['target_sd']}); distance {d}; effective sample "
            f"size {sel['effective_sample_size']} ({int(sel['ess_fraction']*100)}% of "
            f"in-band - low means the pool barely covers the target)",
            bounds=bounds, **sel, best_distance_so_far=session.get("vpop_best_dist"))

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
        mean = s.get("weighted_mean", s.get("accepted_mean"))
        sd = s.get("weighted_sd", s.get("accepted_sd"))
        extra = (f"ESS {s.get('effective_sample_size')}" if match.get("via") == "numeric"
                 else f"yield {s.get('yield_pct')}%")
        return ToolResult.success(
            f"committed vpop design ({match.get('via', 'raw')}): {sev} {mean}/{sd}, "
            f"distance {s.get('distribution_distance')}, {extra}", bounds=bounds, score=s)

    # -- multi-anchor selection (the paper's method: match SEVERAL clinical anchors) -- #
    anchors_cfg = cfg.vpop_anchors or {}

    def select_multi(args, s):
        bounds = args.get("bounds") or {}
        if not bounds:
            return ToolResult.error("no bounds given - pass {bounds:{param:[lo,hi,scale]}}.")
        arms = anchors_cfg.get("arms") or {}
        rate_targets = anchors_cfg.get("rate_targets") or {}
        spec = qsp_tasks.build_sample_spec(bounds)
        arms_spec = ";;".join(f"{lab}:{dose}" for lab, dose in arms.items())
        pool_n = int(args.get("n_pool") or n_pool)
        readout_day = cfg.timeline.get("first_line_readout_day", 284.0)
        r = sb.cohort_multi_arm(spec, arms_spec, baseline_day, readout_day, pool_n,
                                seed, states=cfg.readout_states or None)
        cols = r.get("columns") or {}
        sevs = cols.get("sev_base", [])
        candidates = []
        for i in range(len(sevs)):
            c = {"severity": sevs[i]}
            for lab in arms:
                col = cols.get(lab, [])
                if i < len(col):
                    c[lab] = col[i]
            candidates.append(c)
        anchors = [{"key": "severity", "mean": target["mean"], "sd": target.get("sd")}]
        anchors += [{"key": lab, "target": rate} for lab, rate in rate_targets.items()]
        sel = qsp_tasks.select_multi_anchor(candidates, anchors)
        if not sel.get("ok"):
            return ToolResult.error(f"multi-anchor selection failed: {sel.get('reason')}", **sel)
        s.put("vpop_history", (s.get("vpop_history") or [])
              + [{"bounds": bounds, "score": sel, "via": "multi_anchor"}])
        s.put("vpop_best_bounds", bounds)
        return ToolResult.success(
            f"multi-anchor selection over {sel['n']} candidates matched "
            f"{len(anchors)} anchors: {sel['anchors']}; total error {sel['total_error']}, "
            f"ESS {sel['effective_sample_size']} ({int(sel['ess_fraction']*100)}%).",
            bounds=bounds, **sel)

    registry.register(Tool(
        name="vpop_inspect",
        description=("OBSERVE the virtual-population task: the disease-driver parameters "
                     "(each with its nominal value and observed span), the clinical "
                     "target baseline distribution, and the active-disease band. Call first."),
        input_schema={"type": "object", "properties": {}},
        handler=inspect, phase="observe"))
    registry.register(Tool(
        name="vpop_sample",
        description=("ACT: sample a candidate population at the given bounds and report "
                     "the RAW baseline distribution (yield, mean/sd, distance). Use this "
                     "to PROBE. For the actual population, prefer vpop_select."),
        input_schema={"type": "object", "properties": {
            "bounds": {"type": "object", "description": "{param: [lo, hi, scale]}"}}},
        handler=sample, phase="act"))
    if enable_select:
        registry.register(Tool(
            name="vpop_select",
            description=("ACT: build the population NUMERICALLY. You choose which drivers "
                         "to vary and WIDE bounds that span the target; a prevalence-"
                         "weighting routine (the standard QSP method) reweights the pool "
                         "to the target moments and returns the reweighted mean/sd, the "
                         "distance, and the EFFECTIVE SAMPLE SIZE. A low ESS means the "
                         "pool barely covers the target (widen the bounds)."),
            input_schema={"type": "object", "properties": {
                "bounds": {"type": "object"}, "n_pool": {"type": "number"}}},
            handler=select, phase="act"))
    if anchors_cfg.get("arms"):
        registry.register(Tool(
            name="vpop_select_multi",
            description=("ACT: build the population to match SEVERAL clinical anchors at "
                         "once (the paper's method) - the baseline severity distribution "
                         "AND each therapy arm's response rate. Samples a cohort (each "
                         "candidate simulated under every arm in MATLAB), then optimizes "
                         "selection weights so the weighted population matches all "
                         "anchors. Returns the achieved value per anchor + the effective "
                         "sample size. Choose drivers + WIDE bounds spanning the targets."),
            input_schema={"type": "object", "properties": {
                "bounds": {"type": "object"}, "n_pool": {"type": "number"}}},
            handler=select_multi, phase="act"))
    registry.register(Tool(
        name="vpop_finalize",
        description=("COMMIT the sampling design you recommend (already run). Scores the "
                     "generated population against the clinical target. Call once."),
        input_schema={"type": "object", "properties": {"bounds": {"type": "object"}}},
        handler=finalize, phase="evaluate"))
