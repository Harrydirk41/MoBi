"""Model-agnostic loop tools for the drug-DESIGN task: design a new pathway inhibitor.

The agent designs a drug from scratch - it chooses which disease-driver PATHWAY to
inhibit and how hard (efficacy) - and the harness edits the model structurally (adds
the drug as a time-gated suppression of that pathway's driver parameter), simulates,
and reports the clinical response. The targetable pathways come off the QSPTaskConfig.

  * ``design_inspect``  (observe) - the targetable pathways, the readout, the objective.
  * ``design_try``      (act)     - build a drug {target, efficacy}, edit + simulate.
  * ``design_finalize`` (evaluate)- commit the drug design to be scored.

Each design_try reloads the model (clean slate) so designs never contaminate.
"""

from __future__ import annotations

from ..engines.simbiology import SimBiologyEngine
from ..engines import qsp_tasks
from ..engines.qsp_config import QSPTaskConfig
from .registry import Tool, ToolRegistry, ToolResult


def register_qsp_design_loop_tools(registry: ToolRegistry, config, ctx: dict) -> None:
    """ctx: {cfg, sb, sbproj, vpop, background?, start_day, readout_day, limit}."""
    cfg: QSPTaskConfig = ctx["cfg"]
    sb: SimBiologyEngine = ctx["sb"]
    sbproj: str = ctx["sbproj"]
    vpop: str = ctx["vpop"]
    targets = cfg.design_targets
    background: str = ctx.get("background", cfg.design_background)
    start_day: float = float(ctx.get("start_day")
                             or cfg.timeline.get("baseline_day", 200.0))
    readout_day: float = float(ctx.get("readout_day")
                               or cfg.timeline.get("first_line_readout_day", 284.0))
    limit: int = int(ctx.get("limit") or 40)

    def _targets() -> list[dict]:
        return [{"target_param": name, "pathway": d.get("pathway"),
                 "real_analogue": d.get("analogue"), "note": d.get("note")}
                for name, d in targets.items()]

    def inspect(args: dict, session) -> ToolResult:
        return ToolResult.success(
            f"{cfg.name} drug-design task: design a new drug by choosing a disease-"
            "driver pathway to inhibit and how strongly; the model is edited to add "
            "the drug and the clinical response is simulated.",
            objective=f"design the drug (pathway + efficacy) that gives the best "
                      f"{cfg.readout_desc} on a background of {background or 'none'}",
            background_therapy=background or "none (monotherapy)",
            readout=f"{cfg.readout_desc} at day {readout_day:g}",
            baseline="run efficacy 0 for the background-alone baseline to measure the "
                     "drug's added benefit",
            targetable_pathways=_targets(),
            edit_spec_help={
                "target": "the pathway's driver parameter, e.g. "
                          f"'{next(iter(targets), 'F_X')}'",
                "efficacy": "fractional pathway suppression in [0,1] (0 = no drug, "
                            "1 = full blockade)"})

    def _run(target: str, efficacy: float) -> dict:
        sb.load_project(sbproj)
        if efficacy > 0 and target:
            sb.add_drug(target, efficacy, start_day)
        r = sb.run_vpop(vpop, dose=background, stop_time=400.0, baseline_day=start_day,
                        readout_day=readout_day, limit=limit)
        fl = cfg.summarize_run(r)["first_line"]
        if efficacy > 0 and target:
            sb.load_project(sbproj)
        return fl

    def try_design(args: dict, session) -> ToolResult:
        target = args.get("target")
        efficacy = args.get("efficacy")
        if efficacy is None:
            return ToolResult.error("give efficacy (0..1) and a target pathway.")
        if efficacy > 0 and not target:
            return ToolResult.error("give a target pathway for a non-zero efficacy. "
                                    "Call design_inspect for the targetable pathways.")
        if target and target not in targets:
            return ToolResult.error(f"unknown target '{target}'. Options: "
                                    f"{', '.join(targets)}")
        fl = _run(target or "", float(efficacy))
        primary = next(iter(cfg.run_columns.get("first_line", {"ACR20": 1})), "ACR20")
        score = fl.get(primary)

        hist = session.get("design_history") or []
        hist.append({"target": target, "efficacy": efficacy, "response": fl})
        session.put("design_history", hist)
        best = session.get("design_best_score")
        if score is not None and (best is None or score > best):
            session.put("design_best_score", score)
            session.put("design_best", {"target": target, "efficacy": efficacy})

        label = (f"{targets.get(target, {}).get('pathway', target)} @ efficacy "
                 f"{efficacy:g}" if target and efficacy > 0 else "background-alone baseline")
        return ToolResult.success(
            f"{label}: {fl} (best {primary} so far {session.get('design_best_score')})",
            target=target, efficacy=efficacy, response=fl,
            best_score_so_far=session.get("design_best_score"), iteration=len(hist))

    def finalize(args: dict, session) -> ToolResult:
        target = args.get("target") or (session.get("design_best") or {}).get("target")
        efficacy = args.get("efficacy")
        if efficacy is None:
            efficacy = (session.get("design_best") or {}).get("efficacy")
        hist = session.get("design_history") or []
        match = next((h for h in reversed(hist)
                      if h.get("target") == target and h.get("efficacy") == efficacy), None)
        if match is None:
            return ToolResult.error("that exact design was not run - "
                                    "design_try it first, then finalize.")
        session.put("design_final", match)
        d = targets.get(target, {})
        return ToolResult.success(
            f"committed drug design: anti-{d.get('pathway', target)} (analogue "
            f"{d.get('analogue', 'n/a')}) at efficacy {efficacy:g}; response "
            f"{match['response']}", target=target, efficacy=efficacy,
            response=match["response"])

    registry.register(Tool(
        name="design_inspect",
        description=("OBSERVE the drug-design task: the pathways you can target (each "
                     "with mechanism and real-world analogue), the readout, and the "
                     "objective. Call this first."),
        input_schema={"type": "object", "properties": {}},
        handler=inspect, phase="observe"))
    registry.register(Tool(
        name="design_try",
        description=("ACT: design a drug by choosing a target pathway and efficacy; the "
                     "model is edited to add it and simulated, returning the response. "
                     "Run efficacy 0 once for the background-alone baseline, then screen "
                     "pathways, then tune. Each call is an independent design."),
        input_schema={"type": "object", "properties": {
            "target": {"type": "string"}, "efficacy": {"type": "number"}}},
        handler=try_design, phase="act"))
    registry.register(Tool(
        name="design_finalize",
        description=("COMMIT the drug design you recommend (target + efficacy, already "
                     "run with design_try). Scores the design and names its real-world "
                     "analogue. Call once."),
        input_schema={"type": "object", "properties": {
            "target": {"type": "string"}, "efficacy": {"type": "number"}}},
        handler=finalize, phase="evaluate"))
