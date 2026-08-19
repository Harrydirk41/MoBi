"""LLM-loop tools for the RA Stage-2 task: DESIGN a new anti-cytokine biologic.

The agent designs a drug from scratch - it chooses which cytokine PATHWAY to inhibit
and how hard (efficacy) - and the harness edits the model structurally (adds the drug
as a time-gated suppression of that pathway's disease-driver parameter), simulates,
and reports the clinical response. The agent then screens pathways and tunes the drug,
discovering which target is the best RA drug (and should recover that IL-6 / TNF
dominate, matching real biology).

  * ``design_inspect``  (observe) - the cytokine pathways available to target (with
    mechanism and real-world analogue), the readout, and the objective.
  * ``design_try``      (act)     - build a drug {target, efficacy}, edit the model,
    simulate, and return the ACR response (vs the no-drug / MTX-alone baseline).
  * ``design_finalize`` (evaluate)- commit the drug design to be scored.

Each design_try reloads the model (clean slate), adds the designed drug, and runs -
so designs never contaminate each other.
"""

from __future__ import annotations

from typing import Any

from ..engines.simbiology import SimBiologyEngine
from ..engines import osp_ra_trial
from .registry import Tool, ToolRegistry, ToolResult


def register_ra_design_loop_tools(registry: ToolRegistry, config, ctx: dict) -> None:
    """ctx: {sb, sbproj, vpop, background, start_day, readout_day, limit}."""
    sb: SimBiologyEngine = ctx["sb"]
    sbproj: str = ctx["sbproj"]
    vpop: str = ctx["vpop"]
    background: str = ctx.get("background", "MTX_15mg_Q1W_SC_t200")
    start_day: float = float(ctx.get("start_day") or 200.0)
    readout_day: float = float(ctx.get("readout_day") or 284.0)
    limit: int = int(ctx.get("limit") or 40)

    def _targets() -> list[dict]:
        out = []
        for name, d in osp_ra_trial.DESIGN_TARGETS.items():
            out.append({"target_param": name, "pathway": d["pathway"],
                        "real_analogue": d["analogue"], "note": d["note"]})
        return out

    # -- observe -------------------------------------------------------- #
    def inspect(args: dict, session) -> ToolResult:
        return ToolResult.success(
            "RA drug-design task: design a new biologic by choosing a cytokine "
            "pathway to inhibit and how strongly; the model is edited to add the "
            "drug and the clinical response is simulated.",
            objective="design the drug (pathway + efficacy) that gives the best ACR "
                      "response on a background of methotrexate",
            background_therapy=background or "none (monotherapy)",
            readout=f"first-line ACR20/50/70 + remission at day {readout_day:g}",
            baseline="run efficacy 0 (or any target at efficacy 0) for the MTX-alone "
                     "baseline to measure the drug's added benefit",
            targetable_pathways=_targets(),
            edit_spec_help={
                "target": "the pathway's driver parameter, e.g. 'F_IL6' (IL-6), "
                          "'F_TNFa' (TNF), 'F_IL17' (IL-17)",
                "efficacy": "fractional pathway suppression in [0,1] - the drug's "
                            "effect size (0 = no drug, 1 = full blockade)"},
        )

    def _run(target: str, efficacy: float) -> dict:
        sb.load_project(sbproj)                      # clean slate each design
        if efficacy > 0 and target:
            sb.add_drug(target, efficacy, start_day)
        r = sb.run_vpop(vpop, dose=background, stop_time=400.0,
                        baseline_day=start_day, readout_day=readout_day, limit=limit)
        return osp_ra_trial.summarize_run(r)["first_line"]

    # -- act ------------------------------------------------------------ #
    def try_design(args: dict, session) -> ToolResult:
        target = args.get("target")
        efficacy = args.get("efficacy")
        if efficacy is None:
            return ToolResult.error("give efficacy (0..1) and a target pathway.")
        if efficacy > 0 and not target:
            return ToolResult.error(
                "give a target pathway (e.g. 'F_IL6') for a non-zero efficacy. "
                "Call design_inspect for the targetable pathways.")
        if target and target not in osp_ra_trial.DESIGN_TARGETS:
            return ToolResult.error(
                f"unknown target '{target}'. Options: "
                f"{', '.join(osp_ra_trial.DESIGN_TARGETS)}")
        fl = _run(target or "", float(efficacy))
        acr20 = fl.get("ACR20")

        hist = session.get("design_history") or []
        hist.append({"target": target, "efficacy": efficacy, "response": fl})
        session.put("design_history", hist)
        best = session.get("design_best_acr20")
        if acr20 is not None and (best is None or acr20 > best):
            session.put("design_best_acr20", acr20)
            session.put("design_best", {"target": target, "efficacy": efficacy})

        label = (f"{osp_ra_trial.DESIGN_TARGETS.get(target, {}).get('pathway', target)} "
                 f"@ efficacy {efficacy:g}" if target and efficacy > 0
                 else "MTX-alone baseline")
        return ToolResult.success(
            f"{label}: ACR20 {acr20}% ACR50 {fl.get('ACR50')}% ACR70 "
            f"{fl.get('ACR70')}% remission {fl.get('remission')}% "
            f"(best ACR20 so far {session.get('design_best_acr20')})",
            target=target, efficacy=efficacy, response=fl,
            best_acr20_so_far=session.get("design_best_acr20"), iteration=len(hist))

    # -- evaluate ------------------------------------------------------- #
    def finalize(args: dict, session) -> ToolResult:
        target = args.get("target") or (session.get("design_best") or {}).get("target")
        efficacy = args.get("efficacy")
        if efficacy is None:
            efficacy = (session.get("design_best") or {}).get("efficacy")
        hist = session.get("design_history") or []
        match = next((h for h in reversed(hist)
                      if h.get("target") == target and h.get("efficacy") == efficacy), None)
        if match is None:
            return ToolResult.error(
                "that exact design was not run - design_try it first, then finalize.")
        session.put("design_final", match)
        d = osp_ra_trial.DESIGN_TARGETS.get(target, {})
        return ToolResult.success(
            f"committed drug design: anti-{d.get('pathway', target)} "
            f"(analogue {d.get('analogue', 'n/a')}) at efficacy {efficacy:g}; "
            f"ACR20 {match['response'].get('ACR20')}%",
            target=target, efficacy=efficacy, response=match["response"])

    registry.register(Tool(
        name="design_inspect",
        description=(
            "OBSERVE the drug-design task: the cytokine pathways you can target "
            "(each with its mechanism and real-world drug analogue), the ACR readout, "
            "and the objective. Call this first."),
        input_schema={"type": "object", "properties": {}},
        handler=inspect, phase="observe"))

    registry.register(Tool(
        name="design_try",
        description=(
            "ACT: design a drug by choosing a target pathway and efficacy; the model "
            "is edited to add it and simulated, returning the ACR response. Run "
            "efficacy 0 once for the MTX-alone baseline, then screen pathways at a "
            "meaningful efficacy to see which cytokine is the best target, then tune. "
            "Each call is an independent design (the model is reset)."),
        input_schema={"type": "object", "properties": {
            "target": {"type": "string", "description": "pathway driver, e.g. 'F_IL6'"},
            "efficacy": {"type": "number", "description": "suppression in [0,1]"}}},
        handler=try_design, phase="act"))

    registry.register(Tool(
        name="design_finalize",
        description=(
            "COMMIT the drug design you recommend (target + efficacy, already run "
            "with design_try). Scores the design and names its real-world analogue. "
            "Call once before finishing."),
        input_schema={"type": "object", "properties": {
            "target": {"type": "string"},
            "efficacy": {"type": "number"}}},
        handler=finalize, phase="evaluate"))
