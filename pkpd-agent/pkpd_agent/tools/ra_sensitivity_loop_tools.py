"""LLM-loop tools for the sensitivity task: pick the most influential parameters.

The agent is shown a pool of real model parameters (the GSA top-20 hidden among
distractors) and must select/rank the ones that most drive DAS28-CRP variance. Scored by
overlap with the paper's global sensitivity analysis, vs a blind-pick random baseline.

  * ``sens_inspect``  (observe) - the candidate pool, the readout, the objective. No key.
  * ``sens_rank``     (act)     - submit an ordered list of the most sensitive parameters.
  * ``sens_finalize`` (evaluate)- score overlap + rank correlation vs the GSA. Terminal.
"""

from __future__ import annotations

from ..engines import ra_sensitivity as SEN
from .registry import Tool, ToolRegistry, ToolResult


def register_ra_sensitivity_loop_tools(registry: ToolRegistry, config, ctx: dict) -> None:
    def inspect(args: dict, session) -> ToolResult:
        return ToolResult.success(
            "RA sensitivity task: from a pool of model parameters, pick the ones that most "
            "drive the DAS28-CRP disease-severity readout (global sensitivity). This tests "
            "which knobs matter - dynamical leverage, not biology recall.",
            readout="DAS28-CRP (overall disease severity)",
            objective="select and RANK (most sensitive first) about 20 parameters from the "
                      "pool that contribute most to variance in DAS28-CRP",
            hint="baseline growth/influx rates of the abundant cells and the fractional "
                 "disease-driver terms tend to have high leverage; downstream effect "
                 "strengths and drug-binding constants usually less. Rank order matters.",
            candidate_pool=SEN.pool())

    def rank(args: dict, session) -> ToolResult:
        picks = [str(p) for p in (args.get("ranked") or args.get("parameters") or [])]
        session.put("sens_ranked", picks)
        pool = set(SEN.pool())
        off = [p for p in picks if p not in pool]
        return ToolResult.success(
            f"recorded {len(picks)} ranked parameter(s)."
            + (f" {len(off)} not in the pool (ignored at scoring): {off[:5]}" if off else ""),
            n_ranked=len(picks))

    def finalize(args: dict, session) -> ToolResult:
        picks = session.get("sens_ranked") or []
        if not picks:
            return ToolResult.error("rank the parameters before finalizing.")
        sc = SEN.score_sensitivity(picks)
        session.put("sens_final", sc)
        return ToolResult.success(
            f"scored: {sc['hit']}/20 GSA parameters recovered (recall {sc['recall']}, "
            f"precision {sc['precision']}). Random blind-pick baseline recall "
            f"{sc['random_baseline_recall']}; beats random: {sc['beats_random']}. "
            f"Rank correlation on hits: {sc['spearman_on_hits']}.",
            **sc)

    registry.register(Tool(
        name="sens_inspect",
        description=("OBSERVE the sensitivity task: the candidate parameter pool, the "
                     "DAS28-CRP readout, and the objective. No answer key. Call first."),
        input_schema={"type": "object", "properties": {}},
        handler=inspect, phase="observe"))

    registry.register(Tool(
        name="sens_rank",
        description=("ACT: submit the parameters you judge most sensitive for DAS28-CRP as "
                     "an ORDERED list (most sensitive first), ~20 of them, drawn from the "
                     "pool."),
        input_schema={"type": "object", "properties": {
            "ranked": {"type": "array", "items": {"type": "string"}}}},
        handler=rank, phase="act"))

    registry.register(Tool(
        name="sens_finalize",
        description=("COMMIT your ranking; scored (overlap + rank correlation) against the "
                     "model's global sensitivity analysis. Terminal - call once."),
        input_schema={"type": "object", "properties": {}},
        handler=finalize, phase="evaluate"))
