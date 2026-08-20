"""LLM-loop tools for the readout-mapping task: what drives the DAS28-CRP endpoint?

The agent proposes which model nodes (cells/mediators) the clinical readout DAS28-CRP is a
function of - reconstructing the mechanism->endpoint bridge. Scored node-recovery vs the
species the model's readout rule actually depends on.

  * ``readout_inspect``  (observe) - the readout, what it clinically measures, the objective.
  * ``readout_propose``  (act)     - submit the nodes you think drive DAS28-CRP.
  * ``readout_finalize`` (evaluate)- score vs the model's real readout dependencies.
"""

from __future__ import annotations

from ..engines import ra_readout as RO
from .registry import Tool, ToolRegistry, ToolResult


def register_ra_readout_loop_tools(registry: ToolRegistry, config, ctx: dict) -> None:
    """ctx: {drivers: list[str]} - the species the readout depends on (answer key)."""
    drivers: list = ctx["drivers"]

    def inspect(args: dict, session) -> ToolResult:
        return ToolResult.success(
            "RA readout-mapping task: the model computes DAS28-CRP (disease severity) from "
            "its physiological state. Propose which model nodes (cells / mediators) the "
            "DAS28-CRP readout is built from - the mechanism-to-endpoint mapping.",
            readout="DAS28-CRP - clinically a composite of tender + swollen joint counts, "
                    "CRP (or ESR), and patient global assessment",
            objective="list the model nodes (cell densities and/or mediators) that the "
                      "DAS28-CRP value is a direct function of",
            hint="think what the score physically reflects: the load of infiltrating and "
                 "structural cells in the synovium (swelling) and the systemic acute-phase "
                 "signal (CRP, driven by a particular cytokine). It is a readout formula, "
                 "not the whole network.")

    def propose(args: dict, session) -> ToolResult:
        names = [str(n) for n in (args.get("nodes") or [])]
        acc = list(session.get("readout_nodes") or [])
        acc.extend(names)
        session.put("readout_nodes", acc)
        return ToolResult.success(f"recorded {len(names)} node(s); {len(set(acc))} distinct.",
                                  n_distinct=len(set(acc)))

    def finalize(args: dict, session) -> ToolResult:
        nodes = session.get("readout_nodes") or []
        if not nodes:
            return ToolResult.error("propose the readout drivers before finalizing.")
        sc = RO.score_readout(nodes, drivers)
        session.put("readout_final", sc)
        return ToolResult.success(
            f"scored: {sc['hit']}/{sc['n_truth']} real readout drivers recovered "
            f"(F1 {sc['f1']}, P {sc['precision']} / R {sc['recall']}).", **sc)

    registry.register(Tool(
        name="readout_inspect",
        description=("OBSERVE the readout-mapping task: what DAS28-CRP measures and the "
                     "objective. No answer key. Call first."),
        input_schema={"type": "object", "properties": {}},
        handler=inspect, phase="observe"))

    registry.register(Tool(
        name="readout_propose",
        description=("ACT: submit the model nodes you think DAS28-CRP is computed from "
                     "(a list of names)."),
        input_schema={"type": "object", "properties": {
            "nodes": {"type": "array", "items": {"type": "string"}}}},
        handler=propose, phase="act"))

    registry.register(Tool(
        name="readout_finalize",
        description=("COMMIT the proposed readout drivers; scored (precision/recall/F1) "
                     "vs the species the model's readout rule depends on. Terminal."),
        input_schema={"type": "object", "properties": {}},
        handler=finalize, phase="evaluate"))
