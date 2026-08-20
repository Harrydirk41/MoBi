"""LLM-loop tools for the Stage-2 task: choose the model's SCOPE (cast of nodes).

The agent is given the disease and the modeling goal and must propose the cell types and
soluble mediators to include - and, implicitly, which to leave out. Scored against the
model's real 26-node cast at finalize. Precision is the real test: over-including the whole
RA textbook is penalised, so this measures scope discipline, not just biology recall.

  * ``scope_inspect``  (observe) - the disease, the modeling goal, the node format. No key.
  * ``scope_propose``  (act)     - submit cell/mediator names; accumulate across calls.
  * ``scope_finalize`` (evaluate)- score the cast (precision/recall/F1) vs the model.
"""

from __future__ import annotations

from ..engines import ra_scope as S
from .registry import Tool, ToolRegistry, ToolResult


def register_ra_scope_loop_tools(registry: ToolRegistry, config, ctx: dict) -> None:
    """ctx: {} - the answer key (the model's cast) lives in ra_scope, never shown."""

    def inspect(args: dict, session) -> ToolResult:
        return ToolResult.success(
            "RA scope-selection task: choose the CAST of a rheumatoid-arthritis QSP model "
            "- which cell types and soluble mediators to include. This is the scoping "
            "decision made before any wiring; you are scored against the model's real cast.",
            disease="rheumatoid arthritis (inflamed synovial joint)",
            modeling_goal="simulate late-phase (Phase 2/3) trials of DMARDs and biologics "
                          "and reproduce ACR20/50/70 and DAS28-CRP endpoints",
            instruction="propose the cell types AND the soluble mediators (cytokines, "
                        "chemokines, growth factors, adhesion molecules, autoantibodies) "
                        "that belong in the model",
            scoring="precision/recall/F1 vs the model's actual cast. Recall is easy; the "
                    "test is PRECISION - a good QSP model is parsimonious, so including "
                    "every RA-associated mediator is penalised. Include what drives the "
                    "modeled biology and the endpoints, and leave out the rest.",
            node_format="a flat list of names - each entry is one cell type or one "
                        "soluble mediator, given by its standard name")

    def propose(args: dict, session) -> ToolResult:
        names = args.get("nodes") or []
        acc = list(session.get("scope_nodes") or [])
        acc.extend(str(n) for n in names)
        session.put("scope_nodes", acc)
        return ToolResult.success(
            f"accepted {len(names)} name(s); {len(set(map(str.lower, acc)))} distinct so far.",
            n_distinct=len(set(map(str.lower, acc))))

    def finalize(args: dict, session) -> ToolResult:
        nodes = session.get("scope_nodes") or []
        if not nodes:
            return ToolResult.error("propose the cast before finalizing.")
        sc = S.score_scope(nodes)
        session.put("scope_final", vars(sc))
        return ToolResult.success(
            f"committed a cast of {sc.hit + sc.extra} nodes. F1 {sc.f1} "
            f"(P {sc.precision} / R {sc.recall}). Recovered {sc.hit}/"
            f"{sc.hit + sc.missed} real nodes; {sc.extra} not in the model "
            f"({len(sc.extra_known_mediators)} of them real RA mediators it excludes).",
            **vars(sc))

    registry.register(Tool(
        name="scope_inspect",
        description=("OBSERVE the scope-selection task: the disease, the modeling goal, "
                     "and the node format. No answer key. Call first."),
        input_schema={"type": "object", "properties": {}},
        handler=inspect, phase="observe"))

    registry.register(Tool(
        name="scope_propose",
        description=("ACT: submit cell/mediator names for the model's cast (a list of "
                     "strings). Accumulates across calls. Be deliberate about scope - "
                     "precision is scored, so do not dump every RA-associated molecule."),
        input_schema={"type": "object", "properties": {
            "nodes": {"type": "array", "items": {"type": "string"}}}},
        handler=propose, phase="act"))

    registry.register(Tool(
        name="scope_finalize",
        description=("COMMIT the proposed cast; scored (precision/recall/F1) against the "
                     "model's real 26-node cast. Terminal - call once."),
        input_schema={"type": "object", "properties": {}},
        handler=finalize, phase="evaluate"))
