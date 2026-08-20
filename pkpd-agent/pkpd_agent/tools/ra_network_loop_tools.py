"""LLM-loop tools for the Stage-1 task: RECONSTRUCT the RA disease network.

The agent is given only the cast (cells + cytokines) and must propose the regulatory
network - directed, signed edges "source (up/down)-regulates target". No simulation:
this is pure biological reasoning, scored against the real model's wiring (the answer
key) only at finalize, so the agent cannot hill-climb against the key.

  * ``network_inspect``  (observe) - the cast, the edge format, the objective. No key.
  * ``network_propose``  (act)     - submit a batch of edges; get STRUCTURAL feedback
    (counts, nodes still with no edges) that never leaks the truth. Call repeatedly
    to build/refine the draft.
  * ``network_finalize`` (evaluate)- commit the accumulated edge set; it is scored
    (precision/recall/F1, sign-aware and topology-only) against the model. Terminal.
"""

from __future__ import annotations

from ..engines import ra_network as N
from .registry import Tool, ToolRegistry, ToolResult


def register_ra_network_loop_tools(registry: ToolRegistry, config, ctx: dict) -> None:
    """ctx: {truth: list[Edge]}  (parsed from the model - the answer key, never shown)."""
    truth: list = ctx["truth"]

    def inspect(args: dict, session) -> ToolResult:
        return ToolResult.success(
            "RA network-reconstruction task: propose the regulatory wiring of a "
            "published RA QSP model from immunology alone. No simulation - you are "
            "drawing the model's structure, scored against its real edges at finalize.",
            objective="propose every directed, signed regulatory edge among the cast: "
                      "which cell/cytokine UP- or DOWN-regulates the secretion / "
                      "proliferation / influx of which cytokine/cell",
            cells=N.CELLS,
            cytokines=N.CYTOKINES,
            edge_format={"source": "a node from the cast", "target": "a node from the cast",
                         "sign": "+1 promote / -1 inhibit"},
            guidance=[
                "Think mechanistically: which cells SECRETE each cytokine (a source->"
                "cytokine edge), and which cytokines DRIVE each cell's proliferation / "
                "influx / survival (a cytokine->cell edge).",
                "Include the well-established RA axes: TNF-a and IL-6 from macrophages "
                "and FLS; Th17/IL-17; IFN-g from Th1/CTL; regulatory (Treg, anti-"
                "inflammatory) edges as negative signs.",
                "Aim for coverage AND precision - guessing every pair is penalised "
                "(precision), missing real edges is penalised (recall)."],
            hint="propose in batches with network_propose; refine using the structural "
                 "feedback; then network_finalize once.")

    def propose(args: dict, session) -> ToolResult:
        items = args.get("edges") or []
        new = N.edges_from_proposal(items)
        acc = {e.signed(): e for e in (session.get("net_edges") or [])}
        added = 0
        for e in new:
            if e.signed() not in acc:
                acc[e.signed()] = e
                added += 1
        edges = list(acc.values())
        session.put("net_edges", edges)
        # structural feedback only - never references the truth
        covered = {e.source for e in edges} | {e.target for e in edges}
        no_out = [n for n in N.NODES if n not in {e.source for e in edges}]
        no_in = [n for n in N.NODES if n not in {e.target for e in edges}]
        return ToolResult.success(
            f"accepted {added} new edge(s); {len(edges)} total covering "
            f"{len(covered)}/{len(N.NODES)} nodes.",
            total_edges=len(edges),
            nodes_with_no_outgoing=no_out,
            nodes_with_no_incoming=no_in)

    def finalize(args: dict, session) -> ToolResult:
        edges = session.get("net_edges") or []
        if not edges:
            return ToolResult.error("propose edges before finalizing.")
        sa = N.score_network(edges, truth, sign_aware=True)
        topo = N.score_network(edges, truth, sign_aware=False)
        session.put("net_final", {"sign_aware": sa, "topology": topo,
                                  "n_edges": len(edges)})
        return ToolResult.success(
            f"committed {len(edges)} edges. Topology F1 {topo['f1']} "
            f"(P {topo['precision']} / R {topo['recall']}); sign-aware F1 {sa['f1']}. "
            f"Recovered {topo['hit']}/{topo['n_truth']} real interactions, "
            f"{topo['extra']} not in the model.",
            sign_aware=sa, topology=topo)

    registry.register(Tool(
        name="network_inspect",
        description=("OBSERVE the network-reconstruction task: the cast of cells and "
                     "cytokines, the signed-edge format, and the objective. No answer "
                     "key is given. Call first."),
        input_schema={"type": "object", "properties": {}},
        handler=inspect, phase="observe"))

    registry.register(Tool(
        name="network_propose",
        description=("ACT: submit a batch of proposed regulatory edges. Each edge is "
                     "{source, target, sign(+1/-1)}. Returns structural feedback "
                     "(coverage, nodes still unconnected) - never the truth. Call "
                     "repeatedly to build and refine the network draft."),
        input_schema={"type": "object", "properties": {
            "edges": {"type": "array", "items": {"type": "object", "properties": {
                "source": {"type": "string"}, "target": {"type": "string"},
                "sign": {"type": "number", "description": "+1 promote, -1 inhibit"}}}}}},
        handler=propose, phase="act"))

    registry.register(Tool(
        name="network_finalize",
        description=("COMMIT the accumulated network; it is scored against the real "
                     "model (precision/recall/F1, sign-aware and topology-only). "
                     "Terminal - call once when the draft is complete."),
        input_schema={"type": "object", "properties": {}},
        handler=finalize, phase="evaluate"))
