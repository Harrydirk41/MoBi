"""Model-AGNOSTIC topology loop tools: same task as ra_network, driven by a QSPModel.

Instead of the hardcoded RA cast + answer key, these read everything off a QSPModel derived
from network.json + a spec. Point at any QSP model and the topology benchmark runs unchanged.
"""

from __future__ import annotations

from ..engines import qsp_core as N
from ..engines.qsp_model import QSPModel
from .registry import Tool, ToolRegistry, ToolResult


def register_qsp_topology_loop_tools(registry: ToolRegistry, config, ctx: dict) -> None:
    """ctx: {model: QSPModel}."""
    model: QSPModel = ctx["model"]

    def inspect(args: dict, session) -> ToolResult:
        return ToolResult.success(
            f"{model.spec.name}: propose the regulatory wiring of this QSP model from "
            "biology alone - directed, signed edges among the nodes. Scored vs the "
            "model's real edges at finalize.",
            objective="propose every directed, signed regulatory edge among the nodes",
            nodes=model.nodes,
            edge_format={"source": "a node", "target": "a node", "sign": "+1 / -1"})

    def propose(args: dict, session) -> ToolResult:
        new = model.edges_from_proposal(args.get("edges") or [])
        acc = {e.signed(): e for e in (session.get("qsp_edges") or [])}
        added = 0
        for e in new:
            if e.signed() not in acc:
                acc[e.signed()] = e
                added += 1
        session.put("qsp_edges", list(acc.values()))
        return ToolResult.success(f"accepted {added} new edge(s); {len(acc)} total.",
                                  total_edges=len(acc))

    def finalize(args: dict, session) -> ToolResult:
        edges = session.get("qsp_edges") or []
        if not edges:
            return ToolResult.error("propose edges before finalizing.")
        sa = N.score_network(edges, model.edges, sign_aware=True)
        topo = N.score_network(edges, model.edges, sign_aware=False)
        session.put("qsp_topo_final", {"sign_aware": sa, "topology": topo,
                                       "n_edges": len(edges)})
        return ToolResult.success(
            f"committed {len(edges)} edges. Topology F1 {topo['f1']} "
            f"(P {topo['precision']} / R {topo['recall']}); sign-aware F1 {sa['f1']}. "
            f"Recovered {topo['hit']}/{topo['n_truth']} real interactions.",
            sign_aware=sa, topology=topo)

    registry.register(Tool(
        name="network_inspect",
        description="OBSERVE the topology task: the nodes and edge format. Call first.",
        input_schema={"type": "object", "properties": {}},
        handler=inspect, phase="observe"))
    registry.register(Tool(
        name="network_propose",
        description=("ACT: submit regulatory edges {source, target, sign(+1/-1)}. Call "
                     "repeatedly to build the draft."),
        input_schema={"type": "object", "properties": {
            "edges": {"type": "array", "items": {"type": "object", "properties": {
                "source": {"type": "string"}, "target": {"type": "string"},
                "sign": {"type": "number"}}}}}},
        handler=propose, phase="act"))
    registry.register(Tool(
        name="network_finalize",
        description="COMMIT the network; scored vs the model. Terminal.",
        input_schema={"type": "object", "properties": {}},
        handler=finalize, phase="evaluate"))
