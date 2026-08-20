"""LLM-loop tools for isolated SIGN prediction (the weakest Stage-1 layer, on its own).

In the full topology task, signs and edges are entangled. Here we hand the agent the
model's TRUE edges (unsigned) and ask only: is each activating or inhibiting? Scored as
sign accuracy vs the model, against the majority-class baseline (guessing the more common
sign for every edge). This isolates whether the LLM can get the direction right when it
already knows the interaction exists.

  * ``sign_inspect``  (observe) - the true (unsigned) edges to sign, and the objective.
  * ``sign_predict``  (act)     - submit {source, target, sign} for the edges.
  * ``sign_finalize`` (evaluate)- score sign accuracy vs the model. Terminal.
"""

from __future__ import annotations

from ..engines import ra_network as N
from .registry import Tool, ToolRegistry, ToolResult


def register_ra_sign_loop_tools(registry: ToolRegistry, config, ctx: dict) -> None:
    """ctx: {truth: list[Edge]} - signs hidden; only (source, target) shown."""
    truth: list = ctx["truth"]

    def inspect(args: dict, session) -> ToolResult:
        edges = [{"source": e.source, "target": e.target} for e in truth]
        return ToolResult.success(
            "RA sign task: these are the model's real regulatory edges, unsigned. For each, "
            "decide whether the source ACTIVATES (+1) or INHIBITS (-1) the target's "
            "secretion/proliferation/influx.",
            objective="predict the sign (+1 activate / -1 inhibit) of every edge",
            note="the model includes negative self-feedback edges (a cell down-regulating "
                 "its OWN cytokine - a saturation/boundedness term), plus the expected "
                 "anti-inflammatory edges (TGF-b, IL-10, Treg). Most edges are activating.",
            edges=edges)

    def predict(args: dict, session) -> ToolResult:
        acc = dict(session.get("sign_pred") or {})
        n = 0
        for it in args.get("edges") or []:
            src = N.canon_node(str(it.get("source", "")))
            tgt = N.canon_node(str(it.get("target", "")))
            if src is None or tgt is None:
                continue
            s = it.get("sign", 1)
            if isinstance(s, str):
                s = -1 if s.lower() in ("-", "-1", "inhibit", "suppress", "anti", "down") else 1
            acc[(src, tgt)] = 1 if s >= 0 else -1
            n += 1
        session.put("sign_pred", acc)
        return ToolResult.success(
            f"recorded {n} sign(s); {len(acc)}/{len(truth)} edges signed.",
            n_signed=len(acc), n_remaining=len(truth) - len(acc))

    def finalize(args: dict, session) -> ToolResult:
        pred = session.get("sign_pred") or {}
        if not pred:
            return ToolResult.error("predict signs before finalizing.")
        sc = N.score_signs(pred, truth)
        session.put("sign_final", sc)
        return ToolResult.success(
            f"sign accuracy {sc['accuracy']} ({sc['correct']}/{sc['n']}); majority-class "
            f"baseline {sc['majority_baseline']} ({sc['frac_positive']} of edges are "
            f"activating). Beats majority: {sc['beats_majority']}.",
            **sc)

    registry.register(Tool(
        name="sign_inspect",
        description=("OBSERVE the sign task: the model's real edges (unsigned) and the "
                     "objective. Call first."),
        input_schema={"type": "object", "properties": {}},
        handler=inspect, phase="observe"))

    registry.register(Tool(
        name="sign_predict",
        description=("ACT: submit signs for the edges as a list of {source, target, "
                     "sign(+1/-1)}. Sign every edge before finalizing."),
        input_schema={"type": "object", "properties": {
            "edges": {"type": "array", "items": {"type": "object", "properties": {
                "source": {"type": "string"}, "target": {"type": "string"},
                "sign": {"type": "number"}}}}}},
        handler=predict, phase="act"))

    registry.register(Tool(
        name="sign_finalize",
        description=("COMMIT your signs; scored as accuracy vs the model, against the "
                     "majority-class baseline. Terminal - call once."),
        input_schema={"type": "object", "properties": {}},
        handler=finalize, phase="evaluate"))
