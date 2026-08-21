"""Model-agnostic loop tools for every Stage-1 benchmark, driven by a QSPModel.

One register function per benchmark; each reads its cast / answer key off the QSPModel
(derived from network.json + a spec) instead of hardcoded RA vocab. Together with
qsp_topology_loop_tools these let the whole suite run on any QSP model.
"""

from __future__ import annotations

from ..engines import ra_network as N
from ..engines import ra_params as RP
from ..engines.qsp_model import QSPModel
from .registry import Tool, ToolRegistry, ToolResult


# ----------------------------------------------------------------- scope -- #
def register_qsp_scope_loop_tools(registry: ToolRegistry, config, ctx: dict) -> None:
    model: QSPModel = ctx["model"]

    def inspect(a, s):
        return ToolResult.success(
            f"{model.spec.name} scope task: propose the cast (cells and mediators) this "
            "QSP model should include. Precision is scored - a good model is parsimonious.",
            objective="list the cell types and mediators that belong in the model",
            readout=model.spec.readout_name,
            format="ONE entity per list entry, by its standard short name (e.g. 'Th1', "
                   "'IL-6', 'macrophage') - do not group several into one string or add "
                   "parenthetical aliases")

    def propose(a, s):
        acc = list(s.get("qsp_scope") or []) + [str(n) for n in (a.get("nodes") or [])]
        s.put("qsp_scope", acc)
        return ToolResult.success(f"{len(set(map(str.lower, acc)))} distinct so far.")

    def finalize(a, s):
        names = s.get("qsp_scope") or []
        if not names:
            return ToolResult.error("propose the cast first.")
        sc = model.score_node_set(names, model.nodes)
        s.put("qsp_scope_final", sc)
        return ToolResult.success(f"F1 {sc['f1']} (P {sc['precision']} R {sc['recall']}); "
                                  f"{sc['hit']}/{sc['n_truth']} nodes.", **sc)

    _reg3(registry, "scope", inspect, propose, finalize, "nodes")


# ------------------------------------------------------------------ signs -- #
def register_qsp_sign_loop_tools(registry: ToolRegistry, config, ctx: dict) -> None:
    model: QSPModel = ctx["model"]

    def inspect(a, s):
        return ToolResult.success(
            f"{model.spec.name} sign task: for each real edge, is it activating (+1) or "
            "inhibiting (-1)?",
            edges=[{"source": e.source, "target": e.target} for e in model.edges])

    def predict(a, s):
        acc = dict(s.get("qsp_signs") or {})
        for it in a.get("edges") or []:
            src = model.resolve(str(it.get("source", "")))
            tgt = model.resolve(str(it.get("target", "")))
            if src and tgt:
                v = it.get("sign", 1)
                if isinstance(v, str):
                    v = -1 if v.lower() in ("-", "-1", "inhibit", "suppress", "down") else 1
                acc[f"{src}|{tgt}"] = 1 if v >= 0 else -1
        s.put("qsp_signs", acc)
        return ToolResult.success(f"{len(acc)}/{len(model.edges)} signed.")

    def finalize(a, s):
        raw = s.get("qsp_signs") or {}
        if not raw:
            return ToolResult.error("predict signs first.")
        pred = {tuple(k.split("|")): v for k, v in raw.items()}
        sc = N.score_signs(pred, model.edges)
        s.put("qsp_sign_final", sc)
        return ToolResult.success(f"accuracy {sc['accuracy']} (majority {sc['majority_baseline']}, "
                                  f"beats {sc['beats_majority']}).", **sc)

    _reg3(registry, "sign", inspect, predict, finalize, "edges", act_name="sign_predict")


# ---------------------------------------------------------------- readout -- #
def register_qsp_readout_loop_tools(registry: ToolRegistry, config, ctx: dict) -> None:
    model: QSPModel = ctx["model"]

    def inspect(a, s):
        return ToolResult.success(
            f"{model.spec.name} readout task: which model nodes is the readout "
            f"({model.spec.readout_name}) computed from?",
            objective="list the MODEL NODES (cells/mediators) the readout is a direct "
                      "function of",
            hint="the readout is a formula over the model's own state variables (cell "
                 "densities and mediators), NOT the clinical instrument's inputs - so "
                 "propose model nodes (which cells/mediators drive severity), not clinical "
                 "sub-scores like joint counts, CRP, or patient global assessment",
            format="one model node per entry, standard short name")

    def propose(a, s):
        acc = list(s.get("qsp_readout") or []) + [str(n) for n in (a.get("nodes") or [])]
        s.put("qsp_readout", acc)
        return ToolResult.success(f"{len(set(acc))} distinct.")

    def finalize(a, s):
        names = s.get("qsp_readout") or []
        if not names:
            return ToolResult.error("propose the drivers first.")
        sc = model.score_node_set(names, model.readout_drivers)
        s.put("qsp_readout_final", sc)
        return ToolResult.success(f"F1 {sc['f1']} (P {sc['precision']} R {sc['recall']}); "
                                  f"{sc['hit']}/{sc['n_truth']} drivers.", **sc)

    _reg3(registry, "readout", inspect, propose, finalize, "nodes")


# ----------------------------------------------------------------- params -- #
def register_qsp_params_loop_tools(registry: ToolRegistry, config, ctx: dict) -> None:
    model: QSPModel = ctx["model"]

    def inspect(a, s):
        # a large model can have hundreds of params; present an evenly-spaced, deterministic
        # subset so the agent can actually finish (order-of-magnitude scoring is over
        # whatever it predicts). Prefer the units where physiology helps.
        ps = model.params
        if len(ps) > 150:
            phys = [p for p in ps if p.physiological()]
            rest = [p for p in ps if not p.physiological()]
            step = max(1, len(rest) // (150 - len(phys))) if len(phys) < 150 else 1
            ps = phys + rest[::step]
        return ToolResult.success(
            f"{model.spec.name} parameter task: predict each parameter's value from its "
            f"name and units (order-of-magnitude scoring). {len(ps)} shown.",
            parameters=[{"name": p.name, "units": p.units} for p in ps])

    def estimate(a, s):
        acc = dict(s.get("qsp_params") or {})
        acc.update(RP.clean_predictions(a.get("predictions") or []))
        s.put("qsp_params", acc)
        return ToolResult.success(f"{len(acc)}/{len(model.params)} estimated.")

    def finalize(a, s):
        pred = dict(s.get("qsp_params") or {})
        if not pred:
            return ToolResult.error("estimate values first.")
        sc = RP.score_params(pred, model.params)
        s.put("qsp_params_final", sc)
        ph = sc["physiological"]
        return ToolResult.success(
            f"physiological median {ph.get('median_log10_err')}, beats fair baseline "
            f"{sc['beats_physiological_baseline']}.", **sc)

    _reg3(registry, "param", inspect, estimate, finalize, "predictions", act_name="param_estimate")


# ------------------------------------------------------------ sensitivity -- #
def register_qsp_sensitivity_loop_tools(registry: ToolRegistry, config, ctx: dict) -> None:
    model: QSPModel = ctx["model"]

    def inspect(a, s):
        return ToolResult.success(
            f"{model.spec.name} sensitivity task: from the pool, rank the parameters that "
            f"most drive {model.spec.readout_name}.",
            candidate_pool=model.sensitivity_pool())

    def rank(a, s):
        picks = [str(p) for p in (a.get("ranked") or [])]
        s.put("qsp_sens", picks)
        return ToolResult.success(f"recorded {len(picks)}.")

    def finalize(a, s):
        picks = s.get("qsp_sens") or []
        if not picks:
            return ToolResult.error("rank the parameters first.")
        sc = model.score_sensitivity(picks)
        s.put("qsp_sens_final", sc)
        return ToolResult.success(f"recall {sc['recall']} ({sc['hit']}), random "
                                  f"{sc['random_baseline_recall']}, beats {sc['beats_random']}.",
                                  **sc)

    _reg3(registry, "sens", inspect, rank, finalize, "ranked", act_name="sens_rank")


# ----------------------------------------------------------------- helper -- #
def _reg3(registry, name, inspect, act, finalize, act_key, act_name=None):
    act_name = act_name or f"{name}_propose"
    registry.register(Tool(name=f"{name}_inspect", description=f"OBSERVE the {name} task.",
                           input_schema={"type": "object", "properties": {}},
                           handler=inspect, phase="observe"))
    registry.register(Tool(
        name=act_name, description=f"ACT: submit for the {name} task.",
        input_schema={"type": "object", "properties": {
            act_key: {"type": "array", "items": {"type": "object"}
                      if act_key in ("edges", "predictions") else {"type": "string"}}}},
        handler=act, phase="act"))
    registry.register(Tool(name=f"{name}_finalize", description=f"COMMIT the {name} task.",
                           input_schema={"type": "object", "properties": {}},
                           handler=finalize, phase="evaluate"))
