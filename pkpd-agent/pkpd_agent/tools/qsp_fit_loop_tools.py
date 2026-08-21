"""Model-agnostic loop tools for the CALIBRATION task: fit PD parameters to data.

The inverse of the virtual-trial tools. Here the protocol is FIXED (a given drug arm)
and the unknown is a model PD parameter: the agent estimates it so the model
reproduces an OBSERVED clinical response. Parameter catalog, arm and target come off
the QSPTaskConfig.

  * ``fit_inspect``  (observe) - the fixed arm, the parameter(s) to calibrate, target.
  * ``fit_try``      (act)     - set values, run the arm, return predicted response + error.
  * ``fit_optimize`` (act)     - fit ONE parameter numerically (agent sets up, scipy minimizes).
  * ``fit_finalize`` (evaluate)- commit the fitted parameter set to be scored.
"""

from __future__ import annotations

from ..engines.simbiology import SimBiologyEngine
from ..engines import qsp_tasks
from ..engines.qsp_config import QSPTaskConfig
from .registry import Tool, ToolRegistry, ToolResult

_ENDPOINTS = ("ACR20", "ACR50", "ACR70")


def register_qsp_fit_loop_tools(registry: ToolRegistry, config, ctx: dict) -> None:
    """ctx: {cfg, sb, vpop, arm, target, fit_params?, limit, stop_time,
    enable_optimize?}."""
    cfg: QSPTaskConfig = ctx["cfg"]
    sb: SimBiologyEngine = ctx["sb"]
    vpop: str = ctx["vpop"]
    arm: str = ctx["arm"]
    target: dict = ctx["target"]
    fit_params: list = ctx.get("fit_params") or list(cfg.fit_params.keys())
    endpoints = [k for k in _ENDPOINTS if k in target] or list(target.keys())
    limit: int = int(ctx.get("limit") or 50)
    stop_time: float = float(ctx.get("stop_time") or 700.0)
    enable_optimize: bool = bool(ctx.get("enable_optimize", True))

    def _param_specs() -> list[dict]:
        out = []
        for name in fit_params:
            p = cfg.fit_params.get(name) or {}
            out.append({"name": name, "unit": p.get("unit"), "meaning": p.get("meaning"),
                        "plausible_range": p.get("search_range"),
                        "search_on_log_scale": p.get("log_scale", False)})
        return out

    def inspect(args: dict, session) -> ToolResult:
        return ToolResult.success(
            f"{cfg.name} calibration task: estimate the model PD parameter(s) so the "
            "fixed drug arm reproduces the OBSERVED clinical response.",
            objective=f"fit {', '.join(fit_params)} to match the observed response",
            fixed_arm=arm,
            readout="second-line response among the subgroup (second readout)",
            observed_target=target,
            parameters_to_fit=_param_specs(),
            edit_spec_help={"overrides": "{param_name: value} - the value(s) to try. "
                            "Search a wide range first, then home in; a log-scale "
                            "parameter moves the response over orders of magnitude."})

    def _run(overrides: dict) -> dict:
        spec = qsp_tasks.build_override_spec(overrides)
        r = sb.run_vpop(vpop, dose=arm, stop_time=stop_time, limit=limit,
                        param_overrides=spec)
        return {"spec": spec, "summary": cfg.summarize_run(r),
                "matlab_log": (r.get("matlab_log") or "")}

    def try_fit(args: dict, session) -> ToolResult:
        overrides = args.get("overrides") or {}
        if not overrides:
            return ToolResult.error("no overrides given - pass {overrides:{param:value}}. "
                                    "Call fit_inspect for the parameter names and ranges.")
        out = _run(overrides)
        sl = out["summary"]["second_line"]
        pred = {k: sl.get(k) for k in endpoints}
        score = qsp_tasks.score_flagship(pred, target)
        mae = score.get("mae_pp")

        hist = session.get("fit_history") or []
        hist.append({"overrides": overrides, "predicted": pred, "mae": mae,
                     "n_subgroup": sl.get("n_subgroup")})
        session.put("fit_history", hist)
        best = session.get("fit_best_mae")
        if mae is not None and (best is None or mae < best):
            session.put("fit_best_mae", mae)
            session.put("fit_best_overrides", overrides)

        return ToolResult.success(
            f"{out['spec']}  ->  subgroup n={sl.get('n_subgroup')}: {pred}  |  "
            f"target {target}  |  MAE {mae} pp (best so far {session.get('fit_best_mae')})",
            overrides=overrides, predicted=pred, target=target, mae_pp=mae,
            per_endpoint=score.get("per_endpoint"),
            best_mae_so_far=session.get("fit_best_mae"), iteration=len(hist))

    def optimize(args: dict, session) -> ToolResult:
        param = args.get("param") or (fit_params[0] if fit_params else None)
        lo, hi = args.get("lo"), args.get("hi")
        if not param or lo is None or hi is None:
            return ToolResult.error("give {param, lo, hi} - the parameter to fit and its "
                                    "search bounds. Call fit_inspect for plausible ranges.")
        p = cfg.fit_params.get(param) or {}
        log = bool(args.get("log", p.get("log_scale", True)))
        max_evals = min(int(args.get("max_evals") or 12), 14)
        cache: dict = {}

        def evaluate(val):
            out = _run({param: val})
            sl = out["summary"]["second_line"]
            pred = {k: sl.get(k) for k in endpoints}
            mae = qsp_tasks.score_flagship(pred, target).get("mae_pp")
            cache[val] = pred
            print(f"     [opt {len(cache):>2}/{max_evals}] {param}={val:.3g} -> "
                  f"{endpoints[0]} {pred.get(endpoints[0])} MAE {mae} pp", flush=True)
            return mae

        res = qsp_tasks.numeric_fit_1d(evaluate, float(lo), float(hi), log=log,
                                       max_evals=max_evals)
        fitted = res["fitted"]
        pred = cache.get(fitted, {})

        hist = session.get("fit_history") or []
        hist.append({"overrides": {param: fitted}, "predicted": pred,
                     "mae": res["error"], "via": "numeric"})
        session.put("fit_history", hist)
        best = session.get("fit_best_mae")
        if res["error"] is not None and (best is None or res["error"] < best):
            session.put("fit_best_mae", res["error"])
            session.put("fit_best_overrides", {param: fitted})

        prof = [f"{t['value']:.3g}:{t['error']}" for t in res["trace"]]
        return ToolResult.success(
            f"numeric fit of {param} in [{lo:g},{hi:g}] ({'log' if log else 'lin'}): "
            f"fitted {fitted:.4g}, MAE {res['error']} pp in {res['n_evals']} "
            f"evaluations. Predicted {pred}.",
            param=param, fitted=fitted, mae_pp=res["error"], n_evals=res["n_evals"],
            predicted=pred, profile=prof, best_mae_so_far=session.get("fit_best_mae"))

    def finalize(args: dict, session) -> ToolResult:
        overrides = args.get("overrides") or session.get("fit_best_overrides")
        if not overrides:
            return ToolResult.error("nothing to finalize - run fit_try first.")
        hist = session.get("fit_history") or []
        match = next((h for h in reversed(hist) if h.get("overrides") == overrides), None)
        if match is None:
            return ToolResult.error("that exact override set was not run - "
                                    "fit_try it first, then finalize.")
        refs = {n: (cfg.fit_params.get(n) or {}).get("reference") for n in overrides}
        score = qsp_tasks.score_fit(match["predicted"], target, overrides, refs)
        session.put("fit_final", {"overrides": overrides,
                                  "predicted": match["predicted"], "score": score})
        return ToolResult.success(
            f"committed fit {overrides}: MAE {score.get('acr_mae_pp')} pp vs observed; "
            f"parameter vs literature {score.get('parameters')}",
            overrides=overrides, score=score)

    registry.register(Tool(
        name="fit_inspect",
        description=("OBSERVE the calibration task: the fixed drug arm, the model PD "
                     "parameter(s) to estimate (name, unit, meaning, plausible range), "
                     "and the OBSERVED target the fitted model must reproduce. Call first."),
        input_schema={"type": "object", "properties": {}},
        handler=inspect, phase="observe"))
    registry.register(Tool(
        name="fit_try",
        description=("ACT: run the fixed arm at ONE explicit parameter value and return "
                     "the predicted response + error vs target. Use this to PROBE - not "
                     "to hand-run a minimization. For the actual fit, prefer fit_optimize."),
        input_schema={"type": "object", "properties": {
            "overrides": {"type": "object", "description": "{param_name: value}"}}},
        handler=try_fit, phase="act"))
    if enable_optimize:
        registry.register(Tool(
            name="fit_optimize",
            description=("ACT: fit a parameter NUMERICALLY. You set up the problem - the "
                         "parameter, its bounds [lo,hi], and whether it is log-scale - and "
                         "a bounded optimizer does the minimization, returning the fitted "
                         "value, the error, the evaluation count and the (value,error) "
                         "profile. Your job is choosing what to fit and interpreting the "
                         "result (identifiability, comparison to the literature)."),
            input_schema={"type": "object", "properties": {
                "param": {"type": "string"}, "lo": {"type": "number"},
                "hi": {"type": "number"}, "log": {"type": "boolean"},
                "max_evals": {"type": "number"}}, "required": ["param", "lo", "hi"]},
            handler=optimize, phase="act"))
    registry.register(Tool(
        name="fit_finalize",
        description=("COMMIT the parameter set you are calibrating to (already run with "
                     "fit_try). Scores the fit against the observed response and reports "
                     "each fitted value against its literature reference. Call once."),
        input_schema={"type": "object", "properties": {"overrides": {"type": "object"}}},
        handler=finalize, phase="evaluate"))
