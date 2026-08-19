"""LLM-loop tools for the RA Stage-4 CALIBRATION task: fit PD parameters to data.

The inverse of the virtual-trial tools. Here the protocol is FIXED (a given drug
arm) and the unknown is a model PD parameter: the agent estimates it so the model
reproduces an OBSERVED clinical response - the same parameter-identification loop
as the OSP PK fitting, but on the RA QSP model.

  * ``fit_inspect`` (observe) - the fixed drug arm, the parameter(s) to calibrate
    (name, unit, meaning, plausible range), and the observed target ACR.
  * ``fit_try``     (act)     - set the parameter value(s), run the arm across the
    virtual population, and return the predicted ACR and its error vs the target.
  * ``fit_finalize``(evaluate)- commit the fitted parameter set to be scored.

The model structure and Vpop are GIVEN and fixed; the unknown is the parameter.
"""

from __future__ import annotations

from typing import Any

from ..engines.simbiology import SimBiologyEngine
from ..engines import osp_ra_trial
from .registry import Tool, ToolRegistry, ToolResult


def register_ra_fit_loop_tools(registry: ToolRegistry, config, ctx: dict) -> None:
    """ctx: {sb, vpop, arm (dose string), target (observed ACR dict), fit_params
    (list of param names from osp_ra_trial.FIT_PARAMS), limit, stop_time}."""
    sb: SimBiologyEngine = ctx["sb"]
    vpop: str = ctx["vpop"]
    arm: str = ctx["arm"]
    target: dict = ctx["target"]
    fit_params: list = ctx.get("fit_params") or ["KD_TCZ"]
    limit: int = int(ctx.get("limit") or 50)
    stop_time: float = float(ctx.get("stop_time") or 700.0)
    enable_optimize: bool = bool(ctx.get("enable_optimize", True))

    def _param_specs() -> list[dict]:
        out = []
        for name in fit_params:
            p = osp_ra_trial.FIT_PARAMS.get(name) or {}
            out.append({
                "name": name, "unit": p.get("unit"), "meaning": p.get("meaning"),
                "plausible_range": p.get("search_range"),
                "search_on_log_scale": p.get("log_scale", False),
            })
        return out

    # -- observe -------------------------------------------------------- #
    def inspect(args: dict, session) -> ToolResult:
        return ToolResult.success(
            "RA calibration task: estimate the model PD parameter(s) so the fixed "
            "drug arm reproduces the OBSERVED clinical response.",
            objective=f"fit {', '.join(fit_params)} to match the observed ACR",
            fixed_arm=arm,
            readout="second-line response among MTX non-responders (day 600): "
                    "ACR20/50/70 in %",
            observed_target=target,
            parameters_to_fit=_param_specs(),
            edit_spec_help={
                "overrides": "{param_name: value} - the parameter value(s) to try, "
                             "e.g. {\"KD_TCZ\": 2.5e-12}. Search a wide range first, "
                             "then home in; a parameter marked log-scale moves the "
                             "response over orders of magnitude, so step by factors."},
        )

    def _run(overrides: dict) -> dict:
        spec = osp_ra_trial.build_override_spec(overrides)
        r = sb.run_vpop(vpop, dose=arm, stop_time=stop_time, limit=limit,
                        param_overrides=spec)
        summary = osp_ra_trial.summarize_run(r)
        return {"spec": spec, "summary": summary,
                "matlab_log": (r.get("matlab_log") or "")}

    # -- act ------------------------------------------------------------ #
    def try_fit(args: dict, session) -> ToolResult:
        overrides = args.get("overrides") or {}
        if not overrides:
            return ToolResult.error(
                "no overrides given - pass {overrides:{param:value}}. "
                "Call fit_inspect for the parameter names and ranges.")
        out = _run(overrides)
        sl = out["summary"]["second_line"]
        pred = {k: sl.get(k) for k in ("ACR20", "ACR50", "ACR70")}
        score = osp_ra_trial.score_flagship(pred, target)
        mae = score.get("mae_pp")

        hist = session.get("fit_history") or []
        hist.append({"overrides": overrides, "predicted": pred, "mae": mae,
                     "n_MTX_IR": sl.get("n_MTX_IR")})
        session.put("fit_history", hist)
        best = session.get("fit_best_mae")
        if mae is not None and (best is None or mae < best):
            session.put("fit_best_mae", mae)
            session.put("fit_best_overrides", overrides)

        return ToolResult.success(
            f"{out['spec']}  ->  MTX-IR n={sl.get('n_MTX_IR')}: "
            f"ACR20 {pred.get('ACR20')}% ACR50 {pred.get('ACR50')}% "
            f"ACR70 {pred.get('ACR70')}%  |  target {target}  |  MAE {mae} pp "
            f"(best so far {session.get('fit_best_mae')})",
            overrides=overrides, predicted=pred, target=target,
            mae_pp=mae, per_endpoint=score.get("per_endpoint"),
            best_mae_so_far=session.get("fit_best_mae"), iteration=len(hist))

    # -- act via NUMERICAL optimizer (agent sets up, scipy minimizes) --- #
    def optimize(args: dict, session) -> ToolResult:
        param = args.get("param") or (fit_params[0] if fit_params else None)
        lo, hi = args.get("lo"), args.get("hi")
        if not param or lo is None or hi is None:
            return ToolResult.error(
                "give {param, lo, hi} - the parameter to fit and its search bounds. "
                "Call fit_inspect for names and plausible ranges.")
        p = osp_ra_trial.FIT_PARAMS.get(param) or {}
        log = bool(args.get("log", p.get("log_scale", True)))
        cache: dict = {}

        def evaluate(val):
            out = _run({param: val})
            sl = out["summary"]["second_line"]
            pred = {k: sl.get(k) for k in ("ACR20", "ACR50", "ACR70")}
            mae = osp_ra_trial.score_flagship(pred, target).get("mae_pp")
            cache[val] = pred
            return mae

        res = osp_ra_trial.numeric_fit_1d(
            evaluate, float(lo), float(hi), log=log,
            max_evals=int(args.get("max_evals") or 20))
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
            f"fitted {fitted:.4g}, ACR MAE {res['error']} pp in {res['n_evals']} "
            f"evaluations. Predicted {pred}.",
            param=param, fitted=fitted, mae_pp=res["error"],
            n_evals=res["n_evals"], predicted=pred,
            profile=prof, best_mae_so_far=session.get("fit_best_mae"))

    # -- evaluate (commit) ---------------------------------------------- #
    def finalize(args: dict, session) -> ToolResult:
        overrides = args.get("overrides") or session.get("fit_best_overrides")
        if not overrides:
            return ToolResult.error("nothing to finalize - run fit_try first.")
        hist = session.get("fit_history") or []
        match = next((h for h in reversed(hist) if h.get("overrides") == overrides), None)
        if match is None:
            return ToolResult.error(
                "that exact override set was not run - fit_try it first, then finalize.")
        refs = {n: (osp_ra_trial.FIT_PARAMS.get(n) or {}).get("reference")
                for n in overrides}
        score = osp_ra_trial.score_fit(match["predicted"], target, overrides, refs)
        session.put("fit_final", {"overrides": overrides,
                                  "predicted": match["predicted"], "score": score})
        return ToolResult.success(
            f"committed fit {overrides}: ACR MAE {score.get('acr_mae_pp')} pp vs "
            f"observed; parameter vs literature {score.get('parameters')}",
            overrides=overrides, score=score)

    registry.register(Tool(
        name="fit_inspect",
        description=(
            "OBSERVE the calibration task: the fixed drug arm, the model PD "
            "parameter(s) to estimate (name, unit, biological meaning, plausible "
            "range), and the OBSERVED target ACR the fitted model must reproduce. "
            "Call this first."),
        input_schema={"type": "object", "properties": {}},
        handler=inspect, phase="observe"))

    registry.register(Tool(
        name="fit_try",
        description=(
            "ACT: run the fixed arm at ONE explicit parameter value and return the "
            "predicted ACR + error vs target. Use this to PROBE - check a value, "
            "profile the response around an optimum, or confirm a direction - NOT to "
            "hand-run a minimization. For the actual fit, prefer fit_optimize."),
        input_schema={"type": "object", "properties": {
            "overrides": {"type": "object",
                          "description": "{param_name: value}, e.g. {\"KD_TCZ\": 2.5e-12}"}}},
        handler=try_fit, phase="act"))

    if enable_optimize:
        registry.register(Tool(
            name="fit_optimize",
            description=(
                "ACT: fit a parameter NUMERICALLY. You set up the problem - the "
                "parameter, its search bounds [lo,hi], and whether it is log-scale - "
                "and a bounded numerical optimizer does the minimization, returning "
                "the fitted value, the ACR error there, the evaluation count, and the "
                "(value,error) profile. This is the right tool for the actual fit; "
                "your job is choosing what to fit and interpreting the result "
                "(identifiability from the profile, whether residuals are structural, "
                "how the value compares to the literature)."),
            input_schema={"type": "object", "properties": {
                "param": {"type": "string", "description": "parameter name to fit"},
                "lo": {"type": "number", "description": "lower search bound"},
                "hi": {"type": "number", "description": "upper search bound"},
                "log": {"type": "boolean", "description": "search on log scale (default true)"},
                "max_evals": {"type": "number", "description": "max model evaluations (default 20)"}},
                "required": ["param", "lo", "hi"]},
            handler=optimize, phase="act"))

    registry.register(Tool(
        name="fit_finalize",
        description=(
            "COMMIT the parameter set you are calibrating to (the one you recommend, "
            "which must already have been run with fit_try). Scores the fit against "
            "the observed ACR and reports each fitted value against its literature "
            "reference. Call once before finishing."),
        input_schema={"type": "object", "properties": {
            "overrides": {"type": "object"}}},
        handler=finalize, phase="evaluate"))
