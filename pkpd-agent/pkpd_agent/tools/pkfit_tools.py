"""Tools for the real pkfit engine (numpy/scipy MLE estimation).

These are the *real* tools - they perform actual computation in-process. The
dataset and fit objects live in the session's artifact store so later tools
(VPC, comparisons) can reuse them.
"""

from __future__ import annotations

from typing import Any

from ..engines.pkfit import PKFitEngine, simulate_dataset
from .registry import Tool, ToolRegistry, ToolResult

_ENGINE = PKFitEngine()


def _fit_id(model: str, covariate: str | None) -> str:
    return f"{model}" + (f"+{covariate}" if covariate else "")


def register_pkfit_tools(registry: ToolRegistry, config) -> None:

    def load_data(args: dict[str, Any], session) -> ToolResult:
        source = args.get("source", "builtin")
        if source == "builtin":
            data, truth = simulate_dataset(n_subjects=int(args.get("n_subjects", 12)))
            session.put("dataset", data)
            session.put("truth", truth)
            return ToolResult.success(
                "loaded builtin dataset (simulated from a known 1-compartment "
                "oral truth with allometric WT on CL, so parameter recovery can "
                "be checked)",
                **data.summary(),
            )
        raise NotImplementedError("only the builtin dataset is wired in this skeleton")

    def nca(args: dict[str, Any], session) -> ToolResult:
        data = session.get("dataset")
        if data is None:
            return ToolResult.error("no dataset loaded - call pkfit_load_data first")
        return ToolResult.success("NCA complete", **_ENGINE.nca(data))

    def fit(args: dict[str, Any], session) -> ToolResult:
        data = session.get("dataset")
        if data is None:
            return ToolResult.error("no dataset loaded - call pkfit_load_data first")
        model = args.get("model", "1cpt_oral")
        covariate = None
        if args.get("covariate_param"):
            covariate = {
                "param": args["covariate_param"],
                "cov": args.get("covariate", "WT"),
                "ref": float(args.get("covariate_ref", 70.0)),
            }
        result = _ENGINE.fit(data, model=model, covariate=covariate)
        fid = _fit_id(result["model"], result.get("covariate"))
        session.put(fid, result)  # full result (incl _theta) for later VPC
        content = {k: v for k, v in result.items() if k != "_theta"}
        content["fit_id"] = fid
        return ToolResult.success(f"fit complete: {fid}", **content)

    def vpc(args: dict[str, Any], session) -> ToolResult:
        data = session.get("dataset")
        fit_obj = session.get(args["fit_id"])
        if data is None or fit_obj is None:
            return ToolResult.error("need a loaded dataset and a valid fit_id")
        return ToolResult.success("VPC complete", **_ENGINE.vpc(data, fit_obj))

    registry.register(Tool(
        name="pkfit_load_data",
        description=(
            "Load a concentration-time dataset for real estimation. source="
            "'builtin' gives a simulated-from-known-truth PK dataset (1-cpt "
            "oral, allometric WT on CL) so you can validate recovery. OBSERVE."
        ),
        input_schema={
            "type": "object",
            "properties": {
                "source": {"type": "string", "enum": ["builtin"]},
                "n_subjects": {"type": "integer"},
            },
        },
        handler=load_data,
        phase="observe",
    ))

    registry.register(Tool(
        name="pkfit_nca",
        description=(
            "Real non-compartmental analysis of the loaded dataset's median "
            "profile: Cmax, Tmax, AUC, terminal half-life, apparent clearance. "
            "The model-free first pass. OBSERVE."
        ),
        input_schema={"type": "object", "properties": {}},
        handler=nca,
        phase="observe",
    ))

    registry.register(Tool(
        name="pkfit_fit",
        description=(
            "Fit a compartmental PK model to the loaded data by maximum "
            "likelihood (naive-pooled, proportional error). model is "
            "'1cpt_oral' or '2cpt_oral'. Optionally add a covariate by giving "
            "covariate_param (e.g. 'CL'), covariate (e.g. 'WT') and "
            "covariate_ref. Returns OFV (-2LL), AIC, BIC, parameter estimates, "
            "relative standard errors, condition number, and whether "
            "minimization succeeded, plus a fit_id. ACT.\n\n"
            "Good practice: compare structural models by AIC (lower is better; "
            "prefer the simpler model unless the fit clearly improves). Test a "
            "covariate with a likelihood-ratio test: keep it only if it drops "
            "OFV by more than 3.84 (chi-square, 1 df, p<0.05)."
        ),
        input_schema={
            "type": "object",
            "properties": {
                "model": {"type": "string", "enum": ["1cpt_oral", "2cpt_oral"]},
                "covariate_param": {"type": "string", "description": "e.g. 'CL' or 'V'"},
                "covariate": {"type": "string", "description": "e.g. 'WT'"},
                "covariate_ref": {"type": "number"},
            },
            "required": ["model"],
        },
        handler=fit,
        phase="act",
    ))

    registry.register(Tool(
        name="pkfit_vpc",
        description=(
            "Monte-Carlo visual predictive check for a fitted model (by "
            "fit_id): fraction of observations inside the 90% prediction "
            "interval. EVALUATE - judge model adequacy."
        ),
        input_schema={
            "type": "object",
            "properties": {"fit_id": {"type": "string"}},
            "required": ["fit_id"],
        },
        handler=vpc,
        phase="evaluate",
    ))
