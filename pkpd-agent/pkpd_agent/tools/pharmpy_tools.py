"""pharmpy engine adapter + tools (population PK/PD, NLME estimation).

Real calls target pharmpy's documented API:
    pharmpy.modeling.read_model / create_basic_pk_model
    pharmpy.tools.fit / run_amd
    result.ofv / result.parameter_estimates / result.standard_errors

In ``mock`` mode (default) the handlers return deterministic synthetic results
so the loop runs with no pharmpy install and no estimation backend (NONMEM /
nlmixr2). The real code paths are written out and guarded so flipping
``config.mock = False`` exercises them.
"""

from __future__ import annotations

import hashlib
from typing import Any

from .registry import Tool, ToolRegistry, ToolResult


class PharmpyEngine:
    def __init__(self, config) -> None:
        self.config = config

    @property
    def available(self) -> bool:
        try:
            import pharmpy  # noqa: F401
            return True
        except Exception:  # noqa: BLE001
            return False

    # -- helpers --------------------------------------------------------- #
    @staticmethod
    def _seed(*parts: str) -> int:
        h = hashlib.sha256("|".join(parts).encode()).hexdigest()
        return int(h[:8], 16)

    # -- operations ------------------------------------------------------ #
    def load_model(self, path: str) -> dict[str, Any]:
        if self.config.mock or not self.available:
            return {
                "model_id": f"model::{path}",
                "parameters": ["CL", "V", "KA"],
                "structure": "1-compartment oral, first-order absorption",
                "source": "mock",
            }
        from pharmpy.modeling import read_model  # type: ignore
        model = read_model(path)
        return {
            "model_id": model.name,
            "parameters": [p.name for p in model.parameters],
            "structure": str(model.statements.ode_system),
            "source": "pharmpy",
        }

    def fit(self, model_id: str) -> dict[str, Any]:
        if self.config.mock or not self.available:
            s = self._seed("fit", model_id)
            ofv = 1234.5 + (s % 500) / 10.0
            return {
                "model_id": model_id,
                "ofv": round(ofv, 3),
                "minimization_successful": (s % 7) != 0,
                "condition_number": 12.0 + (s % 90),
                "parameter_estimates": {
                    "CL": round(3.0 + (s % 30) / 10, 3),
                    "V": round(30.0 + (s % 100) / 10, 3),
                    "KA": round(1.0 + (s % 20) / 10, 3),
                },
                "relative_standard_errors": {
                    "CL": round(0.08 + (s % 50) / 1000, 3),
                    "V": round(0.11 + (s % 40) / 1000, 3),
                    "KA": round(0.19 + (s % 60) / 1000, 3),
                },
                "source": "mock",
            }
        from pharmpy.modeling import read_model  # type: ignore
        from pharmpy.tools import fit  # type: ignore

        model = read_model(model_id)
        res = fit(model)
        return {
            "model_id": model_id,
            "ofv": float(res.ofv),
            "minimization_successful": bool(
                getattr(res, "minimization_successful", True)
            ),
            "parameter_estimates": {
                k: float(v) for k, v in res.parameter_estimates.items()
            },
            "relative_standard_errors": {
                k: float(v) for k, v in getattr(res, "relative_standard_errors", {}).items()
            },
            "source": "pharmpy",
        }

    def run_amd(self, path: str, model_type: str) -> dict[str, Any]:
        if self.config.mock or not self.available:
            s = self._seed("amd", path, model_type)
            return {
                "final_model": f"amd::{model_type}::{path}",
                "search_summary": {
                    "structural": "1-compartment" if s % 2 else "2-compartment",
                    "iiv_on": ["CL", "V"],
                    "covariates_added": ["WT_on_CL"] if s % 3 else [],
                },
                "ofv": round(1180.0 + (s % 300) / 10, 3),
                "source": "mock",
            }
        from pharmpy.modeling import read_model  # type: ignore
        from pharmpy.tools import run_amd  # type: ignore

        model = read_model(path)
        res = run_amd(model, modeltype=model_type)
        return {
            "final_model": getattr(res, "final_model", str(res)),
            "source": "pharmpy",
        }

    def vpc(self, model_id: str) -> dict[str, Any]:
        if self.config.mock or not self.available:
            s = self._seed("vpc", model_id)
            inside = 90 + (s % 8)
            return {
                "model_id": model_id,
                "pct_observations_within_90_pi": inside,
                "n_bins": 10,
                "source": "mock",
            }
        # Real VPC generation is backend-specific; left as an explicit gap.
        raise NotImplementedError("real VPC generation not wired in this skeleton")


# --------------------------------------------------------------------------- #
# Tool registration
# --------------------------------------------------------------------------- #

def register_pharmpy_tools(registry: ToolRegistry, config) -> None:
    engine = PharmpyEngine(config)

    def load_model(args: dict[str, Any], session) -> ToolResult:
        info = engine.load_model(args["path"])
        session.put(info["model_id"], {"kind": "pharmpy_model", **info})
        return ToolResult.success(f"loaded {info['model_id']}", **info)

    def fit_model(args: dict[str, Any], session) -> ToolResult:
        res = engine.fit(args["model_id"])
        session.put(f"fit::{args['model_id']}", res)
        return ToolResult.success("fit complete", **res)

    def run_amd(args: dict[str, Any], session) -> ToolResult:
        res = engine.run_amd(args["path"], args.get("model_type", "pk"))
        session.put(res["final_model"], {"kind": "amd_result", **res})
        return ToolResult.success("AMD complete", **res)

    def run_vpc(args: dict[str, Any], session) -> ToolResult:
        res = engine.vpc(args["model_id"])
        return ToolResult.success("VPC complete", **res)

    registry.register(Tool(
        name="pharmpy_load_model",
        description=(
            "Load a population PK/PD model file (NONMEM/nlmixr2) into the "
            "session. Use this to OBSERVE an existing model before deciding "
            "how to modify or fit it. Returns model id, parameter names, and "
            "the structural description."
        ),
        input_schema={
            "type": "object",
            "properties": {
                "path": {"type": "string", "description": "Path to the model file."}
            },
            "required": ["path"],
        },
        handler=load_model,
        phase="observe",
    ))

    registry.register(Tool(
        name="pharmpy_fit",
        description=(
            "Estimate a loaded model's parameters against the data by NLME "
            "(FOCE/SAEM via the configured backend). This is an ACT step: it "
            "changes model state. Returns OFV, whether minimization succeeded, "
            "the condition number, parameter estimates and their relative "
            "standard errors."
        ),
        input_schema={
            "type": "object",
            "properties": {
                "model_id": {"type": "string", "description": "Id from pharmpy_load_model."}
            },
            "required": ["model_id"],
        },
        handler=fit_model,
        phase="act",
    ))

    registry.register(Tool(
        name="pharmpy_run_amd",
        description=(
            "Run pharmpy's Automatic Model Development to build structural, "
            "IIV, covariate and error-model components automatically. "
            "model_type is one of pk, pkpd, tmdd, drug_metabolite. An ACT step."
        ),
        input_schema={
            "type": "object",
            "properties": {
                "path": {"type": "string"},
                "model_type": {
                    "type": "string",
                    "enum": ["pk", "pkpd", "tmdd", "drug_metabolite"],
                },
            },
            "required": ["path", "model_type"],
        },
        handler=run_amd,
        phase="act",
    ))

    registry.register(Tool(
        name="pharmpy_vpc",
        description=(
            "Generate a visual predictive check for a fitted model and report "
            "the fraction of observations falling inside the 90% prediction "
            "interval. An EVALUATE step used to judge model adequacy."
        ),
        input_schema={
            "type": "object",
            "properties": {"model_id": {"type": "string"}},
            "required": ["model_id"],
        },
        handler=run_vpc,
        phase="evaluate",
    ))
