"""NCA tool (non-compartmental analysis) - a gap-filling binding.

Neither pharmpy nor OSP does NCA, so this is one of the thin external bindings
the agent needs. Real implementation would shell out to R's PKNCA; here the
mock computes a simple descriptive summary so the loop can use NCA output as
the first-pass observation before model-based work.
"""

from __future__ import annotations

from typing import Any

from .registry import Tool, ToolRegistry, ToolResult


def _trapz_auc(times: list[float], concs: list[float]) -> float:
    auc = 0.0
    for i in range(1, len(times)):
        dt = times[i] - times[i - 1]
        auc += (concs[i] + concs[i - 1]) / 2.0 * dt
    return auc


class NCAEngine:
    def __init__(self, config) -> None:
        self.config = config

    def analyze(self, times: list[float], concs: list[float]) -> dict[str, Any]:
        if not times or len(times) != len(concs):
            raise ValueError("times and concentrations must be non-empty and equal length")
        c_max = max(concs)
        t_max = times[concs.index(c_max)]
        auc = _trapz_auc(times, concs)
        # crude terminal-slope half-life from the last two points, if declining
        thalf = None
        if len(concs) >= 2 and concs[-1] > 0 and concs[-2] > concs[-1]:
            import math
            k = (math.log(concs[-2]) - math.log(concs[-1])) / (times[-1] - times[-2])
            if k > 0:
                thalf = round(math.log(2) / k, 3)
        return {
            "c_max": round(c_max, 4),
            "t_max": round(t_max, 4),
            "auc_trapezoidal": round(auc, 4),
            "t_half_terminal": thalf,
            "n_points": len(times),
            "source": "builtin-nca",
        }


def register_nca_tools(registry: ToolRegistry, config) -> None:
    engine = NCAEngine(config)

    def nca(args: dict[str, Any], session) -> ToolResult:
        res = engine.analyze(args["times"], args["concentrations"])
        return ToolResult.success("NCA complete", **res)

    registry.register(Tool(
        name="nca_analyze",
        description=(
            "Non-compartmental analysis of a concentration-time profile. "
            "OBSERVE step - the model-free first pass (Cmax, Tmax, AUC, "
            "terminal half-life) before any model-based work. pharmpy and OSP "
            "do not do NCA, so this is a dedicated binding."
        ),
        input_schema={
            "type": "object",
            "properties": {
                "times": {"type": "array", "items": {"type": "number"}},
                "concentrations": {"type": "array", "items": {"type": "number"}},
            },
            "required": ["times", "concentrations"],
        },
        handler=nca,
        phase="observe",
    ))
