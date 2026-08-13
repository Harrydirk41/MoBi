"""The tool registry: names, JSON schemas, handlers, and dispatch.

The registry produces two things:
  * ``to_anthropic_schema()`` -> the ``tools`` array for the Messages API
  * ``dispatch(name, args, session)`` -> executes the handler, returns a
    ``ToolResult``

Tools are tagged with a ``phase`` (observe / act / evaluate) so the loop knows
which ones to run verification gates against (act + evaluate).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Literal

Phase = Literal["observe", "act", "evaluate"]


@dataclass
class ToolResult:
    ok: bool
    data: dict[str, Any] = field(default_factory=dict)
    message: str = ""

    @classmethod
    def success(cls, message: str = "", **data: Any) -> "ToolResult":
        return cls(ok=True, data=data, message=message)

    @classmethod
    def error(cls, message: str, **data: Any) -> "ToolResult":
        return cls(ok=False, data=data, message=message)

    def to_content(self) -> dict[str, Any]:
        """Serializable payload handed back to the model / recorded."""
        return {"ok": self.ok, "message": self.message, **self.data}


# A handler receives parsed arguments and the live session, returns a ToolResult.
Handler = Callable[[dict[str, Any], "ModelingSession"], ToolResult]  # noqa: F821


@dataclass
class Tool:
    name: str
    description: str
    input_schema: dict[str, Any]
    handler: Handler
    phase: Phase = "act"
    destructive: bool = False

    def anthropic_schema(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "description": self.description,
            "input_schema": self.input_schema,
        }


class ToolRegistry:
    def __init__(self) -> None:
        self._tools: dict[str, Tool] = {}

    def register(self, tool: Tool) -> None:
        if tool.name in self._tools:
            raise ValueError(f"duplicate tool: {tool.name}")
        self._tools[tool.name] = tool

    def get(self, name: str) -> Tool:
        return self._tools[name]

    def __contains__(self, name: str) -> bool:
        return name in self._tools

    def __len__(self) -> int:
        return len(self._tools)

    def names(self) -> list[str]:
        return sorted(self._tools)

    def to_anthropic_schema(self) -> list[dict[str, Any]]:
        return [t.anthropic_schema() for t in self._tools.values()]

    def dispatch(self, name: str, args: dict[str, Any], session) -> ToolResult:
        tool = self._tools.get(name)
        if tool is None:
            return ToolResult.error(f"unknown tool: {name}")
        try:
            return tool.handler(args, session)
        except Exception as exc:  # noqa: BLE001 - surface engine errors to the model
            return ToolResult.error(f"{type(exc).__name__}: {exc}")


def build_default_registry(config) -> ToolRegistry:
    """Assemble the full tool surface from every engine."""
    from .pharmpy_tools import register_pharmpy_tools
    from .osp_tools import register_osp_tools
    from .nca_tools import register_nca_tools
    from .pkfit_tools import register_pkfit_tools
    from .nlmixr2_tools import register_nlmixr2_tools

    registry = ToolRegistry()
    register_pkfit_tools(registry, config)     # real, runs-here (numpy/scipy)
    register_nlmixr2_tools(registry, config)   # real NLME via R (nlmixr2)
    register_pharmpy_tools(registry, config)   # popPK/PD via pharmpy (Python/backend)
    register_osp_tools(registry, config)       # mechanistic PBPK/QSP via OSP (R/.NET)
    register_nca_tools(registry, config)       # generic NCA binding
    return registry
