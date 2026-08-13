"""Tool layer: the agent's action space.

Each engine (pharmpy, OSP, NCA) contributes a set of tools to a shared
registry. A tool is the unit the LLM can call; its handler runs the real (or
mocked) engine work and returns a ``ToolResult``.
"""

from .registry import Tool, ToolRegistry, ToolResult, build_default_registry

__all__ = ["Tool", "ToolRegistry", "ToolResult", "build_default_registry"]
