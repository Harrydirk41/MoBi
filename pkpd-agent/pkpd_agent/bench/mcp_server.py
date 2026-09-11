r"""MCP server exposing the pkpd-bench agent-loop tools as NATIVE tools.

This is a thin wrapper over exactly the same functions the `pkpd-bench` CLI uses (one source
of truth) — so a coding agent can drive the model build through structured tool calls instead
of Bash strings, and can run with NO shell access (a tighter sandbox around the answers).

Run it as an MCP stdio server; point a client (Claude Code, etc.) at it via `.mcp.json`:

    {
      "mcpServers": {
        "pkpd-bench": {
          "command": "python",
          "args": ["-m", "pkpd_agent.bench.mcp_server"],
          "env": {
            "PKPD_WORKSPACE": "C:/sandbox/Triazolam/workspace",
            "PKPD_PKSIM_CLI": "C:/Program Files/Open Systems Pharmacology/PK-Sim 12.3/PKSim.CLI.exe"
          }
        }
      }
    }

Needs `pip install mcp` (>=1.0). The setup/grading side (seal / verify / judge) stays on the
`pkpd-bench` CLI — those are the orchestrator's job, not the sandboxed agent's.

CRITICAL: an MCP stdio server speaks the protocol on stdout, so nothing else may print there.
Every tool runs through `_dispatch`, which redirects the engines' progress prints to stderr —
keep it that way, or a stray print corrupts the wire protocol.
"""

from __future__ import annotations

import os
from typing import Any

try:                                                        # mcp >= 2 (FastMCP -> MCPServer)
    from mcp.server.mcpserver import MCPServer as _Server
except ImportError:                                         # mcp < 2
    from mcp.server.fastmcp import FastMCP as _Server        # type: ignore

from .cli import _PersistentSession, _build_registry, _dispatch, _SESSION_FILE

server = _Server("pkpd-bench")


def _ws(workspace: str) -> str:
    return os.path.abspath(workspace or os.environ.get("PKPD_WORKSPACE") or ".")


def _pksim() -> str | None:
    return os.environ.get("PKPD_PKSIM_CLI")


def _session(workspace: str) -> _PersistentSession:
    return _PersistentSession(os.path.join(_ws(workspace), _SESSION_FILE))


def _call(tool: str, args: dict, workspace: str) -> dict:
    reg = _build_registry(_ws(workspace), _pksim())
    return _dispatch(reg, tool, args, _session(workspace)).to_content()


@server.tool()
def inspect(workspace: str = "") -> dict:
    """Read the task: objective, known biology, physchem priors, the building clinical data,
    and the current model (each parameter's value_status: given = measured, trust it;
    placeholder = unknown, determine it). Call this first. workspace defaults to $PKPD_WORKSPACE."""
    return _call("osp_inspect", {}, workspace)


@server.tool()
def options(workspace: str = "") -> dict:
    """The authoritative action space: every editable parameter with its role/range, the legal
    distribution & permeability methods, the molecules the model expresses, and the mechanisms
    you may add. Read this for exact names and legal choices before deciding structure."""
    return _call("osp_options", {}, workspace)


@server.tool()
def optimize(estimate: dict[str, Any], fix: dict[str, Any] | None = None,
             structure: dict[str, Any] | None = None, max_evals: int | None = None,
             workspace: str = "") -> dict:
    """Fit the parameters you choose. estimate = {parameter: [lo, hi]} (2-4 identifiable,
    uncertain ones). fix = {parameter: value} pins literature values. structure =
    {calculation_methods, processes, add_processes}. Returns the fitted values, GMFE, per-route
    bias, params_at_bound, and identifiability recommendations. Watch the IVIVE trap: a given
    in-vitro Vmax does NOT mean clearance is set — the in-vivo kcat usually still needs fitting."""
    args: dict[str, Any] = {"estimate": estimate}
    if fix:
        args["fix"] = fix
    if structure:
        args["structure"] = structure
    if max_evals is not None:
        args["max_evals"] = max_evals
    return _call("osp_optimize", args, workspace)


@server.tool()
def sweep(estimate: dict[str, Any], fix: dict[str, Any] | None = None,
          structure: dict[str, Any] | None = None, max_evals: int | None = None,
          workspace: str = "") -> dict:
    """Choose the distribution method deterministically: try every partition x permeability
    method, re-fit your physchem (estimate = {param: [lo, hi]}) under each, and adopt the best.
    Use once when distribution / Vd is off. Pass the same structure you use with optimize so the
    mechanism is held fixed while only the methods vary."""
    args: dict[str, Any] = {"estimate": estimate}
    if fix:
        args["fix"] = fix
    if structure:
        args["structure"] = structure
    if max_evals is not None:
        args["max_evals"] = max_evals
    return _call("osp_sweep_methods", args, workspace)


@server.tool()
def try_model(edits: dict[str, Any], workspace: str = "") -> dict:
    """Apply an edit spec and score it without a fit. edits =
    {parameters:{name:value}, calculation_methods:{...}, processes:{Molecule:true/false}}.
    Returns GMFE and per-route bias. Use to probe a specific setting; use optimize to fit."""
    return _call("osp_try_model", {"edits": edits}, workspace)


@server.tool()
def best(workspace: str = "") -> dict:
    """The best model found so far this session (GMFE, edits, history), read from
    workspace/.osp_session.json. No model run."""
    s = _session(workspace)
    return {"ok": True, "best_gmfe": s.get("osp_best_gmfe"),
            "best_edits": s.get("osp_best_edits"), "history": s.get("osp_history") or []}


def main() -> None:
    server.run("stdio")


if __name__ == "__main__":
    main()
