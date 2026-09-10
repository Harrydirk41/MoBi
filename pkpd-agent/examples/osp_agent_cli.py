r"""The tool bridge for a CONTINUOUSLY-RUNNING modeling agent (e.g. Claude Code).

Instead of the scripted decision loop, this exposes the same OSP loop tools as plain CLI
subcommands the agent calls via its shell. Run from inside a SEALED WORKSPACE (see
examples.seal_case): the workspace holds only `model.blanked.json` + `task.input.json`, so
the agent - confined to this directory - cannot reach the reference model or held-out data.

The agent drives the loop itself, in any order it likes, until it judges the model done:

    python -m examples.osp_agent_cli inspect
    python -m examples.osp_agent_cli options
    python -m examples.osp_agent_cli sweep    --estimate '{"Lipophilicity":[0,5]}'
    python -m examples.osp_agent_cli optimize --estimate '{"Intrinsic clearance":[0.01,100]}' \
                                              --structure '{"add_processes":[...]}'
    python -m examples.osp_agent_cli try      --edits '{"parameters":{...}}'
    python -m examples.osp_agent_cli best      # the best model found so far

State (best-so-far, sweep cache, history) persists across calls in `.osp_session.json`, so
the agent can stop and resume. Every subcommand prints a JSON object to stdout.

Env: PKPD_PKSIM_CLI = path to PKSim.CLI.exe (or pass --pksim).
"""

from __future__ import annotations

import argparse
import contextlib
import json
import os
import sys
from typing import Any

from pkpd_agent.config import AgentConfig
from pkpd_agent.engines.osp_cli import OSPCli
from pkpd_agent.tools.registry import ToolRegistry
from pkpd_agent.tools.osp_loop_tools import register_osp_loop_tools

_SESSION_FILE = ".osp_session.json"


class _PersistentSession:
    """A minimal get/put session backed by a JSON file, so tool state survives across the
    separate CLI processes an external agent spawns. The loop tools use only get/put."""

    def __init__(self, path: str) -> None:
        self._path = path
        self.artifacts: dict[str, Any] = {}
        if os.path.exists(path):
            try:
                with open(path, encoding="utf-8") as fh:
                    self.artifacts = json.load(fh)
            except (json.JSONDecodeError, OSError):
                self.artifacts = {}

    def get(self, name: str, default: Any = None) -> Any:
        return self.artifacts.get(name, default)

    def put(self, name: str, value: Any) -> None:
        self.artifacts[name] = value
        tmp = self._path + ".tmp"
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(self.artifacts, fh, ensure_ascii=False)
        os.replace(tmp, self._path)


def _build_registry(workspace: str, pksim: str | None):
    snap = os.path.join(workspace, "model.blanked.json")
    task = os.path.join(workspace, "task.input.json")
    for p in (snap, task):
        if not os.path.exists(p):
            sys.exit(f"not a sealed workspace (missing {os.path.basename(p)}): {workspace}\n"
                     f"seal one first:  python -m examples.seal_case --model <Name> --out <dir>")
    with open(task, encoding="utf-8") as fh:
        inp = json.load(fh)
    observed = ((inp.get("given_data") or {}).get("clinical_observed_data")) or []

    cfg = AgentConfig(mock=False)
    cfg.stream_optimizer = False        # a machine parses stdout as JSON - no per-eval progress lines
    cli = OSPCli(pksim_cli_path=pksim or os.environ.get("PKPD_PKSIM_CLI") or cfg.pksim_cli_path,
                 timeout_s=cfg.pksim_timeout_s)
    registry = ToolRegistry()
    register_osp_loop_tools(registry, cfg, {
        "cli": cli, "snapshot_path": snap, "observed": observed, "input": inp})
    return registry, cli


def _dispatch(registry, tool: str, call: dict, session):
    """Run a tool, sending its progress prints to stderr so stdout stays pure JSON (a
    Bash-driving agent parses stdout as one JSON object)."""
    with contextlib.redirect_stdout(sys.stderr):
        return registry.dispatch(tool, call, session)


def _emit(result) -> None:
    payload = result.to_content() if hasattr(result, "to_content") else result
    print(json.dumps(payload, ensure_ascii=False, indent=1))
    if hasattr(result, "ok") and not result.ok:
        sys.exit(2)


def _json_arg(s: str | None, what: str) -> dict:
    if not s:
        return {}
    try:
        v = json.loads(s)
    except json.JSONDecodeError as exc:
        sys.exit(f"--{what} is not valid JSON: {exc}")
    if not isinstance(v, dict):
        sys.exit(f"--{what} must be a JSON object")
    return v


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command",
                    choices=["inspect", "options", "optimize", "sweep", "try", "best"])
    ap.add_argument("--workspace", default=".", help="sealed workspace dir (default: cwd)")
    ap.add_argument("--pksim", default=None)
    ap.add_argument("--estimate", help="JSON {param:[lo,hi]} to fit (optimize/sweep)")
    ap.add_argument("--fix", help="JSON {param:value} pinned at literature")
    ap.add_argument("--structure", help="JSON {calculation_methods, processes, add_processes}")
    ap.add_argument("--link-scale", help="JSON [{members:[..],bounds:[lo,hi]}]")
    ap.add_argument("--edits", help="JSON edit spec (try)")
    ap.add_argument("--max-evals", type=int, default=None)
    args = ap.parse_args()

    workspace = os.path.abspath(args.workspace)
    session = _PersistentSession(os.path.join(workspace, _SESSION_FILE))

    if args.command == "best":
        _emit({"ok": True, "best_gmfe": session.get("osp_best_gmfe"),
               "best_edits": session.get("osp_best_edits"),
               "history": session.get("osp_history") or []})
        return

    registry, _cli = _build_registry(workspace, args.pksim)

    if args.command == "inspect":
        _emit(_dispatch(registry, "osp_inspect", {}, session))
    elif args.command == "options":
        _emit(_dispatch(registry, "osp_options", {}, session))
    elif args.command == "try":
        _emit(_dispatch(registry, "osp_try_model",
                        {"edits": _json_arg(args.edits, "edits")}, session))
    elif args.command in ("optimize", "sweep"):
        call: dict[str, Any] = {"estimate": _json_arg(args.estimate, "estimate")}
        if args.fix:
            call["fix"] = _json_arg(args.fix, "fix")
        if args.structure:
            call["structure"] = _json_arg(args.structure, "structure")
        if args.max_evals is not None:
            call["max_evals"] = args.max_evals
        if args.command == "optimize" and args.link_scale:
            try:
                call["link_scale"] = json.loads(args.link_scale)
            except json.JSONDecodeError as exc:
                sys.exit(f"--link-scale is not valid JSON: {exc}")
        tool = "osp_optimize" if args.command == "optimize" else "osp_sweep_methods"
        _emit(registry.dispatch(tool, call, session))


if __name__ == "__main__":
    main()
