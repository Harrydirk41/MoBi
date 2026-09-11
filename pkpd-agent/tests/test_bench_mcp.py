"""The MCP server is a thin wrapper over the same bench functions. When the mcp SDK is
present, check that the tools are registered and that a no-run tool (`best`) works against a
sealed workspace. Skipped where mcp isn't installed (the wrapper's logic is covered by
test_bench_cli / test_osp_sandbox regardless)."""

from __future__ import annotations

import glob
import os

import pytest

pytest.importorskip("mcp", reason="mcp SDK not installed")

from pkpd_agent.bench import cli, mcp_server  # noqa: E402

_LIB = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..",
                                    "OSP-PBPK-Model-Library"))
_HAS_TRIAZOLAM = bool(glob.glob(os.path.join(_LIB, "Triazolam", "benchmark", "*blanked*.json")))


def test_server_named_and_tools_defined():
    # the module-level functions decorated as tools exist and are callable objects
    for name in ("inspect", "options", "optimize", "sweep", "try_model", "best"):
        assert callable(getattr(mcp_server, name))


@pytest.mark.skipif(not _HAS_TRIAZOLAM, reason="OSP model library not available")
def test_best_tool_reads_session(tmp_path):
    out = str(tmp_path / "sb")
    cli.main(["seal", "--model", "Triazolam", "--out", out])
    ws = os.path.join(out, "Triazolam", "workspace")
    res = mcp_server.best(workspace=ws)          # no PK-Sim needed
    assert res["ok"] is True
    assert res["best_gmfe"] is None              # fresh workspace, nothing fitted yet
    assert res["history"] == []
