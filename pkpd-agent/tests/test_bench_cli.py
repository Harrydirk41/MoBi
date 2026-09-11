"""The agent-facing CLI (pkpd-bench) — parser wiring + a seal→verify roundtrip.

The whole point of this surface is that a coding agent drives ONE clean command. These
tests pin that the command exists with the expected subcommands and that seal + verify work
end to end without PK-Sim (the parts that need the engine are covered in test_osp_sandbox)."""

from __future__ import annotations

import glob
import json
import os

import pytest

from pkpd_agent.bench import cli

_LIB = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..",
                                    "OSP-PBPK-Model-Library"))
_HAS_TRIAZOLAM = bool(glob.glob(os.path.join(_LIB, "Triazolam", "benchmark", "*blanked*.json")))


def test_parser_exposes_the_expected_subcommands():
    parser = cli.build_parser()
    # argparse stores subcommand names on the subparsers action
    sub = [a for a in parser._actions if hasattr(a, "choices") and a.choices]
    names = set(sub[0].choices) if sub else set()
    assert {"seal", "verify", "judge", "inspect", "options",
            "optimize", "sweep", "try", "best"} <= names


def test_run_edits_from_best_merges_fix_into_parameters():
    be = {"parameters": {"kcat": 17.5}, "fix": {"GFR fraction": 0.0},
          "calculation_methods": {"partition": "Rodgers and Rowland"}}
    run = cli._run_edits_from_best(be)
    assert run["parameters"] == {"GFR fraction": 0.0, "kcat": 17.5}
    assert run["calculation_methods"]["partition"] == "Rodgers and Rowland"
    assert "fix" not in run


@pytest.mark.skipif(not _HAS_TRIAZOLAM, reason="OSP model library not available")
def test_seal_then_verify_roundtrip(tmp_path):
    out = str(tmp_path / "sb")
    # seal (cmd_seal prints and returns; no exit on success)
    cli.main(["seal", "--model", "Triazolam", "--out", out])
    sandbox = os.path.join(out, "Triazolam")
    assert os.path.exists(os.path.join(sandbox, "workspace", "model.blanked.json"))
    assert os.path.exists(os.path.join(sandbox, "workspace", "PROMPT.md"))

    # verify exits 0 when the seal is clean
    with pytest.raises(SystemExit) as ei:
        cli.main(["verify", "--sandbox", sandbox])
    assert ei.value.code == 0

    # planting a held-out concentration token into a workspace file makes verify exit non-zero
    tokens = json.load(open(os.path.join(sandbox, "judge", "forbidden.json")))["heldout_conc_tokens"]
    assert tokens, "expected distinctive held-out concentration tokens"
    with open(os.path.join(sandbox, "workspace", "leak.txt"), "w") as fh:
        fh.write(f"a held-out value slipped in: {tokens[0]}")
    with pytest.raises(SystemExit) as ei:
        cli.main(["verify", "--sandbox", sandbox])
    assert ei.value.code == 1
