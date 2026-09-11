r"""Grade a model a continuous agent built in a sealed sandbox - the judge-side step.

After the agent (Claude Code) stops, its best model lives in
``workspace/.osp_session.json``. This re-runs that model AND the reference on the HELD-OUT
studies (whose concentrations live only in ``judge/``, never in the workspace) and grades the
agent RELATIVE to the reference: a valid model predicts held-out data about as well as the
expert model, so PASS = agent held-out GMFE <= 1.25x the reference's.

    set PKPD_PKSIM_CLI=C:\...\PKSim.CLI.exe
    python -m examples.judge_case --sandbox ..\sandbox\Triazolam

Writes ``judge/score.json`` and prints the verdict. Needs PK-Sim (re-runs both models).
"""

from __future__ import annotations

import argparse
import json
import os
import sys

from pkpd_agent.config import AgentConfig
from pkpd_agent.engines.osp_cli import OSPCli
from pkpd_agent.engines import osp_score


def _run_edits_from_best(best_edits: dict) -> dict:
    """Reproduce the fitted model: merge FIXED params (e.g. GFR=0) into parameters and carry
    the structure, exactly as the report's re-run does, so the graded model is the one the
    agent actually built."""
    best_edits = best_edits or {}
    fixed = dict(best_edits.get("fix") or {})
    run_edits = {k: v for k, v in best_edits.items() if k != "fix"}
    run_edits["parameters"] = {**fixed, **(best_edits.get("parameters") or {})}
    return run_edits


def _heldout_gmfe(cli, snapshot: str, edits, held_obs: list) -> "float | None":
    res = cli.build_and_run(snapshot, edits=edits)
    if not res.get("ok"):
        return None
    # heuristic matching for BOTH models against the same held-out observed set (the agent's
    # sealed snapshot has no OutputMappings linkage, so score both the same way - symmetric).
    pred, _ = osp_score.map_predictions(res.get("profiles", []), held_obs)
    return osp_score.score_fit(held_obs, pred)["overall"]["gmfe"] if pred else None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sandbox", required=True, help="sealed sandbox dir (workspace/ + judge/)")
    ap.add_argument("--pksim", default=None)
    args = ap.parse_args()

    ws = os.path.join(args.sandbox, "workspace")
    jd = os.path.join(args.sandbox, "judge")
    ws_snap = os.path.join(ws, "model.blanked.json")
    session_f = os.path.join(ws, ".osp_session.json")
    ref = os.path.join(jd, "reference.json")
    heldf = os.path.join(jd, "heldout.observed.json")
    for p in (ws_snap, ref, heldf):
        if not os.path.exists(p):
            sys.exit(f"missing {p} - seal the case and let the agent run first.")

    held_obs = json.load(open(heldf, encoding="utf-8"))
    if not held_obs:
        sys.exit("no held-out data for this case (split was 'none-too-few') - "
                 "grade on building fit instead.")
    best_edits = {}
    if os.path.exists(session_f):
        best_edits = (json.load(open(session_f, encoding="utf-8")) or {}).get("osp_best_edits") or {}
    if not best_edits:
        sys.exit("the agent recorded no model in workspace/.osp_session.json - "
                 "did it run osp_agent_cli optimize/sweep?")

    cfg = AgentConfig(mock=False)
    cli = OSPCli(pksim_cli_path=args.pksim or os.environ.get("PKPD_PKSIM_CLI") or cfg.pksim_cli_path,
                 timeout_s=cfg.pksim_timeout_s)
    if not cli.pksim_cli_path or not os.path.exists(cli.pksim_cli_path):
        sys.exit(f"PKSim.CLI not found at {cli.pksim_cli_path!r}; set PKPD_PKSIM_CLI or --pksim.")

    print(f"grading {os.path.basename(os.path.normpath(args.sandbox))} on "
          f"{len(held_obs)} held-out dataset(s)...", file=sys.stderr)
    a_ho = _heldout_gmfe(cli, ws_snap, _run_edits_from_best(best_edits), held_obs)
    r_ho = _heldout_gmfe(cli, ref, None, held_obs)
    ratio = (a_ho / r_ho) if (a_ho and r_ho) else None
    verdict = {
        "held_out_datasets": len(held_obs),
        "agent_heldout_gmfe": a_ho,
        "reference_heldout_gmfe": r_ho,
        "ratio_agent_over_reference": ratio,
        "pass": bool(ratio is not None and ratio <= 1.25),
        "pass_rule": "agent held-out GMFE <= 1.25 x reference held-out GMFE",
        "best_edits": best_edits,
    }
    with open(os.path.join(jd, "score.json"), "w", encoding="utf-8") as fh:
        json.dump(verdict, fh, ensure_ascii=False, indent=1)
    print(json.dumps(verdict, ensure_ascii=False, indent=1))
    sys.exit(0 if verdict["pass"] else 3)


if __name__ == "__main__":
    main()
