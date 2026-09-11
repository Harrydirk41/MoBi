r"""pkpd-bench — the single agent-facing CLI for the sealed PBPK benchmark.

One command, one clean surface, runs from ANY directory (the logic lives in the installed
package, so there is no `-m examples` / cwd coupling). A continuous coding agent drives the
middle group; a human runs seal / verify / judge around it.

    pkpd-bench seal    --model Triazolam --out ../sandbox      # build workspace + judge tree
    pkpd-bench verify  --sandbox ../sandbox/Triazolam          # prove the wall (CI gate)

    # the agent, working inside ../sandbox/Triazolam/workspace:
    pkpd-bench inspect  --workspace .
    pkpd-bench options  --workspace .
    pkpd-bench sweep    --workspace . --estimate '{"Lipophilicity":[0,5]}'
    pkpd-bench optimize --workspace . --estimate '{"<param>":[lo,hi]}' --structure '{...}'
    pkpd-bench try      --workspace . --edits '{"parameters":{...}}'
    pkpd-bench best     --workspace .

    pkpd-bench judge   --sandbox ../sandbox/Triazolam          # grade held-out prediction

Every subcommand prints one JSON object on stdout (progress goes to stderr). If installed as
an entry point it's `pkpd-bench ...`; otherwise `python -m pkpd_agent.bench.cli ...`.
Env: PKPD_PKSIM_CLI = path to PKSim.CLI.exe (or pass --pksim).
"""

from __future__ import annotations

import argparse
import contextlib
import glob
import json
import os
import sys
from typing import Any

_LIB = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                    "OSP-PBPK-Model-Library"))
_SESSION_FILE = ".osp_session.json"


# --------------------------------------------------------------------------- #
# shared: persistent session + tool registry for the agent-loop commands
# --------------------------------------------------------------------------- #

class _PersistentSession:
    """A minimal get/put session backed by a JSON file, so tool state survives across the
    separate processes an external agent spawns. The loop tools use only get/put."""

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
    from ..config import AgentConfig
    from ..engines.osp_cli import OSPCli
    from ..tools.registry import ToolRegistry
    from ..tools.osp_loop_tools import register_osp_loop_tools

    snap = os.path.join(workspace, "model.blanked.json")
    task = os.path.join(workspace, "task.input.json")
    for p in (snap, task):
        if not os.path.exists(p):
            sys.exit(f"not a sealed workspace (missing {os.path.basename(p)}): {workspace}\n"
                     f"seal one first:  pkpd-bench seal --model <Name> --out <dir>")
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
    return registry


def _dispatch(registry, tool: str, call: dict, session):
    """Run a tool, sending its progress prints to stderr so stdout stays pure JSON."""
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


# --------------------------------------------------------------------------- #
# seal
# --------------------------------------------------------------------------- #

def _resolve(model: str, hard: bool):
    d = os.path.join(_LIB, model)
    tag = "hard_blanked" if hard else "blanked"
    blanked = glob.glob(os.path.join(d, "benchmark", f"*{tag}*.json"))
    itag = ".hard.input.json" if hard else ".input.json"
    inp = [p for p in glob.glob(os.path.join(d, "json_input", "*input.json")) if p.endswith(itag)]
    ref = glob.glob(os.path.join(d, "json", "*-Model.json"))
    if not blanked or not inp:
        return None
    return blanked[0], inp[0], (ref[0] if ref else None)


def _all_models() -> list[str]:
    return [os.path.basename(os.path.dirname(d))
            for d in sorted(glob.glob(os.path.join(_LIB, "*", "benchmark")))]


def cmd_seal(args) -> None:
    from ..engines import osp_sandbox
    models = _all_models() if args.all else ([args.model] if args.model else [])
    if not models:
        sys.exit("give --model <Name> or --all")
    n_ok = n_leak = 0
    for m in models:
        paths = _resolve(m, args.hard)
        if not paths:
            print(f"  skip {m}: no {'hard ' if args.hard else ''}blanked snapshot + input")
            continue
        blanked, inp, ref = paths
        out_dir = os.path.join(args.out, m)
        mf = osp_sandbox.seal_case(blanked, inp, out_dir, reference_path=ref, task_name=m)
        sp, sealed = mf["split"], mf["sealed_leaks"]
        n_ok += 1
        bad = sealed > 0
        n_leak += bad
        flag = ("!! SEALED COPY STILL LEAKS" if bad
                else f"sealed (redacted {mf['redactions']} fitted value(s))")
        print(f"{m:22} {flag}  [split {sp['method']}: build {sp['n_building']}, "
              f"held-out {sp['n_verification']}]  -> {out_dir}")
    print(f"\nsealed {n_ok} case(s); {n_leak} sealed copies still leaked.")


# --------------------------------------------------------------------------- #
# verify
# --------------------------------------------------------------------------- #

def _verify_sandbox(sandbox: str, transcript: str | None) -> bool:
    from ..engines import osp_sandbox
    ws = os.path.join(sandbox, "workspace")
    jd = os.path.join(sandbox, "judge")
    if not os.path.isdir(ws):
        ws, jd = sandbox, (jd if os.path.isdir(jd) else None)
    rep = osp_sandbox.verify_workspace(ws, judge_dir=jd if jd and os.path.isdir(jd) else None,
                                       transcript_path=transcript)
    tag = os.path.basename(os.path.normpath(sandbox))
    ninfo = len(rep.get("info") or [])
    if rep["ok"]:
        extra = f"  ({ninfo} held-out study name(s) present - info only)" if ninfo else ""
        print(f"PASS  {tag}   checks={rep['checks']}{extra}")
    else:
        print(f"FAIL  {tag}   {len(rep['leaks'])} leak(s):")
        for lk in rep["leaks"][:12]:
            print(f"        [{lk['kind']}] {lk['where']}: {lk['detail']}")
    return rep["ok"]


def _verify_library() -> bool:
    from ..engines import osp_sandbox
    ok_all = True
    files = sorted(glob.glob(os.path.join(_LIB, "*", "benchmark", "*blanked*.json")))
    for f in files:
        with open(f, encoding="utf-8") as fh:
            leaks = osp_sandbox.fitted_leaks(json.load(fh))
        rel = os.path.relpath(f, _LIB)
        if leaks:
            ok_all = False
            print(f"FAIL  {rel}   {len(leaks)} fitted-value leak(s)")
        else:
            print(f"PASS  {rel}")
    print(f"\n{'ALL CLEAN' if ok_all else 'LEAKS FOUND'} across {len(files)} blanked snapshot(s).")
    return ok_all


def cmd_verify(args) -> None:
    ok = True
    if args.library:
        ok = _verify_library() and ok
    if args.sandbox:
        ok = _verify_sandbox(args.sandbox, args.transcript) and ok
    if not args.library and not args.sandbox:
        sys.exit("give --sandbox <dir> and/or --library")
    sys.exit(0 if ok else 1)


# --------------------------------------------------------------------------- #
# judge
# --------------------------------------------------------------------------- #

def _run_edits_from_best(best_edits: dict) -> dict:
    best_edits = best_edits or {}
    fixed = dict(best_edits.get("fix") or {})
    run_edits = {k: v for k, v in best_edits.items() if k != "fix"}
    run_edits["parameters"] = {**fixed, **(best_edits.get("parameters") or {})}
    return run_edits


def _heldout_gmfe(cli, snapshot: str, edits, held_obs: list):
    from ..engines import osp_score
    res = cli.build_and_run(snapshot, edits=edits)
    if not res.get("ok"):
        return None
    pred, _ = osp_score.map_predictions(res.get("profiles", []), held_obs)
    return osp_score.score_fit(held_obs, pred)["overall"]["gmfe"] if pred else None


def cmd_judge(args) -> None:
    from ..config import AgentConfig
    from ..engines.osp_cli import OSPCli
    ws = os.path.join(args.sandbox, "workspace")
    jd = os.path.join(args.sandbox, "judge")
    ws_snap = os.path.join(ws, "model.blanked.json")
    session_f = os.path.join(ws, _SESSION_FILE)
    ref = os.path.join(jd, "reference.json")
    heldf = os.path.join(jd, "heldout.observed.json")
    for p in (ws_snap, ref, heldf):
        if not os.path.exists(p):
            sys.exit(f"missing {p} - seal the case and let the agent run first.")
    held_obs = json.load(open(heldf, encoding="utf-8"))
    if not held_obs:
        sys.exit("no held-out data for this case - grade on building fit instead.")
    best_edits = {}
    if os.path.exists(session_f):
        best_edits = (json.load(open(session_f, encoding="utf-8")) or {}).get("osp_best_edits") or {}
    if not best_edits:
        sys.exit("the agent recorded no model in workspace/.osp_session.json.")

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
        "agent_heldout_gmfe": a_ho, "reference_heldout_gmfe": r_ho,
        "ratio_agent_over_reference": ratio,
        "pass": bool(ratio is not None and ratio <= 1.25),
        "pass_rule": "agent held-out GMFE <= 1.25 x reference held-out GMFE",
        "best_edits": best_edits,
    }
    with open(os.path.join(jd, "score.json"), "w", encoding="utf-8") as fh:
        json.dump(verdict, fh, ensure_ascii=False, indent=1)
    print(json.dumps(verdict, ensure_ascii=False, indent=1))
    sys.exit(0 if verdict["pass"] else 3)


# --------------------------------------------------------------------------- #
# agent-loop commands (inspect / options / optimize / sweep / try / best)
# --------------------------------------------------------------------------- #

def cmd_loop(args) -> None:
    workspace = os.path.abspath(args.workspace)
    session = _PersistentSession(os.path.join(workspace, _SESSION_FILE))

    if args.cmd == "best":
        _emit({"ok": True, "best_gmfe": session.get("osp_best_gmfe"),
               "best_edits": session.get("osp_best_edits"),
               "history": session.get("osp_history") or []})
        return

    registry = _build_registry(workspace, args.pksim)
    if args.cmd == "inspect":
        _emit(_dispatch(registry, "osp_inspect", {}, session))
    elif args.cmd == "options":
        _emit(_dispatch(registry, "osp_options", {}, session))
    elif args.cmd == "try":
        _emit(_dispatch(registry, "osp_try_model", {"edits": _json_arg(args.edits, "edits")}, session))
    else:                                       # optimize / sweep
        call: dict[str, Any] = {"estimate": _json_arg(args.estimate, "estimate")}
        if args.fix:
            call["fix"] = _json_arg(args.fix, "fix")
        if args.structure:
            call["structure"] = _json_arg(args.structure, "structure")
        if args.max_evals is not None:
            call["max_evals"] = args.max_evals
        if args.cmd == "optimize" and args.link_scale:
            try:
                call["link_scale"] = json.loads(args.link_scale)
            except json.JSONDecodeError as exc:
                sys.exit(f"--link-scale is not valid JSON: {exc}")
        tool = "osp_optimize" if args.cmd == "optimize" else "osp_sweep_methods"
        _emit(_dispatch(registry, tool, call, session))


# --------------------------------------------------------------------------- #
# parser
# --------------------------------------------------------------------------- #

def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="pkpd-bench", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command", required=True)

    p = sub.add_parser("seal", help="build a sealed workspace + judge tree")
    p.add_argument("--model"); p.add_argument("--all", action="store_true")
    p.add_argument("--hard", action="store_true"); p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_seal)

    p = sub.add_parser("verify", help="assert a sealed workspace leaks no answer")
    p.add_argument("--sandbox"); p.add_argument("--transcript"); p.add_argument("--library", action="store_true")
    p.set_defaults(func=cmd_verify)

    p = sub.add_parser("judge", help="grade the agent's model on held-out prediction")
    p.add_argument("--sandbox", required=True); p.add_argument("--pksim")
    p.set_defaults(func=cmd_judge)

    # agent-loop commands share --workspace / --pksim
    def loop_parser(name, help_):
        q = sub.add_parser(name, help=help_)
        q.add_argument("--workspace", default="."); q.add_argument("--pksim")
        q.set_defaults(func=cmd_loop, cmd=name)
        return q

    loop_parser("inspect", "read the task, priors, data, and current model")
    loop_parser("options", "the authoritative action space")
    loop_parser("best", "the best model found so far")
    loop_parser("try", "apply an edit spec, run, and score").add_argument("--edits")
    for nm, hp in (("optimize", "fit chosen parameters"),
                   ("sweep", "sweep distribution methods, re-fitting physchem under each")):
        q = loop_parser(nm, hp)
        q.add_argument("--estimate"); q.add_argument("--fix"); q.add_argument("--structure")
        q.add_argument("--max-evals", type=int, default=None)
        if nm == "optimize":
            q.add_argument("--link-scale")
    return ap


def main(argv=None) -> None:
    args = build_parser().parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
