r"""The no-cheat seal check: assert a sealed workspace (and optionally an agent transcript)
reveals no answer. Exits non-zero on any leak, so it doubles as a CI gate.

    # verify one sealed sandbox (workspace/ + judge/)
    python -m examples.verify_no_cheat --sandbox ..\sandbox\Triazolam
    python -m examples.verify_no_cheat --sandbox ..\sandbox\Triazolam --transcript session.log

    # structural sweep: every blanked snapshot in the library must be fitted-value-free
    python -m examples.verify_no_cheat --library

The structural check (blanked snapshot carries no ParameterIdentification value) always runs;
the token scans (held-out names/concentrations, reference fitted numbers) need the judge/ tree.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys

from pkpd_agent.engines import osp_sandbox

_HERE = os.path.dirname(os.path.abspath(__file__))
_LIB = os.path.abspath(os.path.join(_HERE, "..", "..", "OSP-PBPK-Model-Library"))


def _verify_sandbox(sandbox: str, transcript: str | None) -> bool:
    ws = os.path.join(sandbox, "workspace")
    jd = os.path.join(sandbox, "judge")
    if not os.path.isdir(ws):
        # allow pointing directly at a workspace dir
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
    """Structural check across every blanked snapshot in the library."""
    ok_all = True
    files = sorted(glob.glob(os.path.join(_LIB, "*", "benchmark", "*blanked*.json")))
    for f in files:
        with open(f, encoding="utf-8") as fh:
            snap = json.load(fh)
        leaks = osp_sandbox.fitted_leaks(snap)
        rel = os.path.relpath(f, _LIB)
        if leaks:
            ok_all = False
            print(f"FAIL  {rel}   {len(leaks)} fitted-value leak(s):")
            for lk in leaks[:6]:
                print(f"        {lk['path']} = {lk['value']} ({(lk['origin'] or {}).get('Source')})")
        else:
            print(f"PASS  {rel}")
    print(f"\n{'ALL CLEAN' if ok_all else 'LEAKS FOUND'} across {len(files)} blanked snapshot(s).")
    return ok_all


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sandbox", help="a sealed sandbox dir (contains workspace/ and judge/), "
                                      "or a workspace dir directly")
    ap.add_argument("--transcript", help="optional agent session log to scan for leaked tokens")
    ap.add_argument("--library", action="store_true",
                    help="structural check on every blanked snapshot in the model library")
    args = ap.parse_args()

    ok = True
    if args.library:
        ok = _verify_library() and ok
    if args.sandbox:
        ok = _verify_sandbox(args.sandbox, args.transcript) and ok
    if not args.library and not args.sandbox:
        ap.error("give --sandbox <dir> and/or --library")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
