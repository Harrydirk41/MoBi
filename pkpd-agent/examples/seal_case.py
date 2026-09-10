r"""Seal a benchmark case into an answer-free workspace for a continuous agent.

Produces  <out>/<Model>/workspace  (blanked snapshot + building data + task; the agent runs
here) and  <out>/<Model>/judge  (reference, held-out data, forbidden tokens; grading only).

    python -m examples.seal_case --model Triazolam --out ..\sandbox
    python -m examples.seal_case --all              --out ..\sandbox

Then verify the seal, and point Claude Code at the workspace:

    python -m examples.verify_no_cheat --sandbox ..\sandbox\Triazolam
    cd ..\sandbox\Triazolam\workspace   &&   (drive examples.osp_agent_cli from here)
"""

from __future__ import annotations

import argparse
import glob
import os

from pkpd_agent.engines import osp_sandbox

_HERE = os.path.dirname(os.path.abspath(__file__))
_LIB = os.path.abspath(os.path.join(_HERE, "..", "..", "OSP-PBPK-Model-Library"))


def _resolve(model: str, hard: bool) -> tuple[str, str, str | None] | None:
    """(blanked, input, reference) paths for a model dir, or None if incomplete."""
    d = os.path.join(_LIB, model)
    tag = "hard_blanked" if hard else "blanked"
    blanked = glob.glob(os.path.join(d, "benchmark", f"*{tag}*.json"))
    itag = ".hard.input.json" if hard else ".input.json"
    inp = [p for p in glob.glob(os.path.join(d, "json_input", "*input.json"))
           if p.endswith(itag)]
    ref = glob.glob(os.path.join(d, "json", "*-Model.json"))
    if not blanked or not inp:
        return None
    return blanked[0], inp[0], (ref[0] if ref else None)


def _all_models() -> list[str]:
    out = []
    for d in sorted(glob.glob(os.path.join(_LIB, "*", "benchmark"))):
        model = os.path.basename(os.path.dirname(d))
        out.append(model)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", help="model directory name (e.g. Triazolam)")
    ap.add_argument("--all", action="store_true", help="seal every model with a blanked snapshot")
    ap.add_argument("--hard", action="store_true", help="use the mechanism-discovery (hard) variant")
    ap.add_argument("--out", required=True, help="output root for the sealed sandbox tree")
    args = ap.parse_args()

    models = _all_models() if args.all else ([args.model] if args.model else [])
    if not models:
        ap.error("give --model <Name> or --all")

    n_ok = n_leak = 0
    for m in models:
        paths = _resolve(m, args.hard)
        if not paths:
            print(f"  skip {m}: no {'hard ' if args.hard else ''}blanked snapshot + input")
            continue
        blanked, inp, ref = paths
        out_dir = os.path.join(args.out, m)
        manifest = osp_sandbox.seal_case(blanked, inp, out_dir, reference_path=ref, task_name=m)
        sp = manifest["split"]
        src, red, sealed = (manifest["source_blank_leaks"], manifest["redactions"],
                            manifest["sealed_leaks"])
        n_ok += 1
        bad = sealed > 0
        if bad:
            n_leak += 1
        flag = ("!! SEALED COPY STILL LEAKS" if bad
                else f"sealed (redacted {red} fitted value(s) from source)")
        print(f"{m:22} {flag}  [split {sp['method']}: build {sp['n_building']}, "
              f"held-out {sp['n_verification']}]  -> {out_dir}")

    print(f"\nsealed {n_ok} case(s); {n_leak} sealed copies still leaked (a bug - the "
          f"redactor missed a fitted-value shape).")


if __name__ == "__main__":
    main()
