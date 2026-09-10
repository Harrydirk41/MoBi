r"""Surgically blank the LEAKED Weibull dissolution (absorption) parameters in every
benchmark snapshot - and NOTHING else.

The value-blanking only reset Compound-level fitted parameters; the Formulation-level
Weibull dissolution ('Dissolution time (50% dissolved)', 'Dissolution shape') kept its
fitted ParameterIdentification value, so the agent was handed the reference's fitted
ABSORPTION for free. This resets exactly those params to the neutral naive prior
(30 min half-dissolution, exponential shape=1) and marks them blanked - touching no
other field, so it does not drag in unrelated generator drift (unlike a full
regenerate). The generator (build_benchmark._blank) is also fixed so a future
regenerate stays correct.

    python -m examples.blank_dissolution            # patch all (dry run: --check)
"""

from __future__ import annotations

import argparse
import glob
import json
import os

_HERE = os.path.dirname(os.path.abspath(__file__))
_LIB = os.path.abspath(os.path.join(_HERE, "..", "..", "OSP-PBPK-Model-Library"))

_NEUTRAL = {"Dissolution time (50% dissolved)": 30.0, "Dissolution shape": 1.0,
            "Lag time": 0.0}


def blank_dissolution(snap: dict) -> list[tuple[str, str, float]]:
    """Reset each FITTED formulation dissolution parameter to its neutral prior.
    Returns [(formulation, param, old_value)] for what changed."""
    changed = []
    for form in snap.get("Formulations") or []:
        for p in form.get("Parameters") or []:
            nm = p.get("Name")
            if nm in _NEUTRAL and isinstance(p.get("Value"), (int, float)) \
                    and (p.get("ValueOrigin") or {}).get("Source") == "ParameterIdentification":
                changed.append((form.get("Name"), nm, p["Value"]))
                p["Value"] = _NEUTRAL[nm]
                p["ValueOrigin"] = {"Source": "Unknown",
                                    "Description": "benchmark naive prior (blanked)"}
    return changed


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="report leaks, write nothing")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(_LIB, "*", "benchmark", "*blanked*.json")))
    n = 0
    for f in files:
        try:
            snap = json.load(open(f, encoding="utf-8"))
        except Exception as e:                              # noqa: BLE001
            print(f"  skip {f}: {e}"); continue
        changed = blank_dissolution(snap)
        if not changed:
            continue
        n += 1
        rel = os.path.relpath(f, _LIB)
        for form, nm, old in changed:
            print(f"  {rel}: {form}/{nm}  {old} -> {_NEUTRAL[nm]}")
        if not args.check:
            json.dump(snap, open(f, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    verb = "would blank" if args.check else "blanked"
    print(f"\n{verb} leaked dissolution in {n} snapshot(s) of {len(files)}.")


if __name__ == "__main__":
    main()
