r"""Re-blank the DISTRIBUTION/PERMEABILITY calculation methods in every benchmark blanked snapshot.

The value-blanking only reset fitted PARAMETER VALUES (ParameterIdentification origin); it never
touched CalculationMethods, so the reference model's partition/permeability method leaked into the
blanked snapshot the agent starts from. That contradicts the task objective ("reference methods are
withheld; choose your own") and let the agent inherit the answer's method for free. This resets the
method to PK-Sim's NEUTRAL DEFAULT (what a freshly-created compound has), so the method is a genuine
choice again. Only the method strings are touched - fitted values, structure and processes are left
exactly as the value-blanking produced them.

    python -m examples.reblank_methods            # re-blank the whole library (dry run: --check)
"""

from __future__ import annotations

import argparse
import glob
import json
import os

_HERE = os.path.dirname(os.path.abspath(__file__))
_LIB = os.path.abspath(os.path.join(_HERE, "..", "..", "OSP-PBPK-Model-Library"))

NEUTRAL_PARTITION = "PK-Sim Standard"       # PK-Sim's default for a new compound
NEUTRAL_PERMEABILITY = "PK-Sim Standard"
_PART_PREFIX = "Cellular partition coefficient method - "
_PERM_PREFIX = "Cellular permeability - "


def _reblank_cms(cms, changes):
    for i, m in enumerate(cms or []):
        s = m if isinstance(m, str) else (m.get("Name") or "")
        low = s.lower()
        if "partition coefficient" in low and not s.endswith(NEUTRAL_PARTITION):
            cms[i] = _PART_PREFIX + NEUTRAL_PARTITION
            changes.append(("partition", s, cms[i]))
        elif "permeability" in low and not s.endswith(NEUTRAL_PERMEABILITY):
            cms[i] = _PERM_PREFIX + NEUTRAL_PERMEABILITY
            changes.append(("permeability", s, cms[i]))


def reblank_snapshot(snap: dict) -> list[tuple[str, str, str]]:
    """Reset each compound's partition/permeability method to the neutral default - on BOTH the
    top-level Compounds AND each Simulation's own Compounds copy (PK-Sim uses the simulation's copy,
    so leaving it holds the answer's method). Returns the list of (kind, old, new) changes."""
    changes = []
    for c in snap.get("Compounds", []) or []:
        if isinstance(c.get("CalculationMethods"), list):
            _reblank_cms(c["CalculationMethods"], changes)
    for s in snap.get("Simulations", []) or []:
        for sc in s.get("Compounds") or []:
            if isinstance(sc.get("CalculationMethods"), list):
                _reblank_cms(sc["CalculationMethods"], changes)
    return changes


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="report leaks, write nothing")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(_LIB, "*", "benchmark", "*blanked*.json")))
    n_leaked = 0
    for f in files:
        try:
            snap = json.load(open(f, encoding="utf-8"))
        except Exception as e:                              # noqa: BLE001
            print(f"  skip {f}: {e}"); continue
        changes = reblank_snapshot(snap)
        if not changes:
            continue
        n_leaked += 1
        rel = os.path.relpath(f, _LIB)
        for kind, old, new in changes:
            print(f"  {rel}: {kind} {old.split(' - ')[-1]} -> {new.split(' - ')[-1]}")
        if not args.check:
            json.dump(snap, open(f, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    verb = "would re-blank" if args.check else "re-blanked"
    print(f"\n{verb} the method in {n_leaked} snapshot(s) that leaked a non-default method "
          f"(of {len(files)} blanked files).")


if __name__ == "__main__":
    main()
