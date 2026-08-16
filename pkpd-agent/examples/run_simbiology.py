r"""Drive the Vantage RA QSP model headlessly: load, inspect, simulate a drug.

Proves the SimBiologyEngine end to end on the real model - loads the .sbproj,
lists its species/doses, runs a baseline simulation and a drug-dosed simulation,
and reports each state's endpoint value (baseline vs dosed) so we can identify the
RA disease-activity readout among the 59 states. All array data crosses via files.

    python -m examples.run_simbiology --sbproj "C:\...\Vantage RA QSP Model v1.0.sbproj"
    python -m examples.run_simbiology --sbproj "..." --dose TCZ8mgkg_Q4W_IV_t200
"""

from __future__ import annotations

import argparse
import os
import sys
import traceback

from pkpd_agent.engines.simbiology import SimBiologyEngine

_LOG = None


def log(msg: str) -> None:
    print(msg, flush=True)
    if _LOG is not None:
        try:
            _LOG.write(msg + "\n")
            _LOG.flush()
        except Exception:                              # noqa: BLE001
            pass


def _endpoint(res: dict, name: str):
    col = (res.get("columns") or {}).get(name) or []
    for v in reversed(col):
        if v is not None:
            return v
    return None


def main() -> None:
    global _LOG
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sbproj", required=True)
    ap.add_argument("--dose", default=None, help="dose name (default: first TCZ/one listed)")
    ap.add_argument("--stop-time", type=float, default=0.0)
    args = ap.parse_args()

    try:
        sys.stdout.reconfigure(line_buffering=True)
    except Exception:                                  # noqa: BLE001
        pass
    _LOG = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "run_simbiology.log"), "w", encoding="utf-8")

    sb = SimBiologyEngine()
    try:
        log("-> start MATLAB engine")
        sb.start()
        log("<- engine started")

        log(f"-> load project {args.sbproj}")
        sb.load_project(args.sbproj)
        log("<- loaded")

        log("-> model_info")
        info = sb.model_info()
        log(f"<- model: {info.get('name')} | "
            f"{info.get('nSpecies')} species, {info.get('nParameters')} params, "
            f"{info.get('nReactions')} reactions")
        doses = info.get("doses") or []
        species = info.get("species") or []
        log(f"   doses ({len(doses)}): {doses}")
        log(f"   species ({len(species)}): {species}")

        # pick a dose to test
        dose = args.dose
        if not dose:
            dose = next((d for d in doses if "TCZ" in d), doses[0] if doses else None)
        log(f"   testing dose: {dose}")

        log("-> simulate baseline (no dose)")
        base = sb.simulate(stop_time=args.stop_time)
        log(f"<- baseline: {len(base['time'])} timepoints, {base.get('n_columns')} columns")

        log(f"-> simulate dosed ({dose})")
        dosed = sb.simulate(dose=dose or "", stop_time=args.stop_time)
        log(f"<- dosed: {len(dosed['time'])} timepoints")

        # report endpoint baseline vs dosed for each state, sorted by |log fold|
        log("\nendpoint (last timepoint) baseline vs dosed, largest change first:")
        rows = []
        for name in species:
            b, d = _endpoint(base, name), _endpoint(dosed, name)
            if isinstance(b, (int, float)) and isinstance(d, (int, float)) and b not in (0, None):
                fold = (d / b) if b else None
                rows.append((name, b, d, fold))
        rows.sort(key=lambda r: abs((r[3] or 1) - 1), reverse=True)
        for name, b, d, fold in rows[:15]:
            fs = f"{fold:.3f}x" if isinstance(fold, (int, float)) else "n/a"
            log(f"   {name:36} base={b:.4g}  dosed={d:.4g}  ({fs})")

        log("\n=== RUN END (reached the end cleanly) ===")
    except Exception as exc:                            # noqa: BLE001
        log(f"[FAIL] {type(exc).__name__}: {exc}")
        log(traceback.format_exc())
    finally:
        sb.stop()
        if _LOG is not None:
            _LOG.close()


if __name__ == "__main__":
    main()
