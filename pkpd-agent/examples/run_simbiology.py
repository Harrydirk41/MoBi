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


def _fmt(v):
    return f"{v:.4g}" if isinstance(v, (int, float)) else "n/a"


def _interp(times, vals, t):
    pts = [(a, b) for a, b in zip(times, vals)
           if isinstance(a, (int, float)) and isinstance(b, (int, float))]
    if not pts:
        return None
    if t <= pts[0][0]:
        return pts[0][1]
    if t >= pts[-1][0]:
        return pts[-1][1]
    for i in range(1, len(pts)):
        x0, y0 = pts[i - 1]
        x1, y1 = pts[i]
        if x0 <= t <= x1:
            return y0 if x1 == x0 else y0 + (t - x0) / (x1 - x0) * (y1 - y0)
    return pts[-1][1]


def _endpoint(res: dict, name: str):
    col = (res.get("columns") or {}).get(name) or []
    for v in reversed(col):
        if v is not None:
            return v
    return None


def _stats(res: dict, name: str):
    """(endpoint, min, max) of a state over the whole trajectory."""
    col = [v for v in ((res.get("columns") or {}).get(name) or [])
           if isinstance(v, (int, float))]
    if not col:
        return None, None, None
    return col[-1], min(col), max(col)


# the model's built-in clinical readouts (disease activity + ACR responses)
_CLINICAL = ["DAS28_CRP", "DAS28_BL", "ACR_Perc", "ACR20", "ACR50", "ACR70",
             "Remission", "Response", "Response_2"]


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

        log(f"   (baseline logged {base.get('n_columns')} columns, "
            f"dosed logged {dosed.get('n_columns')} columns)")

        # sample key states over a shared time grid (baseline vs dosed at MATCHED
        # times) - the states start at 0 and ramp to disease steady state, so we
        # need the trajectory, not min/max. Confirms (a) the drug reaches the
        # synovium (dosing works) and (b) whether it bends IL-6 / DAS28.
        tmax = max((base["time"] or [0])[-1], (dosed["time"] or [0])[-1])
        grid = [round(tmax * f, 3) for f in (0, .05, .1, .15, .2, .3, .5, .75, 1.0)]
        key = [s for s in ("DAS28_CRP", "IL6", "TNFa", "ACR_Perc",
                           "TCZDrug_Synovium", "TCZDrug_Central") if s in species]
        log(f"\ntrajectory at matched times (base -> dosed), tmax={tmax:g}:")
        log("   " + "state".ljust(18) + "  " + "  ".join(f"t={t:g}".rjust(11) for t in grid))
        for name in key:
            cells = []
            for t in grid:
                b = _interp(base["time"], base["columns"].get(name, []), t)
                d = _interp(dosed["time"], dosed["columns"].get(name, []), t)
                cells.append(f"{_fmt(b)}->{_fmt(d)}".rjust(11))
            log(f"   {name:18}  " + "  ".join(cells))

        # biggest movers: dosed vs baseline at the SAME sample times (t>0)
        log("\nlargest dosed-vs-baseline changes at matched times (excludes t=0):")
        rows = []
        for name in species:
            best = None
            for t in grid[1:]:
                b = _interp(base["time"], base["columns"].get(name, []), t)
                d = _interp(dosed["time"], dosed["columns"].get(name, []), t)
                if isinstance(b, (int, float)) and isinstance(d, (int, float)) and b not in (0, None):
                    fold = d / b
                    if best is None or abs(fold - 1) > abs(best[2] - 1):
                        best = (t, b, fold, d)
            if best and abs(best[2] - 1) > 1e-3:
                rows.append((name, best[0], best[1], best[3], best[2]))
        rows.sort(key=lambda r: abs(r[4] - 1), reverse=True)
        for name, t, b, d, fold in rows[:12]:
            log(f"   {name:36} @t={t:g}: {b:.4g} -> {d:.4g}  ({fold:.3f}x)")
        if not rows:
            log("   (no state differs between baseline and dosed - the dose may need "
                "a scenario: a variant, an MTX background, or a different readout time)")

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
