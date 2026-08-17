r"""Run the RA QSP virtual population — the foundation for the virtual-trial agent.

Loads the model + a virtual-population .xlsx (Vpop1/Vpop2 from the Vantage repo),
applies each patient's parameter set, simulates (optionally under a drug regimen),
and reports the population clinical readouts (DAS28-CRP distribution, ACR20/50/70
and remission RATES). Reproducing the paper's baseline DAS28 spread and the
MTX/ADA/TCZ response rates is the "pipeline works" milestone before the agent loop.

    # quick check on 5 patients, no drug (baseline disease):
    python -m examples.run_ra_vpop --sbproj "C:\...\Vantage RA QSP Model v1.0.sbproj" ^
        --vpop "C:\...\Vpop1.xlsx" --limit 5

    # full population under a TCZ regimen:
    python -m examples.run_ra_vpop --sbproj "...sbproj" --vpop "...Vpop1.xlsx" ^
        --dose TCZ8mgkg_Q4W_IV_t200
"""

from __future__ import annotations

import argparse
import os
import statistics
import sys
import traceback

from pkpd_agent.engines.simbiology import SimBiologyEngine

_LOG = None


def log(msg: str) -> None:
    print(msg, flush=True)
    if _LOG is not None:
        try:
            _LOG.write(msg + "\n"); _LOG.flush()
        except Exception:                              # noqa: BLE001
            pass


def _col(res, name):
    return [v for v in (res.get("columns") or {}).get(name, []) if isinstance(v, (int, float))]


def _summ(xs):
    xs = [x for x in xs if isinstance(x, (int, float)) and x == x]
    if not xs:
        return "n/a"
    return (f"n={len(xs)} mean={statistics.mean(xs):.3f} "
            f"sd={statistics.pstdev(xs):.3f} min={min(xs):.3f} max={max(xs):.3f}")


def _rate(xs):
    xs = [x for x in xs if isinstance(x, (int, float)) and x == x]
    if not xs:
        return "n/a"
    return f"{100.0 * sum(1 for x in xs if x >= 0.5) / len(xs):.1f}% (of {len(xs)})"


def main() -> None:
    global _LOG
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sbproj", required=True)
    ap.add_argument("--vpop", required=True, help="Vpop1.xlsx / Vpop2.xlsx")
    ap.add_argument("--dose", default="", help="dose regimen name (default: baseline, no drug)")
    ap.add_argument("--readout-time", type=float, default=0.0,
                    help="time at which to read the clinical endpoints (0 = sim end)")
    ap.add_argument("--limit", type=int, default=5,
                    help="run only the first N patients (default 5; use 300 for the full Vpop)")
    args = ap.parse_args()

    try:
        sys.stdout.reconfigure(line_buffering=True)
    except Exception:                                  # noqa: BLE001
        pass
    _LOG = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "run_ra_vpop.log"), "w", encoding="utf-8")

    sb = SimBiologyEngine()
    try:
        log("-> start engine"); sb.start(); log("<- engine started")
        log(f"-> load {args.sbproj}"); sb.load_project(args.sbproj); log("<- loaded")

        label = args.dose or "BASELINE (no drug)"
        log(f"-> run vpop '{os.path.basename(args.vpop)}' under {label}, "
            f"limit={args.limit} (each patient is one simulation)...")
        res = sb.run_vpop(args.vpop, dose=args.dose,
                          readout_time=args.readout_time, limit=args.limit)
        mlog = (res.get("matlab_log") or "").strip()
        if mlog:
            log("   [MATLAB] " + mlog.replace("\n", "\n   [MATLAB] "))

        pts = _col(res, "patient")
        log(f"<- {len(pts)} patients returned\n")

        log(f"population clinical readouts under {label}:")
        log(f"   DAS28_CRP : {_summ(_col(res, 'DAS28_CRP'))}")
        log(f"   DAS28_BL  : {_summ(_col(res, 'DAS28_BL'))}")
        log(f"   ACR_Perc  : {_summ(_col(res, 'ACR_Perc'))}   (% DAS28 improvement)")

        # RESPONSE from the model's own definition: ACR_Perc thresholds + DAS28
        # remission. This is robust to the event-set ACR flags reading 0.
        acrp = _col(res, "ACR_Perc")
        das = _col(res, "DAS28_CRP")
        n = len(acrp)
        if n:
            log("\n   response rates (thresholding ACR_Perc = % DAS28 improvement):")
            for thr, name in ((20, "ACR20"), (50, "ACR50"), (70, "ACR70")):
                rate = 100.0 * sum(1 for x in acrp if x >= thr) / n
                log(f"      {name}: {rate:.1f}%")
            rem = 100.0 * sum(1 for x in das if x <= 2.6) / len(das) if das else 0.0
            low = 100.0 * sum(1 for x in das if x <= 3.2) / len(das) if das else 0.0
            log(f"      DAS28 remission (<=2.6): {rem:.1f}%    low activity (<=3.2): {low:.1f}%")
        log("\n   (model's event-set flags, for comparison - expected unreliable:)")
        for r in ("ACR20", "ACR50", "ACR70", "Remission", "Response"):
            log(f"      {r:9} state rate: {_rate(_col(res, r))}")

        log("\n=== VPOP RUN END (reached the end cleanly) ===")
    except Exception as exc:                            # noqa: BLE001
        log(f"[FAIL] {type(exc).__name__}: {exc}")
        log(traceback.format_exc())
    finally:
        sb.stop()
        if _LOG is not None:
            _LOG.close()


if __name__ == "__main__":
    main()
