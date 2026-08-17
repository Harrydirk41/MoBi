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


_INF = (float("inf"), float("-inf"))


def _finite(xs):
    return [x for x in xs
            if isinstance(x, (int, float)) and x == x and x not in _INF]


def _summ(xs):
    xs = _finite(xs)
    if not xs:
        return "n/a"
    m = sum(xs) / len(xs)
    sd = (sum((x - m) ** 2 for x in xs) / len(xs)) ** 0.5
    return (f"n={len(xs)} mean={m:.3f} sd={sd:.3f} "
            f"min={min(xs):.3f} max={max(xs):.3f}")


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
                    help="absolute time (days) to read the endpoints (0 = sim end)")
    ap.add_argument("--week", type=float, default=None,
                    help="clinical week POST treatment-start (overrides --readout-time); "
                         "readout day = tx-start + week*7")
    ap.add_argument("--tx-start", type=float, default=200.0,
                    help="day treatment starts (the _t200 doses begin at day 200)")
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

        lim = args.limit
        bl_day = args.tx_start
        week = args.week if args.week is not None else 12.0
        rd_day = bl_day + week * 7.0
        label = args.dose or "PLACEBO (no drug)"
        log(f"trial: {label}; baseline day {bl_day:g} (treatment start), "
            f"readout week {week:g} = day {rd_day:g}")

        log(f"-> run vpop, limit={lim} (each patient = one sim)...")
        r = sb.run_vpop(args.vpop, dose=args.dose, baseline_day=bl_day,
                        readout_day=rd_day, limit=lim)
        ml = (r.get("matlab_log") or "").strip()
        if ml:
            log("   [MATLAB] " + ml.replace("\n", "\n   [MATLAB] "))

        base = _col(r, "DAS28_base")
        read = _col(r, "DAS28_read")
        log(f"<- {len(base)} patients")
        log(f"\nDAS28 at baseline (day {bl_day:g}): {_summ(base)}")
        log(f"DAS28 at readout  (day {rd_day:g}): {_summ(read)}")

        # ACR = % DAS28 improvement from each patient's OWN baseline (model's def)
        acrp = [100.0 * (b - t) / b for b, t in zip(base, read) if b and b > 0]
        n = len(acrp)
        log(f"\nresponse under {label}, n={n} "
            f"(ACR = % DAS28 improvement from own baseline):")
        if n:
            for thr, name in ((20, "ACR20"), (50, "ACR50"), (70, "ACR70")):
                log(f"   {name}: {100.0 * sum(1 for x in acrp if x >= thr) / n:.1f}%")
            log(f"   mean DAS28 improvement: {sum(acrp) / n:.1f}%")
            rem = 100.0 * sum(1 for x in read if x <= 2.6) / len(read)
            low = 100.0 * sum(1 for x in read if x <= 3.2) / len(read)
            log(f"   DAS28 remission (<=2.6): {rem:.1f}%   low activity (<=3.2): {low:.1f}%")

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
