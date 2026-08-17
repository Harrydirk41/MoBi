r"""Run the RA QSP virtual population and report the MODEL's OWN clinical readouts.

The Vantage RA model encodes the whole trial as events, so the response is read
from the model, never recomputed:
  * DAS28_BL captured at day 199, ACR20/ACR50/ACR70/Remission set at day 284
    (first-line, week 12);
  * for MTX-inadequate responders (MTX_NonResp==1), MTX_NonResp_TCZ_ACR20/50/70/Rem
    set at day 600 (the second-line TCZ readout - the paper's held-out validation).

Reproducing the paper's baseline DAS28 spread and the MTX/ADA/TCZ response RATES
is the "pipeline works" milestone before wiring the agent loop on top.

    # first-line MTX, quick check on 20 patients:
    python -m examples.run_ra_vpop --sbproj "C:\...\Vantage RA QSP Model v1.0.sbproj" ^
        --vpop "C:\...\Vpop1.xlsx" --dose MTX_15mg_Q1W_SC_t200 --limit 20

    # full population, MTX first-line + TCZ (captures the second-line TCZ-in-MTX-IR
    # flagship in the same run):
    python -m examples.run_ra_vpop --sbproj "...sbproj" --vpop "...Vpop1.xlsx" ^
        --dose "MTX_15mg_Q1W_SC_t200;TCZ8mgkg_Q4W_IV_t200" --stop-time 700 --limit 300
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


def _rate(flags):
    """Percent of finite flags that equal 1 (a model-set responder flag)."""
    xs = _finite(flags)
    if not xs:
        return None, 0
    return 100.0 * sum(1 for x in xs if x >= 0.5) / len(xs), len(xs)


def main() -> None:
    global _LOG
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sbproj", required=True)
    ap.add_argument("--vpop", required=True, help="Vpop1.xlsx / Vpop2.xlsx")
    ap.add_argument("--dose", default="",
                    help="dose name, or several joined by ';' (default: baseline, no drug). "
                         "A 'name@START' suffix retimes that dose's StartTime (days), e.g. "
                         "'MTX_15mg_Q1W_SC_t200;TCZ8mgkg_Q4W_IV_t200@285' gives MTX from day "
                         "200 and switches TCZ in at day 285 (the sequential MTX-then-TCZ trial)")
    ap.add_argument("--stop-time", type=float, default=700.0,
                    help="force sim end (days); needs >=601 to capture the second-line "
                         "TCZ-in-MTX-IR readout (day 600)")
    ap.add_argument("--baseline-day", type=float, default=200.0,
                    help="day to report DAS28 baseline (treatment start; flags use day 199)")
    ap.add_argument("--readout-day", type=float, default=284.0,
                    help="day to report DAS28 first-line readout (week 12)")
    ap.add_argument("--limit", type=int, default=20,
                    help="run only the first N patients (default 20; use 300 for the full Vpop)")
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

        label = args.dose or "PLACEBO (no drug)"
        log(f"trial: {label}; stop day {args.stop_time:g}; "
            f"DAS28 baseline day {args.baseline_day:g}, first-line readout day {args.readout_day:g}")

        log(f"-> run vpop, limit={args.limit} (each patient = one sim)...")
        r = sb.run_vpop(args.vpop, dose=args.dose, stop_time=args.stop_time,
                        baseline_day=args.baseline_day, readout_day=args.readout_day,
                        limit=args.limit)
        ml = (r.get("matlab_log") or "").strip()
        if ml:
            log("   [MATLAB] " + ml.replace("\n", "\n   [MATLAB] "))

        base = _col(r, "DAS28_base")
        read = _col(r, "DAS28_read")
        n = len(base)
        log(f"<- {n} patients")
        log(f"\nDAS28 baseline (day {args.baseline_day:g}): {_summ(base)}")
        log(f"DAS28 first-line readout (day {args.readout_day:g}): {_summ(read)}")

        # ---- first-line response: read the model's OWN ACR/Remission flags ----
        log(f"\nFIRST-LINE response under {label} (model-computed flags, day 284):")
        for key, name in (("ACR20", "ACR20"), ("ACR50", "ACR50"),
                          ("ACR70", "ACR70"), ("Rem", "Remission (DAS28<=2.6)")):
            rate, k = _rate(_col(r, key))
            log(f"   {name:<24}: " + (f"{rate:.1f}%  (of {k})" if rate is not None else "n/a"))

        # ---- second-line: TCZ in MTX-inadequate responders (the flagship) ----
        # use the RAW (row-aligned) columns so MTX_NonResp pairs with each TCZ flag
        cols = r.get("columns") or {}
        nonresp_raw = cols.get("MTX_NonResp", [])
        def _isnr(v):
            return isinstance(v, (int, float)) and v == v and v >= 0.5
        n_ir = sum(1 for v in nonresp_raw if _isnr(v))
        n_nr_finite = sum(1 for v in nonresp_raw
                          if isinstance(v, (int, float)) and v == v)
        log(f"\nSECOND-LINE (flagship) TCZ in MTX-non-responders "
            f"(MTX_NonResp==1: {n_ir}/{n_nr_finite} patients):")
        if n_ir == 0:
            log("   (no MTX non-responders flagged - either no TCZ dose given, "
                "sim ended before day 600, or the MTX arm resolved everyone)")
        else:
            def _among_ir(key):
                col = cols.get(key, [])
                vals = [v for v, nr in zip(col, nonresp_raw)
                        if isinstance(v, (int, float)) and v == v and _isnr(nr)]
                if not vals:
                    return None, 0
                return 100.0 * sum(1 for v in vals if v >= 0.5) / len(vals), len(vals)
            for key, name in (("TCZ_ACR20", "ACR20"), ("TCZ_ACR50", "ACR50"),
                              ("TCZ_ACR70", "ACR70"), ("TCZ_Rem", "Remission")):
                rate, k = _among_ir(key)
                log(f"   TCZ {name:<10}: " + (f"{rate:.1f}%  (of {k} MTX-IR)"
                                              if rate is not None else "n/a"))

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
