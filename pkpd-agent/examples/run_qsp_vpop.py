r"""Run the QSP virtual population and report the MODEL's OWN clinical readouts.

Model-agnostic (--model loads projects/<name>/tasks.json). The model encodes the
trial as events, so the response is read from the model, never recomputed. This is
the "pipeline works" smoke test before wiring the agent loop on top.

    python -m examples.run_qsp_vpop --model ra ^
        --sbproj "C:\...\Vantage RA QSP Model v1.0.sbproj" ^
        --vpop "C:\...\Vpop1.xlsx" --dose MTX_15mg_Q1W_SC_t200 --limit 20
"""

from __future__ import annotations

import argparse
import os
import sys
import traceback

from pkpd_agent.engines.simbiology import SimBiologyEngine
from pkpd_agent.engines import qsp_config, qsp_tasks

_LOG = None


def log(msg: str) -> None:
    print(msg, flush=True)
    if _LOG is not None:
        try:
            _LOG.write(msg + "\n"); _LOG.flush()
        except Exception:                              # noqa: BLE001
            pass


def _summ(xs):
    xs = qsp_tasks._finite(xs)
    if not xs:
        return "n/a"
    m = sum(xs) / len(xs)
    sd = (sum((x - m) ** 2 for x in xs) / len(xs)) ** 0.5
    return f"n={len(xs)} mean={m:.3f} sd={sd:.3f} min={min(xs):.3f} max={max(xs):.3f}"


def main() -> None:
    global _LOG
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="vantage_ra")
    ap.add_argument("--sbproj", required=True)
    ap.add_argument("--vpop", required=True)
    ap.add_argument("--dose", default="", help="dose name(s) joined by ';' (default: baseline)")
    ap.add_argument("--stop-time", type=float, default=700.0)
    ap.add_argument("--baseline-day", type=float, default=None)
    ap.add_argument("--readout-day", type=float, default=None)
    ap.add_argument("--limit", type=int, default=20)
    args = ap.parse_args()

    tcfg = qsp_config.get(args.model)
    baseline_day = args.baseline_day if args.baseline_day is not None \
        else tcfg.timeline.get("baseline_day", 200.0)
    readout_day = args.readout_day if args.readout_day is not None \
        else tcfg.timeline.get("first_line_readout_day", 284.0)

    try:
        sys.stdout.reconfigure(line_buffering=True)
    except Exception:                                  # noqa: BLE001
        pass
    _LOG = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "run_qsp_vpop.log"), "w", encoding="utf-8")

    sb = SimBiologyEngine()
    try:
        log("-> start engine"); sb.start(); log("<- engine started")
        log(f"-> load {args.sbproj}"); sb.load_project(args.sbproj); log("<- loaded")

        label = args.dose or "PLACEBO (no drug)"
        log(f"trial: {label}; stop day {args.stop_time:g}; baseline day "
            f"{baseline_day:g}, first-line readout day {readout_day:g}")

        log(f"-> run vpop, limit={args.limit} (each patient = one sim)...")
        r = sb.run_vpop(args.vpop, dose=args.dose, stop_time=args.stop_time,
                        baseline_day=baseline_day, readout_day=readout_day,
                        limit=args.limit)
        ml = (r.get("matlab_log") or "").strip()
        if ml:
            log("   [MATLAB] " + ml.replace("\n", "\n   [MATLAB] "))

        cols = r.get("columns") or {}
        dcols = tcfg.run_columns.get("severity") or {}
        log(f"\n{tcfg.severity_readout} baseline: {_summ(cols.get(dcols.get('baseline',''), []))}")
        log(f"{tcfg.severity_readout} first-line readout: "
            f"{_summ(cols.get(dcols.get('readout',''), []))}")

        summary = tcfg.summarize_run(r)
        log(f"\nFIRST-LINE response under {label} (model-computed flags):")
        for role, val in summary["first_line"].items():
            if role == "n":
                continue
            log(f"   {role:<12}: " + (f"{val}%" if val is not None else "n/a"))

        sl = summary["second_line"]
        log(f"\nSECOND-LINE (subgroup n={sl.get('n_subgroup')}):")
        if not sl.get("n_subgroup"):
            log("   (no subgroup patients flagged)")
        else:
            for role, val in sl.items():
                if role == "n_subgroup":
                    continue
                log(f"   {role:<12}: " + (f"{val}%" if val is not None else "n/a"))

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
