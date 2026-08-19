r"""Deterministic smoke test for the Stage-2 structural edit (designed anti-cytokine
drug). NOT the agent loop - this proves the model editing works and the knob is live
before we wrap an agent around it.

For a chosen target cytokine driver, it sweeps the drug efficacy and prints the
resulting first-line ACR (MTX background + the designed drug from day 200, readout
day 284). efficacy 0 = MTX alone; higher efficacy on an effective target should
raise ACR. If ACR moves with efficacy and there are no MATLAB errors, the edit
works and we can build the design-the-drug agent.

    python -m examples.run_ra_designed_drug ^
        --sbproj "..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj" ^
        --vpop   "..\RA-QSP-Model\Vpop1.xlsx" --target F_IL6 --limit 40
"""

from __future__ import annotations

import argparse
import os
import sys
import traceback

from pkpd_agent.engines.simbiology import SimBiologyEngine
from pkpd_agent.engines import osp_ra_trial


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sbproj", required=True)
    ap.add_argument("--vpop", required=True)
    ap.add_argument("--target", default="F_IL6",
                    help="disease-driver parameter the drug suppresses (e.g. F_IL6, "
                         "F_TNFa, F_IL17)")
    ap.add_argument("--efficacies", default="0,0.5,0.9,0.99",
                    help="comma-separated efficacy sweep")
    ap.add_argument("--background", default="MTX_15mg_Q1W_SC_t200",
                    help="background therapy dose (blank for drug monotherapy)")
    ap.add_argument("--start-day", type=float, default=200.0)
    ap.add_argument("--readout-day", type=float, default=284.0)
    ap.add_argument("--limit", type=int, default=40)
    args = ap.parse_args()

    try:
        sys.stdout.reconfigure(line_buffering=True)
    except Exception:                                  # noqa: BLE001
        pass
    effs = [float(x) for x in args.efficacies.split(",") if x.strip()]

    sb = SimBiologyEngine()
    try:
        print("-> start engine"); sb.start(); print("<- started")
        print(f"target driver: {args.target}; background: {args.background or 'none'}; "
              f"readout day {args.readout_day:g}")
        print(f"\n{'efficacy':>9} | {'ACR20':>6} {'ACR50':>6} {'ACR70':>6} "
              f"{'remission':>9} | {'DAS28_read':>10}")
        print("-" * 60)
        for eff in effs:
            # reset the model, then add the designed drug at this efficacy
            sb.load_project(args.sbproj)
            if eff > 0:
                sb.add_drug(args.target, eff, args.start_day)
            r = sb.run_vpop(args.vpop, dose=args.background, stop_time=400.0,
                            baseline_day=args.start_day, readout_day=args.readout_day,
                            limit=args.limit)
            fl = osp_ra_trial.summarize_run(r)["first_line"]
            das = (r.get("columns") or {}).get("DAS28_read", [])
            dm = osp_ra_trial._finite(das)
            dmean = round(sum(dm) / len(dm), 2) if dm else "n/a"
            print(f"{eff:>9.2f} | {fl.get('ACR20'):>6} {fl.get('ACR50'):>6} "
                  f"{fl.get('ACR70'):>6} {fl.get('remission'):>9} | {dmean:>10}")
        print("\n(if ACR rises as efficacy rises, the structural edit + knob work)")
    except Exception as exc:                            # noqa: BLE001
        print(f"[FAIL] {type(exc).__name__}: {exc}")
        print(traceback.format_exc())
    finally:
        sb.stop()


if __name__ == "__main__":
    main()
