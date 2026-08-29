r"""Measure each therapy arm's first-line response rate from the SHIPPED virtual
population - the same way the MTX target (~33.7% ACR20 at day 284) was obtained.

The config's MTX target is labelled 'day284_full_pop_n300': it was measured by running
the shipped Vpop under MTX, not taken from a paper table (the paper only shows Figure 5).
So the ADA/TCZ multi-arm targets can be measured the same way - and reading every arm at
the SAME day makes them self-consistent, sidestepping the per-arm-timepoint problem.

    python -m examples.run_qsp_arm_targets --model ra ^
        --sbproj "..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj" ^
        --vpop "..\RA-QSP-Model\Vpop1.xlsx" --drugs MTX,ADA,TCZ

Prints each arm's ACR20/50/70 rate at the first-line readout day. Feed these into
projects/<name>/tasks.json:vpop_anchors to match multiple arms. Needs the MATLAB engine.
"""

from __future__ import annotations

import argparse
import os

from pkpd_agent.engines.simbiology import SimBiologyEngine
from pkpd_agent.engines import qsp_config


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="ra")
    ap.add_argument("--sbproj", required=True)
    ap.add_argument("--vpop", required=True, help="the shipped Vpop .xlsx")
    ap.add_argument("--drugs", default="MTX,ADA,TCZ", help="comma-separated drug keys")
    ap.add_argument("--limit", type=int, default=0, help="subsample N patients (0 = all)")
    args = ap.parse_args()

    cfg = qsp_config.get(args.model)
    first_line = cfg.run_columns.get("first_line") or {}
    readout_day = cfg.timeline.get("first_line_readout_day", 284.0)
    want = [d.strip() for d in args.drugs.split(",") if d.strip()]

    sb = SimBiologyEngine()
    try:
        print("== starting MATLAB engine =="); sb.start()
        print(f"== loading {os.path.basename(args.sbproj)} =="); sb.load_project(args.sbproj)
        print(f"\n{'drug':<6} {'dose':<28} " +
              " ".join(f"{r:>7}" for r in first_line))
        for d in want:
            doses = (cfg.drugs.get(d) or {}).get("doses") or []
            if not doses:
                print(f"{d:<6} (no dose in config)"); continue
            dose = doses[0]
            r = sb.run_vpop(args.vpop, dose=dose, stop_time=readout_day + 2,
                            baseline_day=cfg.timeline.get("baseline_day", 200.0),
                            readout_day=readout_day, limit=args.limit,
                            states=cfg.readout_states or None)
            cols = r.get("columns") or {}
            rates = {}
            for role, colname in first_line.items():
                col = [v for v in cols.get(colname, []) if isinstance(v, (int, float)) and v == v]
                rates[role] = round(100.0 * sum(1 for v in col if v >= 0.5) / len(col), 1) \
                    if col else None
            print(f"{d:<6} {dose:<28} " +
                  " ".join(f"{str(rates.get(r)):>7}" for r in first_line))
        print(f"\n(rates at day {readout_day:g}; use these as vpop_anchors.rate_targets to "
              "match multiple arms)")
    finally:
        sb.stop()


if __name__ == "__main__":
    main()
