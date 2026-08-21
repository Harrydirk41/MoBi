r"""Deterministic smoke test for the structural edit (a designed pathway inhibitor).

NOT the agent loop - this proves the model editing works and the knob is live before
wrapping an agent around it. Model-agnostic (--model loads projects/<name>/tasks.json).
For a chosen target driver, it sweeps the drug efficacy and prints the resulting
first-line response. If the response moves with efficacy and there are no MATLAB
errors, the edit works.

    python -m examples.run_qsp_designed_drug --model ra ^
        --sbproj "..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj" ^
        --vpop   "..\RA-QSP-Model\Vpop1.xlsx" --target F_IL6 --limit 40
"""

from __future__ import annotations

import argparse
import sys
import traceback

from pkpd_agent.engines.simbiology import SimBiologyEngine
from pkpd_agent.engines import qsp_config, qsp_tasks


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="ra")
    ap.add_argument("--sbproj", required=True)
    ap.add_argument("--vpop", required=True)
    ap.add_argument("--target", default=None, help="disease-driver parameter to suppress")
    ap.add_argument("--efficacies", default="0,0.5,0.9,0.99")
    ap.add_argument("--background", default=None)
    ap.add_argument("--start-day", type=float, default=None)
    ap.add_argument("--readout-day", type=float, default=None)
    ap.add_argument("--limit", type=int, default=40)
    args = ap.parse_args()

    tcfg = qsp_config.get(args.model)
    target = args.target or next(iter(tcfg.design_targets), "F_IL6")
    background = args.background if args.background is not None else tcfg.design_background
    start_day = args.start_day if args.start_day is not None \
        else tcfg.timeline.get("baseline_day", 200.0)
    readout_day = args.readout_day if args.readout_day is not None \
        else tcfg.timeline.get("first_line_readout_day", 284.0)
    dcol = (tcfg.run_columns.get("das28") or {}).get("readout", "DAS28_read")

    try:
        sys.stdout.reconfigure(line_buffering=True)
    except Exception:                                  # noqa: BLE001
        pass
    effs = [float(x) for x in args.efficacies.split(",") if x.strip()]
    roles = [r for r in tcfg.run_columns.get("first_line", {}) if r != "n"]

    sb = SimBiologyEngine()
    try:
        print("-> start engine"); sb.start(); print("<- started")
        print(f"target driver: {target}; background: {background or 'none'}; "
              f"readout day {readout_day:g}")
        header = "efficacy | " + " ".join(f"{r:>8}" for r in roles) + f" | {dcol:>10}"
        print("\n" + header)
        print("-" * len(header))
        for eff in effs:
            sb.load_project(args.sbproj)
            if eff > 0:
                sb.add_drug(target, eff, start_day)
            r = sb.run_vpop(args.vpop, dose=background, stop_time=400.0,
                            baseline_day=start_day, readout_day=readout_day,
                            limit=args.limit)
            fl = tcfg.summarize_run(r)["first_line"]
            dm = qsp_tasks._finite((r.get("columns") or {}).get(dcol, []))
            dmean = round(sum(dm) / len(dm), 2) if dm else "n/a"
            cells = " ".join(f"{str(fl.get(role)):>8}" for role in roles)
            print(f"{eff:>8.2f} | {cells} | {str(dmean):>10}")
        print("\n(if the response rises as efficacy rises, the edit + knob work)")
    except Exception as exc:                            # noqa: BLE001
        print(f"[FAIL] {type(exc).__name__}: {exc}")
        print(traceback.format_exc())
    finally:
        sb.stop()


if __name__ == "__main__":
    main()
