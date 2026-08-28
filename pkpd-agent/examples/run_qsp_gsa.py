r"""Smoke test for sb_gsa.m (native Sobol global sensitivity) - no API, no agent.

Ranks the disease-driver parameters by how much they drive the clinical readout, computed
by MATLAB's own variance-based sensitivity analysis (sbiosobol) -- NOT read from a stored
list. This is the numerical half of "which parameters to vary": an agent would then apply
biological reasoning on top of this ranking. To show it is genuinely derived, the runner
also reports the overlap with the project's stored gsa_top (the paper's list).

    python -m examples.run_qsp_gsa --model ra ^
        --sbproj "..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj" ^
        --observable DAS28_CRP --samples 500
"""

from __future__ import annotations

import argparse
import os

from pkpd_agent.engines.simbiology import SimBiologyEngine
from pkpd_agent.engines import qsp_config, qsp_model


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="ra")
    ap.add_argument("--sbproj", required=True)
    ap.add_argument("--observable", default="DAS28_CRP", help="model readout to rank on")
    ap.add_argument("--day", type=float, default=None, help="readout day (default: first-line)")
    ap.add_argument("--samples", type=int, default=500, help="Sobol sample size")
    args = ap.parse_args()

    cfg = qsp_config.get(args.model)
    drivers = cfg.vpop_drivers
    if not drivers:
        print(f"project '{args.model}' declares no vpop_drivers to screen.")
        return
    param_spec = ";".join(f"{n},{p['span'][0]:g},{p['span'][1]:g}"
                          for n, p in drivers.items())
    day = args.day if args.day is not None else cfg.timeline.get("first_line_readout_day", 284.0)

    sb = SimBiologyEngine()
    try:
        print("== starting MATLAB engine =="); sb.start()
        print(f"== loading {os.path.basename(args.sbproj)} =="); sb.load_project(args.sbproj)
        print(f"== sb_gsa: Sobol over {len(drivers)} drivers, obs {args.observable} "
              f"@ day {day:g}, {args.samples} samples ==", flush=True)
        r = sb.gsa(param_spec, args.observable, day, n_samples=args.samples)
        ml = (r.get("matlab_log") or "").strip()
        if ml:
            print("   [MATLAB] " + ml.replace("\n", "\n   [MATLAB] "))
        cols = r.get("columns") or {}
        names = cols.get("parameter", [])
        first = cols.get("first_order", [])
        total = cols.get("total_order", [])
        if not names:
            print("  no sensitivities returned - check the MATLAB log above.")
            return
        print("\n== parameter ranking (by total-order Sobol index) ==")
        print(f"  {'parameter':<26} {'first':>8} {'total':>8}")
        for nm, f, t in zip(names, first, total):
            print(f"  {nm:<26} {f if f is None else round(f,4):>8} "
                  f"{t if t is None else round(t,4):>8}")

        # show it is DERIVED, not the stored list: overlap with the paper's gsa_top
        try:
            stored = set(qsp_model.get_spec(args.model).gsa_top or [])
        except Exception:
            stored = set()
        if stored:
            hit = [nm for nm in names if nm in stored]
            print(f"\n  computed independently; {len(hit)}/{len(names)} of these also "
                  f"appear in the project's stored gsa_top ({len(stored)} params): {hit}")
    finally:
        sb.stop()


if __name__ == "__main__":
    main()
