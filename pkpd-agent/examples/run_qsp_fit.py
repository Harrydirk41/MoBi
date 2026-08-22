r"""Smoke test for sb_fit.m (SimBiology native sbiofit) - no API, no agent.

Self-contained: it SIMULATES the model once to produce a synthetic data curve for a
severity state, writes it as a data table, then asks sbiofit to estimate one parameter
back against that curve. If sbiofit runs and returns an estimate near the model's own
value (residual ~0 at the true value), the native-calibration plumbing works. This
validates the .m mechanics, not the science.

    python -m examples.run_qsp_fit --model ra ^
        --sbproj "..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj" ^
        --param F_TNFa --state DAS28_CRP --method lsqnonlin
"""

from __future__ import annotations

import argparse
import os
import tempfile

from pkpd_agent.engines.simbiology import SimBiologyEngine
from pkpd_agent.engines import qsp_config


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="ra")
    ap.add_argument("--sbproj", required=True)
    ap.add_argument("--param", default=None, help="parameter to estimate (default: a driver)")
    ap.add_argument("--state", default="DAS28_CRP", help="observed severity state")
    ap.add_argument("--method", default="lsqnonlin",
                    help="sbiofit optimizer: lsqnonlin/fmincon/particleswarm/ga/...")
    ap.add_argument("--stop-time", type=float, default=200.0)
    args = ap.parse_args()

    cfg = qsp_config.get(args.model)
    param = args.param or next(iter(cfg.vpop_drivers), None)
    if not param:
        print("no parameter to fit (pass --param).")
        return
    span = (cfg.vpop_drivers.get(param) or {}).get("span") or [1e-3, 1e3]

    sb = SimBiologyEngine()
    try:
        print("== starting MATLAB engine =="); sb.start()
        print(f"== loading {os.path.basename(args.sbproj)} =="); sb.load_project(args.sbproj)

        # 1. generate a synthetic data curve for the state (baseline, no drug)
        print(f"== simulating to produce a '{args.state}' data curve ==", flush=True)
        sim = sb.simulate(stop_time=args.stop_time)
        t = sim.get("time", [])
        y = (sim.get("columns") or {}).get(args.state, [])
        if not y:
            print(f"state '{args.state}' not in simulation output - try another --state. "
                  f"Available: {list((sim.get('columns') or {}).keys())[:20]}")
            return
        # sample ~8 evenly spaced points
        idx = [int(k) for k in
               (i * (len(t) - 1) // 7 for i in range(8))] if len(t) > 8 else list(range(len(t)))
        data_path = os.path.join(tempfile.gettempdir(), "sb_fit_data.csv")
        with open(data_path, "w", encoding="utf-8") as fh:
            fh.write("Time,obs\n")
            for i in idx:
                fh.write(f"{t[i]},{y[i]}\n")
        print(f"   wrote {len(idx)} data points to {data_path}")

        # 2. fit the parameter back with native sbiofit
        param_spec = f"{param},{span[0]:g},{span[1]:g},log"
        response_map = f"{args.state} = obs"
        print(f"== sb_fit: estimate {param} in [{span[0]:g},{span[1]:g}] via {args.method} ==",
              flush=True)
        r = sb.fit_native(param_spec, data_path, response_map, method=args.method)
        ml = (r.get("matlab_log") or "").strip()
        if ml:
            print("   [MATLAB] " + ml.replace("\n", "\n   [MATLAB] "))
        cols = r.get("columns") or {}
        names = cols.get("parameter", [])
        ests = cols.get("estimate", [])
        print("\n== estimates ==")
        for nm, e in zip(names, ests):
            print(f"  {nm} = {e}")
        if not names:
            print("  (no estimate returned - check the MATLAB log above)")
    finally:
        sb.stop()


if __name__ == "__main__":
    main()
