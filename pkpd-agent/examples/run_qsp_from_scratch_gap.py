r"""Controlled-failure experiment: quantify the CLINICAL cost of the one from-scratch defect we
can inject cleanly - the parameters L1 proved are UNIDENTIFIABLE without dose-response data
(the half-effect concentrations and Hill slopes). A from-scratch modeller, lacking that data,
would set them to priors, not the calibrated values. So set them to priors here, run the same
held-out per-arm clinical readout as the real model, and measure how far the prediction moves.

    python -m examples.run_qsp_from_scratch_gap --model ra ^
        --sbproj "..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj" ^
        --vpop "..\RA-QSP-Model\Vpop1.xlsx" --perturb 3 --seed 1

This is NOT a full from-scratch model (reconstructing all 59 species would itself require the
answer model). It isolates ONE defect - unpinnable shape parameters -> priors - on the real
model, so the clinical degradation is attributable. The topology defect (missed load-bearing
edges) is not injected here because those edges have no clean single knob (benchmark B).
Needs the MATLAB engine.
"""

from __future__ import annotations

import argparse
import os
import random
import re

from pkpd_agent.engines.simbiology import SimBiologyEngine
from pkpd_agent.engines import qsp_config


def _arm_rates(sb, cfg, vpop, first_line, readout_day, limit):
    out = {}
    for d, spec in cfg.drugs.items():
        doses = (spec or {}).get("doses") or []
        if not doses:
            continue
        r = sb.run_vpop(vpop, dose=doses[0], stop_time=readout_day + 2,
                        baseline_day=cfg.timeline.get("baseline_day", 200.0),
                        readout_day=readout_day, limit=limit,
                        states=cfg.readout_states or None)
        cols = r.get("columns") or {}
        rates = {}
        for role, colname in first_line.items():
            col = [v for v in cols.get(colname, []) if isinstance(v, (int, float)) and v == v]
            rates[role] = round(100.0 * sum(1 for v in col if v >= 0.5) / len(col), 1) \
                if col else None
        out[d] = rates
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="ra")
    ap.add_argument("--sbproj", required=True)
    ap.add_argument("--vpop", required=True)
    ap.add_argument("--pattern", default=r"(?i)HalfEffectConc|Slope_",
                    help="params to treat as unidentifiable (set to priors)")
    ap.add_argument("--perturb", type=float, default=3.0,
                    help="prior spread: each shape param x a random factor in [1/p, p]")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    cfg = qsp_config.get(args.model)
    first_line = cfg.run_columns.get("first_line") or {}
    readout_day = cfg.timeline.get("first_line_readout_day", 284.0)

    sb = SimBiologyEngine()
    try:
        print("== starting MATLAB =="); sb.start()
        print(f"== loading {os.path.basename(args.sbproj)} =="); sb.load_project(args.sbproj)
        rgx = re.compile(args.pattern)
        shape = [p["name"] for p in sb.list_parameters().get("parameters", [])
                 if rgx.search(p["name"])]
        print(f"unidentifiable shape parameters (L1): {len(shape)} "
              f"(e.g. {shape[:3]})")

        print("\n== baseline: real (calibrated) model, held-out per-arm ACR ==", flush=True)
        base = _arm_rates(sb, cfg, args.vpop, first_line, readout_day, args.limit)
        for d, r in base.items():
            print(f"  {d:<5} " + " ".join(f"{role}={r.get(role)}" for role in first_line))

        # inject the defect: set every unidentifiable shape param to a PRIOR guess
        rng = random.Random(args.seed)
        saved = {}
        for p in sb.list_parameters().get("parameters", []):
            if rgx.search(p["name"]):
                old = float(p["value"])
                saved[p["name"]] = old
                sb.set_parameter(p["name"], old * args.perturb ** rng.uniform(-1, 1))
        print(f"\n== degraded: {len(saved)} shape params set to priors (x[1/{args.perturb:g}, "
              f"{args.perturb:g}]), same Vpop, same readout ==", flush=True)
        deg = _arm_rates(sb, cfg, args.vpop, first_line, readout_day, args.limit)
        for name, old in saved.items():
            sb.set_parameter(name, old)                # restore

        print(f"\n== held-out clinical prediction: real vs from-scratch-defect ==")
        print(f"  {'arm':<5} {'metric':<7} {'real':>6} {'degraded':>9} {'shift':>7}")
        shifts = []
        for d in base:
            for role in first_line:
                a, b = base[d].get(role), deg[d].get(role)
                if a is not None and b is not None:
                    shifts.append(abs(b - a))
                    print(f"  {d:<5} {role:<7} {a:>6} {b:>9} {b-a:>+7.1f}")
        if shifts:
            print(f"\n  mean absolute shift in ACR rate: {sum(shifts)/len(shifts):.1f} points")
        print("\n  -> this is the clinical cost of ONE from-scratch defect (parameters that "
              "cannot be\n     pinned without dose-response data). A real from-scratch model "
              "also carries the\n     topology defect (missed load-bearing edges), so this is a "
              "LOWER bound on the gap.")
    finally:
        sb.stop()


if __name__ == "__main__":
    main()
