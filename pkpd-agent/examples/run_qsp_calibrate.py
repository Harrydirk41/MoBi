r"""Stage-0 calibration, a scoped honest slice: fit a cell's proliferation rate to its
steady-state target - exactly the step the paper describes ("Rate of proliferation ... were
determined by fitting to obtain the observed steady state number of the cell type").

Two parts, both against the REAL model + the paper's own Supplementary-Data-1 targets:

  1. GAP CHECK    - simulate the shipped (already-calibrated) model to pre-drug steady state,
                    read the cell's value, compare to the MOESM1 target. Shows whether the
                    model reproduces its own target, and surfaces any unit mismatch honestly.
  2. RECOVER DEMO - perturb the cell's baseline proliferation rate away from its shipped
                    value, then use SimBiology's native fit (sbiofit) to recover it from the
                    steady-state target. This is the agent mechanically redoing a slice of
                    Stage-0 calibration: knob + target -> fitted value, compared to the truth.

    python -m examples.run_qsp_calibrate --model ra ^
        --sbproj "..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj" ^
        --targets projects\vantage_ra\data\steady_state_targets.json ^
        --cell FLS --perturb 3.0

The proliferation parameter is DISCOVERED from the loaded model by name pattern (kp/kg/prolif
+ cell); pass --param to override. Needs the MATLAB engine (+ Optimization Toolbox for the fit).
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import tempfile

from pkpd_agent.engines.simbiology import SimBiologyEngine


def _discover_prolif(params: list, cell: str, override: str) -> list:
    if override:
        return [override]
    pats = [re.compile(rf"(?i)^k[pg]_{cell}.*base", ),
            re.compile(rf"(?i){cell}.*prolif.*(base|rate)"),
            re.compile(rf"(?i)^k[pg]_{cell}\b"),
            re.compile(rf"(?i){cell}.*prolif")]
    names = [p["name"] for p in params]
    for pat in pats:
        hit = [n for n in names if pat.search(n)]
        if hit:
            return hit
    return []


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="ra")
    ap.add_argument("--sbproj", required=True)
    ap.add_argument("--targets", required=True, help="steady_state_targets.json")
    ap.add_argument("--cell", default="FLS", help="model species to calibrate (a cell)")
    ap.add_argument("--param", default="", help="override the proliferation parameter name")
    ap.add_argument("--perturb", type=float, default=3.0,
                    help="factor to perturb the rate by before refitting")
    ap.add_argument("--steady-day", type=float, default=199.0,
                    help="pre-drug day at which the disease steady state is read")
    args = ap.parse_args()

    with open(args.targets, encoding="utf-8") as fh:
        targets = json.load(fh)
    tgt = next((t for t in targets if t.get("model_species") == args.cell
                and t.get("target_model_unit") is not None), None)
    if not tgt:
        print(f"no steady-state target for model species {args.cell!r} in {args.targets}")
        return
    target_val = float(tgt["target_model_unit"])
    print(f"target: {args.cell} steady state = {target_val:g} {tgt['target_units']} "
          f"(MOESM1: {tgt['name']})")

    sb = SimBiologyEngine()
    try:
        print("== starting MATLAB =="); sb.start()
        print(f"== loading {os.path.basename(args.sbproj)} =="); sb.load_project(args.sbproj)

        params = sb.list_parameters().get("parameters", [])
        cand = _discover_prolif(params, args.cell, args.param)
        if not cand:
            print(f"could not find a proliferation parameter for {args.cell}; "
                  "pass --param NAME. Parameters containing the cell name:")
            for p in params:
                if args.cell.lower() in p["name"].lower():
                    print(f"    {p['name']} = {p.get('value')}")
            return
        if len(cand) > 1:
            print(f"multiple proliferation-rate candidates for {args.cell}: {cand}")
            print("re-run with --param to pick one.")
            return
        pname = cand[0]
        p0 = next((float(p["value"]) for p in params if p["name"] == pname), None)
        print(f"proliferation parameter: {pname} (shipped value {p0:g})")

        # ---- Phase 1: gap check ---------------------------------------------------
        sim = sb.simulate(stop_time=args.steady_day + 1.0)
        cols = sim.get("columns", {})
        if args.cell not in cols:
            print(f"model has no state named {args.cell!r}; states: "
                  f"{list(cols)[:12]}...")
            return
        s0 = float(cols[args.cell][-1])
        ratio = s0 / target_val if target_val else float("nan")
        print(f"\n== gap check ==")
        print(f"  shipped model {args.cell} steady state = {s0:g}")
        print(f"  MOESM1 target                          = {target_val:g} {tgt['target_units']}")
        print(f"  ratio model/target = {ratio:.3g}  "
              + ("(same unit, close)" if 0.2 < ratio < 5 else
                 "(differ - likely a unit/compartment scaling; recover demo is unit-safe)"))

        # ---- Phase 2: recover demo (unit-safe: fit back to the model's own S0) -----
        print(f"\n== recover demo: perturb {pname} x{args.perturb:g}, then fit it back ==")
        sb.set_parameter(pname, p0 * args.perturb)
        sp = float(sb.simulate(stop_time=args.steady_day + 1.0)["columns"][args.cell][-1])
        print(f"  after perturb: {args.cell} steady state {s0:g} -> {sp:g}")

        data = os.path.join(tempfile.gettempdir(), "calib_target.csv")
        with open(data, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["Time", args.cell])
            w.writerow([args.steady_day, s0])          # the steady-state target (model's own)
        spec = f"{pname},{p0/50:g},{p0*50:g},log"
        print(f"  fitting {pname} in [{p0/50:g}, {p0*50:g}] to reproduce {args.cell}={s0:g} ...",
              flush=True)
        res = sb.fit_native(spec, data, f"{args.cell} = {args.cell}", method="lsqnonlin")
        est = res.get("columns", {}).get("estimate", [])
        fitted = float(est[0]) if est else float("nan")
        ml = (res.get("matlab_log") or "").strip()
        if ml:
            print("   [MATLAB] " + ml.replace("\n", "\n   [MATLAB] "))

        sb.set_parameter(pname, p0)                    # restore
        print(f"\n== result ==")
        print(f"  shipped (truth) {pname} = {p0:g}")
        print(f"  perturbed start         = {p0*args.perturb:g}")
        print(f"  fitted (recovered)      = {fitted:g}")
        if fitted == fitted and p0:
            err = abs(fitted - p0) / p0
            print(f"  recovery error = {err:.1%}  "
                  + ("<- agent recovered the calibrated rate from the target"
                     if err < 0.1 else "(loose - widen bounds / check identifiability)"))
    finally:
        sb.stop()


if __name__ == "__main__":
    main()
