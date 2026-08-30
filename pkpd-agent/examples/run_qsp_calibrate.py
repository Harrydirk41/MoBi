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
import random
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


def _sim_cell(sb, cell, steady_day):
    return float(sb.simulate(stop_time=steady_day + 1.0)["columns"][cell][-1])


def _recover_multi(sb, names, truths, cell, s0, steady_day, restarts, perturb, seed):
    """Jointly fit several coupled parameters to the SAME single steady-state target, from
    several random starts. One scalar target cannot pin many parameters, so each restart lands
    a DIFFERENT parameter vector that reproduces the target equally well - that spread IS the
    non-identifiability. Returns (per-restart fitted vectors, per-restart output value)."""
    rng = random.Random(seed)
    data = os.path.join(tempfile.gettempdir(), "calib_multi_target.csv")
    with open(data, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["Time", cell]); w.writerow([steady_day, s0])
    spec = ";".join(f"{n},{t/50:g},{t*50:g},log" for n, t in zip(names, truths))
    fits, outs = [], []
    for r in range(restarts):
        for n, t in zip(names, truths):                # random start in [1/perturb, perturb]*truth
            sb.set_parameter(n, t * (perturb ** rng.uniform(-1, 1)))
        res = sb.fit_native(spec, data, f"{cell} = {cell}", method="lsqnonlin")
        est = res.get("columns", {}).get("estimate", [])
        fitted = [float(x) for x in est] if len(est) == len(names) else [float("nan")] * len(names)
        for n, v in zip(names, fitted):                # apply the fit, read the output it gives
            if v == v:
                sb.set_parameter(n, v)
        outs.append(_sim_cell(sb, cell, steady_day))
        fits.append(fitted)
        print(f"  restart {r+1}/{restarts}: output {cell} = {outs[-1]:g} "
              f"(target {s0:g})", flush=True)
    for n, t in zip(names, truths):
        sb.set_parameter(n, t)                         # restore
    return fits, outs


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="ra")
    ap.add_argument("--sbproj", required=True)
    ap.add_argument("--targets", required=True, help="steady_state_targets.json")
    ap.add_argument("--cell", default="FLS", help="model species to calibrate (a cell)")
    ap.add_argument("--param", default="", help="override the proliferation parameter name")
    ap.add_argument("--params", default="", help="comma list: JOINTLY fit these (identifiability)")
    ap.add_argument("--regulators", action="store_true",
                    help="auto-pick the cell's proliferation regulators (<cell>Prolif_Maxby*) "
                         "and jointly fit them - the identifiability experiment")
    ap.add_argument("--restarts", type=int, default=4,
                    help="joint fit from N random starts (spread of fits = non-identifiability)")
    ap.add_argument("--seed", type=int, default=1)
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

        # ---- Phase 2b: identifiability experiment (joint fit of coupled params) ----
        if args.params or args.regulators:
            if args.params:
                jnames = [n.strip() for n in args.params.split(",") if n.strip()]
            else:                                      # auto: the cell's proliferation regulators
                rgx = re.compile(rf"(?i){args.cell}.*prolif.*max")
                jnames = [p["name"] for p in params if rgx.search(p["name"])]
            truths = [next((float(p["value"]) for p in params if p["name"] == n), None)
                      for n in jnames]
            jnames = [n for n, t in zip(jnames, truths) if t is not None]
            truths = [t for t in truths if t is not None]
            if len(jnames) < 2:
                print(f"need >=2 joint parameters; found {jnames}. Pass --params a,b,c.")
                return
            print(f"\n== identifiability experiment: jointly fit {len(jnames)} coupled "
                  f"regulators to ONE target ({args.cell}={s0:g}), {args.restarts} restarts ==")
            for n, t in zip(jnames, truths):
                print(f"    {n} (truth {t:g})")
            fits, outs = _recover_multi(sb, jnames, truths, args.cell, s0, args.steady_day,
                                        args.restarts, args.perturb, args.seed)
            print(f"\n== result: output fits, parameters DON'T (non-identifiability) ==")
            out_ok = sum(1 for o in outs if s0 and abs(o - s0) / s0 < 0.05)
            print(f"  output match: {out_ok}/{len(outs)} restarts reproduce {args.cell} "
                  f"within 5%")
            print(f"  {'parameter':32} {'truth':>10}  {'fitted range across restarts':>28}"
                  f"  spread")
            for j, (n, t) in enumerate(zip(jnames, truths)):
                col = [f[j] for f in fits if f[j] == f[j]]
                if not col:
                    continue
                lo, hi = min(col), max(col)
                spread = (hi - lo) / t if t else float("nan")
                print(f"  {n:32} {t:10.3g}  [{lo:10.3g}, {hi:10.3g}]  {spread:5.1f}x")
            print("\n  -> many parameter sets reproduce the one target equally well: a single "
                  "steady-state\n     value cannot pin coupled regulators. THIS is why the full "
                  "527-param fit needs\n     many constraints + biological judgment, not a "
                  "push-button optimizer.")
            return

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
