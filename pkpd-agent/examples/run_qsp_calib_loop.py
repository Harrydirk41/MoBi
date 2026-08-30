r"""Closed-loop calibration: the LLM diagnoses non-identifiability and ACQUIRES the isolating
experiment for each coupled parameter, round by round - the try-and-respond process the
one-shot bounds experiment lacked.

The one-shot run showed 5 FLS-proliferation regulators fit to ONE aggregate FLS target land
~124% off, even with the LLM's (correct) bounds - because bounds add no equation. Here the
loop instead ACQUIRES an equation per parameter: a single-cytokine perturbation (elevate only
IL6 -> the FLS response isolates FLSProlif_MaxbyIL6). Each round the LLM picks the next
experiment; a 1-D solve pins that regulator; recovery error is expected to fall toward ~0 as
experiments accumulate. If it does, the earlier failure was the one-shot DESIGN, not the LLM.

    python -m examples.run_qsp_calib_loop --model ra ^
        --sbproj "..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj" ^
        --targets projects\vantage_ra\data\steady_state_targets.json ^
        --cell FLS --regulators --high-fold 10 --llm

--llm lets the LLM choose each experiment (needs ANTHROPIC_API_KEY); without it the loop
auto-acquires them in order (isolates "does the loop work" from "does the LLM pick well").
The perturbation "data" is generated from the shipped model (the in-silico stand-in for the
in-vitro experiments the modellers ran). Needs the MATLAB engine.
"""

from __future__ import annotations

import argparse
import json
import os
import re

from pkpd_agent.engines.simbiology import SimBiologyEngine
from pkpd_agent.engines import calib_loop as CL


def _mean_err(est: dict, truth: dict) -> float:
    e = [abs(est[n] - truth[n]) / truth[n] for n in truth if truth[n]]
    return sum(e) / len(e) if e else float("nan")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="ra")
    ap.add_argument("--sbproj", required=True)
    ap.add_argument("--targets", required=True)
    ap.add_argument("--cell", default="FLS")
    ap.add_argument("--params", default="")
    ap.add_argument("--regulators", action="store_true")
    ap.add_argument("--high-fold", type=float, default=10.0,
                    help="elevate each cytokine to this multiple of its baseline for the "
                         "isolating experiment")
    ap.add_argument("--readout-day", type=float, default=199.0)
    ap.add_argument("--decouple", action="store_true",
                    help="clamp the OTHER regulators' cytokines to ~0 during each experiment "
                         "(in-vitro condition: the cell + only one cytokine, no feedback)")
    ap.add_argument("--llm", action="store_true", help="let the LLM pick each experiment")
    ap.add_argument("--llm-model", default=None)
    args = ap.parse_args()

    sb = SimBiologyEngine()
    try:
        print("== starting MATLAB =="); sb.start()
        print(f"== loading {os.path.basename(args.sbproj)} =="); sb.load_project(args.sbproj)
        params = sb.list_parameters().get("parameters", [])
        pmap = {p["name"]: float(p["value"]) for p in params}

        if args.params:
            names = [n.strip() for n in args.params.split(",") if n.strip()]
        else:
            rgx = re.compile(rf"(?i){args.cell}.*prolif.*max")
            names = [p["name"] for p in params if rgx.search(p["name"])]
        names = [n for n in names if n in pmap and CL.regulator_cytokine(n)]
        if len(names) < 2:
            print(f"need >=2 regulators with a parseable cytokine; got {names}")
            return
        truth = {n: pmap[n] for n in names}
        cyt = {n: CL.regulator_cytokine(n) for n in names}
        print("coupled regulators and their isolating cytokine:")
        for n in names:
            print(f"    {n}  (truth {truth[n]:g})  <- elevate {cyt[n]}")

        # baseline cytokine levels -> perturbation 'high' value for each isolating experiment
        base = sb.simulate(stop_time=args.readout_day + 1.0).get("columns", {})
        high = {}
        for n in names:
            c = cyt[n]
            if c not in base:
                print(f"cytokine {c} not a model state; skipping {n}")
            else:
                high[n] = float(base[c][-1]) * args.high_fold or args.high_fold

        all_cyt = {cyt[n] for n in names}

        def others_of(n):                              # cytokines to clamp to ~0 (decouple)
            return sorted(all_cyt - {cyt[n]}) if args.decouple else None

        # truth targets: the isolating experiment run on the shipped (calibrated) model
        print(f"\n== generating isolating-experiment data from the shipped model "
              f"(elevate each cytokine {args.high_fold:g}x"
              + (", others clamped to 0 = in-vitro) ==" if args.decouple else ") =="),
              flush=True)
        target = {}
        for n in names:
            if n not in high:
                continue
            target[n] = sb.perturb_response(cyt[n], high[n], args.cell, args.readout_day,
                                            decouple=others_of(n))
            print(f"    {cyt[n]} elevated -> {args.cell} = {target[n]:g}")
        names = [n for n in names if n in target]

        # start all regulators at 1.0 (no measured effect yet); pin them one by one
        for n in names:
            sb.set_parameter(n, 1.0)
        est = {n: 1.0 for n in names}
        pinned, remaining = [], list(names)
        print(f"\n== loop start: all regulators = 1.0 (unpinned); "
              f"mean recovery error {_mean_err(est, truth):.0%} ==")

        picker = None
        if args.llm:
            from pkpd_agent.config import AgentConfig
            from pkpd_agent.engines import llm_tasks as LT
            cfg = AgentConfig(mock=False)
            if args.llm_model:
                cfg.model = args.llm_model
            if not cfg.anthropic_key_present():
                print("ANTHROPIC_API_KEY not set; falling back to auto order.")
            else:
                picker = LT.default_call(cfg)

        rounds = 0
        while remaining and rounds < len(names) + 2:
            rounds += 1
            state = {"remaining": [cyt[n] for n in remaining],
                     "pinned": [cyt[n] for n in pinned],
                     "mean_error": f"{_mean_err(est, truth):.0%}"}
            if picker:
                act = CL.choose_experiment(state, picker)
            else:
                act = {"action": "experiment", "cytokine": cyt[remaining[0]]}
            if act["action"] == "stop":
                print(f"  round {rounds}: LLM chose STOP"); break
            chosen = next((n for n in remaining if cyt[n] == act["cytokine"]), None)
            if chosen is None:
                break

            # isolating solve: pin this regulator from its single-cytokine experiment
            def evaluator(m, _n=chosen):
                sb.set_parameter(_n, m)
                return sb.perturb_response(cyt[_n], high[_n], args.cell, args.readout_day,
                                           decouple=others_of(_n))
            solved = CL.solve_1d(evaluator, target[chosen], lo=0.1, hi=20.0)
            sb.set_parameter(chosen, solved)
            est[chosen] = solved
            pinned.append(chosen); remaining.remove(chosen)
            print(f"  round {rounds}: acquired {cyt[chosen]} experiment -> pinned {chosen} "
                  f"= {solved:.3g} (truth {truth[chosen]:g}); "
                  f"mean recovery error now {_mean_err(est, truth):.0%}", flush=True)

        for n in names:                                # restore
            sb.set_parameter(n, truth[n])
        print(f"\n== result ==")
        print(f"  {'parameter':30} {'truth':>8} {'pinned':>8} {'error':>7}")
        for n in names:
            err = abs(est[n] - truth[n]) / truth[n] if truth[n] else float("nan")
            print(f"  {n:30} {truth[n]:8.3g} {est[n]:8.3g} {err:6.0%}")
        print(f"\n  mean recovery error: one-shot bounds ~124%  ->  closed loop "
              f"{_mean_err(est, truth):.0%}")
        print("  -> if this is small, the earlier failure was the one-shot DESIGN "
              "(no experiment\n     acquisition), not the LLM. Each parameter needs its "
              "own isolating experiment;\n     the loop's job is to realize that and go get them.")
    finally:
        sb.stop()


if __name__ == "__main__":
    main()
