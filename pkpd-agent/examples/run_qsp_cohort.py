r"""Smoke test for sb_cohort.m + multi-anchor Vpop selection - no API, no agent.

Samples a small virtual cohort (each candidate simulated under the project's therapy
arms via sb_cohort.m), then optimizes selection weights to match SEVERAL clinical
anchors at once (baseline severity + each arm's response rate). If the cohort table
prints and the selection reports achieved values near the targets, the .m + the
multi-anchor path work end to end.

    python -m examples.run_qsp_cohort --model ra ^
        --sbproj "..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj" --n 40
"""

from __future__ import annotations

import argparse
import os

from pkpd_agent.engines.simbiology import SimBiologyEngine
from pkpd_agent.engines import qsp_config, qsp_tasks


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="ra")
    ap.add_argument("--sbproj", required=True)
    ap.add_argument("--n", type=int, default=40, help="cohort size")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--narrow", type=float, default=None, metavar="K",
                    help="instead of each driver's full observed span, sample a band "
                         "[nominal/K, nominal*K] (clipped to the span). K>1; smaller K = "
                         "tighter band around the literature nominal. Use it to test the "
                         "hypothesis that narrower/less-severe bounds raise the ESS.")
    ap.add_argument("--low", type=float, default=None, metavar="FRAC",
                    help="sample only the LOWER FRAC fraction of each driver's log-span "
                         "([lo, lo*(hi/lo)^FRAC]). FRAC in (0,1]; the low end = low "
                         "amplification = milder disease = the MTX responders. Overrides "
                         "--narrow.")
    ap.add_argument("--ga", action="store_true",
                    help="also select a Vpop with the NATIVE genetic algorithm "
                         "(sb_select_ga.m) - picks a real subset of candidates, the "
                         "paper's method, to compare against the weighting.")
    ap.add_argument("--pop-target", type=int, default=0,
                    help="desired GA-selected population size (0 = unconstrained).")
    args = ap.parse_args()

    cfg = qsp_config.get(args.model)
    anchors_cfg = cfg.vpop_anchors or {}
    arms = anchors_cfg.get("arms") or {}
    rate_targets = anchors_cfg.get("rate_targets") or {}
    if not arms:
        print(f"project '{args.model}' declares no vpop_anchors.arms - nothing to do.")
        return

    # sample every disease driver over its observed span (log scale) -- or, with
    # --low FRAC, only the lower fraction of the span (milder = the responders); or,
    # with --narrow K, a tight band [nominal/K, nominal*K] around the literature nominal.
    def _bounds(p: dict) -> tuple:
        lo, hi = p["span"][0], p["span"][1]
        if args.low and 0 < args.low <= 1 and hi > lo > 0:
            hi = lo * (hi / lo) ** args.low
        elif args.narrow and args.narrow > 1 and p.get("nominal"):
            nom = float(p["nominal"])
            lo = max(lo, nom / args.narrow)
            hi = min(hi, nom * args.narrow)
        return (lo, hi, "log")

    if args.low:
        print(f"== sampling the lower {args.low:g} fraction of each driver's log-span "
              f"(milder disease) ==")
    elif args.narrow:
        print(f"== narrowing driver bounds to nominal/{args.narrow:g} .. "
              f"nominal*{args.narrow:g} (clipped to span) ==")
    spec = qsp_tasks.build_sample_spec(
        {n: _bounds(p) for n, p in cfg.vpop_drivers.items()})
    arms_spec = ";;".join(f"{lab}:{dose}" for lab, dose in arms.items())
    baseline_day = cfg.timeline.get("baseline_day", 200.0)
    readout_day = cfg.timeline.get("first_line_readout_day", 284.0)

    sb = SimBiologyEngine()
    try:
        print("== starting MATLAB engine =="); sb.start()
        print(f"== loading {os.path.basename(args.sbproj)} =="); sb.load_project(args.sbproj)
        print(f"== sb_cohort: {args.n} candidates x {len(cfg.vpop_drivers)} drivers, "
              f"arms {list(arms)} ==", flush=True)
        r = sb.cohort_multi_arm(spec, arms_spec, baseline_day, readout_day, args.n,
                                args.seed, states=cfg.readout_states or None)
        ml = (r.get("matlab_log") or "").strip()
        if ml:
            print("   [MATLAB] " + ml.replace("\n", "\n   [MATLAB] "))
        cols = r.get("columns") or {}
        sevs = cols.get("sev_base", [])
        print(f"\ncohort table: {len(sevs)} candidates")
        print("  sev_base: n=%d, range %.2f..%.2f" %
              (len(sevs), min(sevs) if sevs else 0, max(sevs) if sevs else 0))
        for lab in arms:
            col = [v for v in cols.get(lab, []) if isinstance(v, (int, float)) and v == v]
            rate = 100.0 * sum(1 for v in col if v >= 0.5) / len(col) if col else None
            print(f"  arm {lab}: raw response rate {rate}% (n={len(col)})")

        candidates = [{"severity": sevs[i],
                       **{lab: cols.get(lab, [None] * len(sevs))[i] for lab in arms}}
                      for i in range(len(sevs))]
        anchors = [{"key": "severity", "mean": cfg.vpop_target["mean"],
                    "sd": cfg.vpop_target.get("sd")}]
        anchors += [{"key": lab, "target": rate_targets[lab]} for lab in rate_targets]
        sel = qsp_tasks.select_multi_anchor(candidates, anchors)
        print("\n== multi-anchor selection ==")
        if not sel.get("ok"):
            print("  failed:", sel.get("reason"))
        else:
            for a in sel["anchors"]:
                print(f"  {a['key']:>10}: target {a['target']}  ->  achieved {a['achieved']}")
            print(f"  total error {sel['total_error']}, ESS {sel['effective_sample_size']} "
                  f"({int(sel['ess_fraction']*100)}% of {sel['n']})")

        # -- native GA subset selection (the paper's method): pick a real Vpop -- #
        if args.ga:
            spec = [f"moment:sev_base:{cfg.vpop_target['mean']}:{cfg.vpop_target.get('sd','')}"]
            spec += [f"rate:{lab}:{rate_targets[lab]}" for lab in rate_targets]
            g = sb.select_ga(cols, ";".join(spec), pop_target=args.pop_target)
            gml = (g.get("matlab_log") or "").strip()
            if gml:
                print("   [MATLAB] " + gml.replace("\n", "\n   [MATLAB] "))
            gcols = g.get("columns") or {}
            gsev = gcols.get("sev_base", [])
            ns = g.get("n_selected", len(gsev))
            print("\n== native GA selection ==")
            if not ns:
                print("  GA selected nobody - check the MATLAB log above.")
            else:
                gm = sum(gsev) / len(gsev) if gsev else 0
                print(f"  selected {ns} / {len(sevs)} candidates")
                print(f"  {'severity':>10}: target {cfg.vpop_target['mean']}  ->  "
                      f"achieved {round(gm, 2)}")
                for lab in rate_targets:
                    gc = [v for v in gcols.get(lab, []) if isinstance(v, (int, float)) and v == v]
                    gr = 100.0 * sum(1 for v in gc if v >= 0.5) / len(gc) if gc else None
                    print(f"  {lab:>10}: target {rate_targets[lab]}  ->  achieved "
                          f"{round(gr, 1) if gr is not None else None}")
    finally:
        sb.stop()


if __name__ == "__main__":
    main()
