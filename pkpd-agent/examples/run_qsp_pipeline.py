r"""End-to-end virtual-population pipeline, from the model's parameters to a Vpop - no
pre-curated driver list. Chains the general pieces:

  1. enumerate EVERY model parameter                         (sb_params)
  2. LLM selects the varied set by inferred category          (propose_vpop_set)
  3. assign default bounds (a +/- fold band around each       (heuristic: ranges are an
     parameter's shipped value; log scale)                     input, so propose a default)
  4. sample a cohort over that (high-dim) set under the arm   (sb_cohort)
  5. select a Vpop to match the clinical anchors + ESS        (weighting, optional GA)

The question it answers: with the FULL (~200-parameter) driver set the LLM picks - rather
than a 7-parameter toy - does the Vpop selection converge (healthy ESS, achieved rate near
target)? Calibration is NOT re-done here: the .sbproj already ships a calibrated reference
patient, so this uses it as-is (automated re-calibration to the clinical means is a
separate, unbuilt step).

    python -m examples.run_qsp_pipeline --model ra ^
        --sbproj "..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj" --n 400 --k 10 --ga

Needs ANTHROPIC_API_KEY and the MATLAB engine.
"""

from __future__ import annotations

import argparse
import collections
import os
import tempfile

from pkpd_agent.config import AgentConfig
from pkpd_agent.engines.simbiology import SimBiologyEngine
from pkpd_agent.engines import qsp_config, qsp_tasks, llm_tasks as LT


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="ra")
    ap.add_argument("--sbproj", required=True)
    ap.add_argument("--n", type=int, default=400, help="cohort size")
    ap.add_argument("--k", type=float, default=10.0,
                    help="default bound fold: vary each param over [v/k, v*k] (log)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--ga", action="store_true", help="also select a Vpop via native GA")
    ap.add_argument("--qualify", action="store_true",
                    help="predict the held-out trial: run the realized Vpop under the "
                         "flagship switch protocol and compare the second-line response to "
                         "the refractory target (an extra ~Vpop-size long simulations).")
    ap.add_argument("--llm-model", default=None)
    args = ap.parse_args()

    cfg = qsp_config.get(args.model)
    anchors_cfg = cfg.vpop_anchors or {}
    arms = anchors_cfg.get("arms") or {}
    rate_targets = anchors_cfg.get("rate_targets") or {}
    if not arms:
        print(f"project '{args.model}' declares no vpop_anchors.arms - nothing to match.")
        return

    cfg_llm = AgentConfig(mock=False)
    if args.llm_model:
        cfg_llm.model = args.llm_model
    if not cfg_llm.anthropic_key_present():
        print("ANTHROPIC_API_KEY not set.")
        return

    sb = SimBiologyEngine()
    try:
        print("== starting MATLAB engine =="); sb.start()
        print(f"== loading {os.path.basename(args.sbproj)} =="); sb.load_project(args.sbproj)

        # 1. enumerate every parameter
        print("== [1] sb_params: enumerating every model parameter ==", flush=True)
        params = (sb.list_parameters().get("parameters") or [])
        print(f"   model exposes {len(params)} parameters")
        if not params:
            print("   no parameters - aborting."); return

        # 2. LLM selects the varied set
        print("== [2] LLM selecting the varied set by inferred category ==", flush=True)
        sel = LT.propose_vpop_set(params, LT.default_call(cfg_llm))
        chosen = set(sel["selected"])
        print(f"   selected {sel['n_selected']} / {sel['n_candidates']} parameters")

        # 3. default bounds: a +/- fold band around each shipped value (log), skipping
        #    non-positive values (structurally inactive / not log-varyable)
        value = {p["name"]: p.get("value") for p in params}
        bounds = {}
        for nm in sel["selected"]:
            v = value.get(nm)
            if isinstance(v, (int, float)) and v > 0:
                bounds[nm] = (v / args.k, v * args.k, "log")
        print(f"   assigned +/-{args.k:g}x log bounds to {len(bounds)} of them "
              f"(dropped {len(chosen) - len(bounds)} non-positive/zero)")
        print("   NOTE: calibration not re-done - using the model's shipped reference patient")
        if not bounds:
            print("   no varyable parameters - aborting."); return

        # 4. sample a high-dimensional cohort under the arm(s)
        spec = qsp_tasks.build_sample_spec(bounds)
        arms_spec = ";;".join(f"{lab}:{dose}" for lab, dose in arms.items())
        baseline_day = cfg.timeline.get("baseline_day", 200.0)
        readout_day = cfg.timeline.get("first_line_readout_day", 284.0)
        print(f"== [3] sb_cohort: {args.n} candidates x {len(bounds)} params, "
              f"arms {list(arms)} ==", flush=True)
        r = sb.cohort_multi_arm(spec, arms_spec, baseline_day, readout_day, args.n,
                                args.seed, states=cfg.readout_states or None)
        ml = (r.get("matlab_log") or "").strip()
        if ml:
            print("   [MATLAB] " + ml.replace("\n", "\n   [MATLAB] "))
        cols = r.get("columns") or {}
        sevs = cols.get("sev_base", [])
        if not sevs:
            print("   cohort returned no rows - check the MATLAB log."); return
        print(f"   cohort: {len(sevs)} candidates, sev_base range "
              f"{min(sevs):.2f}..{max(sevs):.2f}")

        # plausibility gate (the paper's): keep only patients whose baseline severity is in
        # the active-disease band, dropping implausible ones BEFORE matching.
        band = cfg.vpop_target.get("band")
        if band:
            f = qsp_tasks.filter_columns_to_band(cols, "sev_base", band)
            cols = f["columns"]
            sevs = cols.get("sev_base", [])
            print(f"   plausibility gate [{band[0]}, {band[1]}]: kept {f['n_kept']} / "
                  f"{f['n_total']} candidates")
            if not sevs:
                print("   nobody plausible - widen bounds or the band."); return

        for lab in arms:
            col = [v for v in cols.get(lab, []) if isinstance(v, (int, float)) and v == v]
            rate = 100.0 * sum(1 for v in col if v >= 0.5) / len(col) if col else None
            print(f"   arm {lab}: raw response rate {rate}% (of plausible)")

        # 5. select a Vpop to match the anchors (weighting), and optionally native GA.
        #    Match the full response DISTRIBUTION per arm when rate_targets_full is given
        #    (primary threshold = readout_states[0] lives in the '<arm>' column; other
        #    thresholds in '<arm>__<state>'), else just the primary rate per arm.
        full = anchors_cfg.get("rate_targets_full")
        primary_thr = (cfg.readout_states or [None])[0]
        resp = {}   # column name -> target rate
        if full:
            for lab, thrs in full.items():
                for thr, tgt in thrs.items():
                    colname = lab if thr == primary_thr else f"{lab}__{thr}"
                    if colname in cols:
                        resp[colname] = tgt
        else:
            resp = {lab: rate_targets[lab] for lab in rate_targets if lab in cols}
        print(f"   matching {len(resp)} response anchors: {sorted(resp)}")

        anchors = [{"key": "severity", "mean": cfg.vpop_target["mean"],
                    "sd": cfg.vpop_target.get("sd")}]
        anchors += [{"key": c, "target": t} for c, t in resp.items()]
        candidates = [{"severity": sevs[i],
                       **{c: (cols[c][i] if i < len(cols[c]) else None) for c in resp}}
                      for i in range(len(sevs))]
        wsel = qsp_tasks.select_multi_anchor(candidates, anchors)
        print("== [4] Vpop selection (weighting) ==")
        if wsel.get("ok"):
            for a in wsel["anchors"]:
                print(f"   {a['key']:>10}: target {a['target']}  ->  achieved {a['achieved']}")
            print(f"   ESS {wsel['effective_sample_size']} "
                  f"({int(wsel['ess_fraction']*100)}% of {wsel['n']})")
        else:
            print("   failed:", wsel.get("reason"))

        if args.ga:
            spec_a = [f"moment:sev_base:{cfg.vpop_target['mean']}:{cfg.vpop_target.get('sd','')}"]
            spec_a += [f"rate:{lab}:{rate_targets[lab]}" for lab in rate_targets]
            print("== [4b] Vpop selection (native GA) ==")
            g = None
            try:
                g = sb.select_ga(cols, ";".join(spec_a), pop_target=0)
            except Exception as e:
                print(f"   GA unavailable/failed: {e}")
                print("   (needs the Global Optimization Toolbox; the weighting result "
                      "above still stands)")
            if g is not None:
                gml = (g.get("matlab_log") or "").strip()
                if gml:
                    print("   [MATLAB] " + gml.replace("\n", "\n   [MATLAB] "))
                gcols = g.get("columns") or {}
                gsev = gcols.get("sev_base", [])
                ns = g.get("n_selected", len(gsev))
                if ns:
                    gm = sum(gsev) / len(gsev) if gsev else 0
                    print(f"   selected {ns} / {len(sevs)} candidates")
                    print(f"   {'severity':>10}: target {cfg.vpop_target['mean']}  ->  "
                          f"achieved {round(gm,2)}")
                    for lab in rate_targets:
                        gc = [v for v in gcols.get(lab, [])
                              if isinstance(v, (int, float)) and v == v]
                        gr = 100.0 * sum(1 for v in gc if v >= 0.5) / len(gc) if gc else None
                        print(f"   {lab:>10}: target {rate_targets[lab]}  ->  achieved "
                              f"{round(gr,1) if gr is not None else None}")
                else:
                    print("   GA selected nobody - check the MATLAB log.")

        # realize a discrete Vpop from the prevalence weights (the paper's final step)
        if wsel.get("ok") and wsel.get("weights"):
            vp = qsp_tasks.realize_vpop(wsel["weights"], size=300, seed=args.seed)
            idx = vp["indices"]
            if idx:
                rsev = [sevs[i] for i in idx if i < len(sevs)]
                rmean = sum(rsev) / len(rsev) if rsev else 0
                print("== [5] realized Vpop (enrich high-weight patients) ==")
                print(f"   drew {vp['size']} patients ({vp['unique']} distinct) from the "
                      f"weighted pool")
                print(f"   {'severity':>10}: target {cfg.vpop_target['mean']}  ->  "
                      f"realized {round(rmean,2)}")
                for lab in rate_targets:
                    rc = [cols.get(lab, [None]*len(sevs))[i] for i in idx if i < len(sevs)]
                    rc = [v for v in rc if isinstance(v, (int, float)) and v == v]
                    rr = 100.0 * sum(1 for v in rc if v >= 0.5) / len(rc) if rc else None
                    print(f"   {lab:>10}: target {rate_targets[lab]}  ->  realized "
                          f"{round(rr,1) if rr is not None else None}")

            # [6] QUALIFY: run the realized Vpop under the flagship switch protocol and
            # compare its second-line (held-out) response to the refractory target.
            if args.qualify and idx and cfg.refractory_target:
                print("== [6] qualify: predict the held-out trial with this Vpop ==",
                      flush=True)
                try:
                    import csv as _csv
                    counts = collections.Counter(idx)          # weight = times drawn
                    uidx = sorted(counts)
                    pnames = list(bounds.keys())
                    xlsx = os.path.join(tempfile.gettempdir(), "vpop_realized.csv")
                    with open(xlsx, "w", newline="", encoding="utf-8") as fh:
                        wtr = _csv.writer(fh)
                        wtr.writerow(pnames)
                        for i in uidx:
                            wtr.writerow([cols[p][i] if i < len(cols.get(p, [])) else ""
                                          for p in pnames])
                    dose = ";".join(cfg.flagship_protocol.get("first_line", []) +
                                    cfg.flagship_protocol.get("second_line", []))
                    stop = cfg.timeline.get("second_line_readout_day", 600.0) + 100
                    print(f"   running {len(uidx)} distinct patients under '{dose}' to "
                          f"day {stop:g} ...", flush=True)
                    rv = sb.run_vpop(xlsx, dose=dose, stop_time=stop,
                                     baseline_day=baseline_day, readout_day=readout_day,
                                     states=cfg.readout_states or None)
                    rml = (rv.get("matlab_log") or "").strip()
                    if rml:
                        print("   [MATLAB] " + rml.replace("\n", "\n   [MATLAB] "))
                    rc = rv.get("columns") or {}
                    flag = cfg.run_columns.get("subgroup_flag")
                    second = cfg.run_columns.get("second_line") or {}
                    sub = rc.get(flag, [])
                    w = [counts[i] for i in uidx]      # multiplicity, aligned to output rows
                    m = min(len(w), len(sub))
                    print("   second-line response in the subgroup vs the held-out trial:")
                    for role, tgt in cfg.refractory_target.items():
                        colname = second.get(role)
                        if not colname or colname not in rc or not isinstance(tgt, (int, float)):
                            continue
                        resp = rc[colname]
                        mm = min(m, len(resp))
                        den = sum(w[j] * (sub[j] or 0) for j in range(mm))
                        num = sum(w[j] * (sub[j] or 0) * (resp[j] or 0) for j in range(mm))
                        pred = 100.0 * num / den if den else None
                        print(f"   {role:>7}: predicted "
                              f"{round(pred,1) if pred is not None else None}  vs observed {tgt}")
                except Exception as e:
                    print(f"   qualify failed: {e}")

        print("\n== verdict ==")
        ess = wsel.get("ess_fraction") if wsel.get("ok") else 0
        raw = None
        for lab in rate_targets:
            col = [v for v in cols.get(lab, []) if isinstance(v, (int, float)) and v == v]
            raw = 100.0 * sum(1 for v in col if v >= 0.5) / len(col) if col else None
        print(f"   {len(bounds)}-parameter Vpop: raw arm response {raw}%, "
              f"weighting ESS {int((ess or 0)*100)}%. "
              f"{'converging' if (ess or 0) > 0.2 else 'still degenerate - responder pool too thin'}")
    finally:
        sb.stop()


if __name__ == "__main__":
    main()
