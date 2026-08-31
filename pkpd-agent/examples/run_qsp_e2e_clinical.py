r"""C) The agent's hub, dropped into the FULL model: does its choice change the CLINICAL endpoint?

A showed the IL-6 half-effect concentrations are NOT pinnable from steady state - the agent must
pick a sloppy solution (K = each regulator's own level, Hill-0.5). B showed modules compose into a
closed loop. C closes the last gap: take the agent's principled IL-6-hub choice, drop it into the
REAL 59-species sbproj, run the full virtual-population -> ACR clinical pipeline, and measure how far
the week-12 ACR rates move versus the paper's calibrated model. This is the endpoint that matters -
if the agent's unpinned-K choice barely moves ACR, the data bottleneck from A is clinically benign;
if it moves ACR a lot, the dose-response data is genuinely required.

    python -m examples.run_qsp_e2e_clinical --model ra ^
        --sbproj "..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj" ^
        --vpop "..\RA-QSP-Model\Vpop1.xlsx" --limit 200

Needs the MATLAB engine + the real sbproj/Vpop (your machine). The IL-6-secretion half-effect
parameter names are DISCOVERED at runtime by pattern and printed - confirm them on the first run
and refine --pattern if your model names them differently.
"""

from __future__ import annotations

import argparse
import json
import os
import re

from pkpd_agent.engines.simbiology import SimBiologyEngine
from pkpd_agent.engines import qsp_config

_DATA = os.path.join(os.path.dirname(__file__), "..", "projects", "vantage_ra", "data")


def _arm_rates(sb, cfg, vpop, first_line, readout_day, limit, overrides=""):
    out = {}
    for d, spec in cfg.drugs.items():
        doses = (spec or {}).get("doses") or []
        if not doses:
            continue
        r = sb.run_vpop(vpop, dose=doses[0], stop_time=readout_day + 2,
                        baseline_day=cfg.timeline.get("baseline_day", 200.0),
                        readout_day=readout_day, limit=limit,
                        param_overrides=overrides, states=cfg.readout_states or None)
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
    ap.add_argument("--pattern", default=r"(?i)IL6.*HalfEffect|HalfEffect.*IL6|IL6Sec.*HalfEffect",
                    help="regex selecting the IL-6-secretion half-effect params to override")
    ap.add_argument("--limit", type=int, default=200)
    args = ap.parse_args()

    cfg = qsp_config.get(args.model)
    first_line = cfg.run_columns.get("first_line") or {}
    readout_day = cfg.timeline.get("first_line_readout_day", 284.0)
    tg = {t["model_species"]: t for t in json.load(open(os.path.join(_DATA,
          "steady_state_targets.json"))) if t.get("model_species")}
    levels = {c: float(tg[c]["target_model_unit"]) for c in tg
              if tg[c].get("target_model_unit") is not None}

    sb = SimBiologyEngine()
    try:
        print("== starting MATLAB =="); sb.start()
        print(f"== loading {os.path.basename(args.sbproj)} =="); sb.load_project(args.sbproj)

        rgx = re.compile(args.pattern)
        params = sb.list_parameters().get("parameters", [])
        hub = [p for p in params if rgx.search(p["name"])]
        print(f"\n== IL-6-secretion half-effect params discovered ({len(hub)}) ==")
        # agent's principled choice: K = the regulator cytokine's own baseline level (Hill 0.5).
        overrides, mapped = [], []
        for p in hub:
            cyt = next((c for c in levels if re.search(rf"(?i){c}\b", p["name"])), None)
            if cyt is None:
                print(f"    {p['name']}  (no cytokine matched - left at calibrated value)")
                continue
            overrides.append(f"{p['name']}={levels[cyt]:g}")
            mapped.append((p["name"], cyt, float(p["value"]) if p["value"] else None, levels[cyt]))
            print(f"    {p['name']:40} -> K = level[{cyt}] = {levels[cyt]:g}  "
                  f"(calibrated was {p['value']})")
        if not overrides:
            print("  no overridable params matched --pattern; refine it and re-run."); return

        print(f"\n== baseline: paper's calibrated model, week-12 ACR per arm (limit {args.limit}) ==",
              flush=True)
        base = _arm_rates(sb, cfg, args.vpop, first_line, readout_day, args.limit)
        for d, r in base.items():
            print(f"  {d:<5} " + " ".join(f"{role}={r.get(role)}" for role in first_line))

        print(f"\n== agent's IL-6 hub (K = regulator level) dropped into the full model ==",
              flush=True)
        agent = _arm_rates(sb, cfg, args.vpop, first_line, readout_day, args.limit,
                           overrides=";".join(overrides))
        for d, r in agent.items():
            print(f"  {d:<5} " + " ".join(f"{role}={r.get(role)}" for role in first_line))

        print(f"\n== clinical shift: paper-calibrated vs agent's hub choice ==")
        print(f"  {'arm':<5} {'metric':<7} {'paper':>6} {'agent':>6} {'shift':>7}")
        shifts = []
        for d in base:
            for role in first_line:
                a, b = base[d].get(role), agent[d].get(role)
                if a is not None and b is not None:
                    shifts.append(abs(b - a)); print(f"  {d:<5} {role:<7} {a:>6} {b:>6} {b-a:>+7.1f}")
        if shifts:
            print(f"\n  mean absolute ACR shift: {sum(shifts)/len(shifts):.1f} points")
        print("\n  -> this is the CLINICAL cost of the agent's unpinned-K choice in the full closed-")
        print("     loop model. Small shift => the dose-response data A flagged is clinically benign")
        print("     here; large shift => it is genuinely required. Either way it is measured on the")
        print("     real endpoint, not a reduced module.")
    finally:
        sb.stop()


if __name__ == "__main__":
    main()
