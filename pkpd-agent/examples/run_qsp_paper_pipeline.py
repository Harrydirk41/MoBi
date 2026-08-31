r"""The whole paper workflow, one command, narrated end to end.

Runs the Vantage RA QSP process the way the paper does, and reports each step:

  PART I  - from scratch (agent-built subsystem, runs anywhere):
    1. BUILD an ODE          - the agent chooses the IL-6 hub structure from biology
    2. WRITE a simulation file - emit SBML + a runnable SimBiology .m
    3. CALIBRATE to a reference mean - fit the baseline rate so the model reproduces the target
  PART II - the full 59-species clinical model (needs MATLAB + the real sbproj/Vpop):
    4. BUILD / qualify a VIRTUAL POPULATION - check the baseline DAS28 distribution vs the target
    5. CALIBRATION check    - first-line MTX ACR20/50/70 vs the calibrated reference means
    6. QUALIFY against a held-out trial + SIMULATE THE SWITCH - escalate MTX-inadequate responders
       to second-line TCZ, read the second-line ACR, compare to the RADIATE trial
    7. REPORT the whole process

Honest scope seam: steps 1-3 prove the agent can BUILD and calibrate a model from nothing, but on
one subsystem - reconstructing all 59 species from scratch would need the answer model. Steps 4-6
are the population/clinical scale, which only exists in the full model, so they run on the real
sbproj (the paper's model). The command shows both halves as one continuous, honestly-labelled run.

    python -m examples.run_qsp_paper_pipeline --model ra ^
        --sbproj "..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj" ^
        --vpop "..\RA-QSP-Model\Vpop1.xlsx" --limit 300 [--live]

Without --sbproj it runs PART I only (the from-scratch build) and prints how to run PART II.
"""

from __future__ import annotations

import argparse
import json
import os
import statistics
import tempfile

from pkpd_agent.engines import model_assembly as MA
from pkpd_agent.engines.sbml_import import sbml_to_network
from pkpd_agent.engines import qsp_config
from examples.run_qsp_end_to_end import (assemble, integrate, generate_matlab_script,
                                         EXPERIMENTS, _RECORDED_AGENT)

_DATA = os.path.join(os.path.dirname(__file__), "..", "projects", "vantage_ra", "data")


def _load(name):
    return json.load(open(os.path.join(_DATA, name)))


def _frac(cols, flag, mask=None):
    vals = cols.get(flag, [])
    idx = list(range(len(vals)))
    if mask:
        m = cols.get(mask, [])
        idx = [i for i in idx if i < len(m) and isinstance(m[i], (int, float)) and m[i] >= 0.5]
    sub = [vals[i] for i in idx if i < len(vals)
           and isinstance(vals[i], (int, float)) and vals[i] == vals[i]]
    return (round(100.0 * sum(1 for v in sub if v >= 0.5) / len(sub), 1), len(sub)) if sub \
        else (None, 0)


def _cmp(model, ref):
    return "  ".join(f"{k}: model {model.get(k)}  ref {ref.get(k)}  "
                     f"(Δ{(model.get(k) - ref.get(k)):+.1f})"
                     if model.get(k) is not None and ref.get(k) is not None else f"{k}: -"
                     for k in ("ACR20", "ACR50", "ACR70"))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="ra")
    ap.add_argument("--sbproj")
    ap.add_argument("--vpop")
    ap.add_argument("--limit", type=int, default=300)
    ap.add_argument("--live", action="store_true", help="agent builds the structure via the LLM")
    args = ap.parse_args()

    cfg = qsp_config.get(args.model)
    calib = _load("calibration.json")
    valid = _load("validation.json")
    tg = {t["model_species"]: t for t in _load("steady_state_targets.json") if t.get("model_species")}
    levels = {c: float(tg[c]["target_model_unit"]) for c in tg
              if tg[c].get("target_model_unit") is not None}
    truth_maxes = {}
    for pref in ("IL6SecFLS_Maxby", "IL6SecMacro_Maxby"):
        prov = {p["name"]: p for p in _load("param_provenance.json")}
        for n, p in prov.items():
            if n.startswith(pref):
                c = n.split("Maxby")[-1]; v = p.get("value_from_reference")
                if v is not None and c in levels and c not in truth_maxes:
                    truth_maxes[c] = float(v)
    kcl = float({p["name"]: p for p in _load("param_provenance.json")}["kcl_IL6"]
                ["value_from_reference"])
    il6_target = levels["IL6"]

    print("#" * 78)
    print("# PART I - from scratch: build an ODE, write a sim file, calibrate  (subsystem)")
    print("#" * 78)

    # STAGE 1: agent builds the ODE structure
    print("\n== STAGE 1: BUILD an ODE - agent chooses the IL-6 hub structure from biology ==")
    call = None
    if args.live:
        from pkpd_agent.config import AgentConfig
        from pkpd_agent.engines import llm_tasks as LT
        c2 = AgentConfig(mock=False)
        call = LT.default_call(c2) if c2.anthropic_key_present() else None
    cyts = sorted(c for c, t in tg.items() if t.get("kind") == "cytokine" and c in levels)
    if call is not None:
        chosen = MA.propose_regulators("IL6", cyts, "secretion", call)
        motif = MA.propose_motif("IL6", [{"species": r["cytokine"]} for r in chosen], "", call)
        print("  (LIVE) agent's regulators + reasoning:")
        for r in chosen:
            print(f"    {r['cytokine']:6} {(r.get('direction') or ''):4} {r.get('basis') or ''}")
        print(f"  (LIVE) rate-law form: {motif.get('proliferation_order')}/"
              f"{motif.get('combination')} - {motif.get('reason') or ''}")
    else:
        chosen = _RECORDED_AGENT
        motif = {"proliferation_order": "zeroth", "combination": "product", "cap": None}
        print("  (recorded clean-agent choice; --live to call the LLM) "
              f"{[r['cytokine'] for r in chosen]}")

    # STAGE 2: write the simulation file(s)
    print("\n== STAGE 2: WRITE a simulation file - emit SBML + a runnable SimBiology .m ==")
    spec, fd, fp = assemble(chosen, truth_maxes, levels, kcl, il6_target, motif=motif)
    xml = os.path.join(tempfile.gettempdir(), "il6_hub.xml")
    open(xml, "w", encoding="utf-8").write(MA.to_sbml(spec))
    mm = os.path.join(tempfile.gettempdir(), "sim_il6_hub.m")
    open(mm, "w", encoding="utf-8").write(generate_matlab_script("il6_hub.xml", EXPERIMENTS))
    print(f"  SBML  -> {xml}   ({len(spec['species'])} species, {len(spec['reactions'])} reactions)")
    print(f"  MATLAB-> {mm}    (run `sim_il6_hub` in SimBiology)")

    # STAGE 3: calibrate to the reference mean
    print("\n== STAGE 3: CALIBRATE to a reference mean - fit baseline rate to the target ==")
    kg = next(p["value"] for p in spec["parameters"] if p["name"] == "kg_IL6")
    ss = integrate(sbml_to_network(xml), {c: levels[c] for c in levels if c != "IL6"})["IL6"]
    print(f"  fit kg_IL6 = {kg:.4g};  simulated steady state = {ss:.4g}  "
          f"(target mean {il6_target:g}, error {abs(ss-il6_target)/il6_target:.1%})")

    if not (args.sbproj and args.vpop):
        print("\n[PART II skipped: pass --sbproj and --vpop to run the full-model population/"
              "clinical stages on the real 59-species model.]")
        return

    print("\n" + "#" * 78)
    print("# PART II - full 59-species clinical model  (real sbproj; population + trial scale)")
    print("#" * 78)
    from pkpd_agent.engines.simbiology import SimBiologyEngine
    sb = SimBiologyEngine()
    try:
        print("\n== starting MATLAB & loading the model =="); sb.start(); sb.load_project(args.sbproj)
        second_day = cfg.timeline.get("second_line_readout_day", 600.0)
        first_day = cfg.timeline.get("first_line_readout_day", 284.0)

        # STAGE 4: virtual population - baseline DAS28 distribution vs the target
        print("\n== STAGE 4: BUILD/qualify the VIRTUAL POPULATION - baseline DAS28 vs target ==",
              flush=True)
        r0 = sb.run_vpop(args.vpop, dose="MTX_15mg_Q1W_SC_t200", stop_time=first_day + 2,
                         baseline_day=cfg.timeline.get("baseline_day", 200.0),
                         readout_day=first_day, limit=args.limit, states=cfg.readout_states or None)
        cols = r0.get("columns") or {}
        das = [v for v in cols.get("DAS28_BL", []) if isinstance(v, (int, float)) and v == v]
        vt = calib["vpop_target"]
        if das:
            print(f"  baseline DAS28: model mean {statistics.mean(das):.2f} "
                  f"sd {statistics.pstdev(das):.2f}  (target {vt['mean']}±{vt['sd']}, "
                  f"band {vt['band']}); n={len(das)}")

        # STAGE 5: calibration check - first-line MTX ACR vs calibrated means
        print("\n== STAGE 5: CALIBRATION check - first-line MTX ACR vs reference means ==")
        mtx = {k: _frac(cols, col)[0] for k, col in cfg.run_columns["first_line"].items()}
        ref = calib["calibrated_arms"][0]["known_rates_day284_full_pop_n300"]
        print(f"  MTX first-line: {_cmp(mtx, ref)}")

        # STAGE 6: qualify held-out + simulate the switch to second-line TCZ
        print("\n== STAGE 6: SIMULATE THE SWITCH - MTX-IR -> second-line TCZ, read second line ==",
              flush=True)
        rs = sb.run_vpop(args.vpop, dose="MTX_15mg_Q1W_SC_t200;TCZ8mgkg_Q4W_IV_t200@285",
                         stop_time=second_day + 2, baseline_day=cfg.timeline.get("baseline_day", 200.0),
                         readout_day=second_day, limit=args.limit, states=cfg.readout_states or None)
        sc = rs.get("columns") or {}
        sub = cfg.run_columns.get("subgroup_flag", "MTX_NonResp")
        second = {k: _frac(sc, col, mask=sub)[0]
                  for k, col in cfg.run_columns["second_line"].items() if k != "remission"}
        n_ir = _frac(sc, cfg.run_columns["second_line"]["ACR20"], mask=sub)[1]
        rt = valid["refractory_target"]
        print(f"  MTX-inadequate responders escalated to TCZ 8mg/kg (n={n_ir} in subgroup):")
        print(f"  second-line: {_cmp(second, rt)}")
        print(f"  (held-out reference: {rt['trial']})")

        # STAGE 7: report
        print("\n== STAGE 7: REPORT - the whole process, one run ==")
        print("  1 BUILD     agent chose the IL-6 hub structure from biology")
        print(f"  2 SIM FILE  emitted SBML + SimBiology .m")
        print(f"  3 CALIBRATE subsystem hits its steady-state target ({ss:.4g} vs {il6_target:g})")
        print(f"  4 VPOP      baseline DAS28 mean {statistics.mean(das):.2f} within target band "
              f"{vt['band']}" if das else "  4 VPOP      -")
        print(f"  5 CALIB CHK MTX first-line ACR20 {mtx.get('ACR20')} vs ref {ref.get('ACR20')}")
        print(f"  6 SWITCH    second-line TCZ ACR20 {second.get('ACR20')} vs RADIATE "
              f"{rt.get('ACR20')}")
        print("  -> one continuous run: build -> sim file -> calibrate -> Vpop -> qualify -> "
              "switch -> second-line.")
    finally:
        sb.stop()


if __name__ == "__main__":
    main()
