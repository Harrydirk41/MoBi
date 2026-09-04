r"""The WHOLE process, one command: from structured data to a clinical trial readout, then graded
against the ground-truth model. No project-specific scaffolding - every step reads the model's own
conventions and the given trial data.

Stages (each reported):
  1 BUILD MODEL    - the agent proposes the immune network's structure (topology + rate form) from
                     biology (--live), or the model's own edges offline; assemble + steady-state
                     calibrate every free rate. Reports topology recall/precision and calibration.
  2 SIM FILE       - transplant the agent network into the paper's clinical shell (given DAS28/PK/
                     dose model), re-wiring a small-molecule PD; emit one runnable sbproj. Reports
                     reactions swapped and baseline DAS28 finiteness.
  3 CALIBRATE      - the joint steady state (immune) is pinned to targets; the population baseline
                     DAS28 is qualified against the trial's observed distribution (stage 4).
  4 VIRTUAL POP    - sample the agent's severity knobs, keep patients whose baseline DAS28 is in the
                     trial band. Reports the qualified baseline DAS28 vs target.
  5 QUALIFY        - first-line MTX ACR20/50/70 at week 12 vs the calibrated reference.
  6 SIMULATE END   - the switch: MTX-inadequate responders -> second-line TCZ, ACR at the second
                     readout, vs the held-out RADIATE trial.
  7 GROUND TRUTH   - run the SAME first/second-line readout on the PAPER's own model + Vpop, so the
                     agent process is graded against the validated model, not only the trials.
  8 REPORT         - a full, honest comparison (agent vs ground-truth-model vs trial) written to a
                     markdown file, with the process metrics and the localized gaps.

    python -m examples.run_qsp_full_pipeline --sbproj "..\RA-QSP-Model\...sbproj" ^
        --modeldir "..\RA-QSP-Model" --paper-vpop "..\RA-QSP-Model\Vpop1.xlsx" ^
        --model ra --live --prune --rewire-mtx --n 300 --limit 150 --out report.md

Needs MATLAB + SimBiology + matlab.engine (the run_qsp_paper_pipeline --matlab setup); --live also
needs ANTHROPIC_API_KEY. Runs the pure build stage even without MATLAB.
"""

from __future__ import annotations

import argparse
import json
import os
import statistics

from pkpd_agent.engines import cell_lifecycle as CL, model_assembly as MA, network_assembly as NA
from examples.run_qsp_agent_vpop import (select_severity_params, _frac, _write_vpop_csv,
                                         _marginal_cells, _load_targets)
from examples.run_qsp_build_general import _project_dir


def _load(model):
    root = os.path.join(_project_dir(model), "data")
    prov = {p["name"]: p for p in json.load(open(os.path.join(root, "param_provenance.json")))}
    targets = json.load(open(os.path.join(root, "steady_state_targets.json")))
    levels = {t["model_species"]: float(t["target_model_unit"]) for t in targets
              if t.get("target_model_unit") is not None and t.get("model_species")}
    aliases = CL.load_cell_aliases(_project_dir(model))
    cells = CL.discover_cells(prov, targets, aliases)
    return prov, levels, cells, aliases


def build_model(model, live, prune, force_influx, cfg=None, stabilize=0, log=print):
    """STAGE 1: build + steady-state calibrate the agent immune network. Returns (spec, meta, topo)
    where topo summarizes agent topology recall/precision (or 'model edges' offline). With
    ``stabilize`` > 0, run the GUARDED inner loop: the agent revises its own structure from the
    dynamics (which species diverge) - never from the answer - before assembly."""
    prov, levels, cells, aliases = _load(model)
    for c in [x.strip() for x in (force_influx or "").split(",") if x.strip()]:
        if c in cells:
            CL.synthesize_influx(cells[c], levels)
    topo = {"mode": "model edges (offline)"}
    if live and cfg and cfg.anthropic_key_present():
        from pkpd_agent.engines import llm_tasks as LT
        call = LT.default_call(cfg)
        sec_struct, cell_struct, sc = NA.propose_structure(prov, levels, cells, aliases, call)
        model_cells = cells
        if stabilize > 0:
            log("  STABILIZE LOOP (agent revises its own structure from the dynamics, guarded):")
            sec_struct, cell_struct, hist = NA.stabilize_loop(
                prov, levels, cells, aliases, sec_struct, cell_struct, sc["confidence"],
                max_iters=stabilize, call=call, log=log)
            topo["stabilize_history"] = hist
        if prune:
            sec_struct, cell_struct, _ = NA.prune_structure(sec_struct, cell_struct,
                                                            sc["confidence"], prov, aliases,
                                                            model_cells)

        def agg(d):
            have = [v for v in d.values() if v["truth"]]
            return (round(sum(v["recall"] for v in have) / len(have), 2),
                    round(sum(v["precision"] for v in have) / len(have), 2)) if have else (0, 0)
        model_sec = NA.discover_secretion(prov, NA.cell_token_map(aliases))
        sec2, cells2 = NA.apply_structure(model_sec, cells, sec_struct, cell_struct)
        spec, meta = NA.assemble_network(prov, levels, cells2, aliases, sec_override=sec2)
        topo = {"mode": "agent (--live)" + (" +prune" if prune else ""),
                "secreting_cells": agg(sc["secreting_cells"]),
                "secretion_mods": agg(sc["secretion_mods"]), "cell_flux": agg(sc["cell_flux"])}
    else:
        spec, meta = NA.assemble_network(prov, levels, cells, aliases)

    targ = {s["name"]: s["initial"] for s in spec["species"]}
    ss = NA.integrate_network(spec, t_end=5.0, dt=5e-3)
    drift = max(abs(ss[k] - targ[k]) / targ[k] for k in targ)
    marg = sorted(c for c, (kp, m) in meta["free_kprolif"].items() if m)
    return spec, meta, {**topo, "calibration_drift": drift, "marginal_cells": marg,
                        "species": len(spec["species"]), "reactions": len(spec["reactions"])}


def _acr(cols, run_columns, line, mask=None):
    return {k: _frac(cols, col, mask=mask)[0] for k, col in run_columns[line].items()
            if k != "remission" or line == "first_line"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sbproj", required=True, help="the paper's real .sbproj (clinical shell + ground truth)")
    ap.add_argument("--modeldir", required=True, help="folder with the model's helper functions (MM.m)")
    ap.add_argument("--paper-vpop", dest="paper_vpop", default=None,
                    help="the paper's Vpop .xlsx, for the ground-truth-model comparison (stage 7)")
    ap.add_argument("--model", default="ra")
    ap.add_argument("--live", action="store_true", help="agent proposes the structure (needs a key)")
    ap.add_argument("--prune", action="store_true", help="drop low-confidence uncited agent edges")
    ap.add_argument("--stabilize", type=int, default=0, metavar="N",
                    help="run the guarded inner loop up to N iterations: the agent revises its own "
                         "structure from the dynamics (which species diverge), never from the answer")
    ap.add_argument("--rewire-mtx", action="store_true", dest="rewire_mtx",
                    help="re-attach MTX's PD onto the agent's secretion/influx reactions")
    ap.add_argument("--force-influx", default=None, dest="force_influx",
                    help="what-if: give these marginal cells a synthesized influx arm")
    ap.add_argument("--n", type=int, default=300)
    ap.add_argument("--span", type=float, default=2.0)
    ap.add_argument("--limit", type=int, default=150)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default="qsp_full_report.md")
    args = ap.parse_args()

    calib, valid, tasks = _load_targets(args.model)
    vt = calib["vpop_target"]
    ref = calib["calibrated_arms"][0]["known_rates_day284_full_pop_n300"]
    rt = valid["refractory_target"]
    rc = tasks["run_columns"]
    tl = tasks.get("timeline", {})
    b_day = float(tl.get("baseline_day", 200.0))
    r1 = float(tl.get("first_line_readout_day", 284.0))
    r2 = float(tl.get("second_line_readout_day", 600.0))
    mtx = tasks["drugs"]["MTX"]["doses"][0]
    tcz = tasks["drugs"]["TCZ"]["doses"][0]
    sub = rc.get("subgroup_flag", "MTX_NonResp")
    rep = {"stages": {}, "targets": {"vpop": vt, "first_line_ref": ref, "radiate": rt}}

    # ---- STAGE 1: BUILD MODEL (pure Python) ----
    from pkpd_agent.config import AgentConfig
    cfg = AgentConfig(mock=False)
    print("== STAGE 1: BUILD MODEL (agent structure + steady-state calibration) ==", flush=True)
    spec, meta, topo = build_model(args.model, args.live, args.prune, args.force_influx, cfg,
                                   stabilize=args.stabilize)
    net_xml = os.path.abspath("mynet.xml")
    open(net_xml, "w", encoding="utf-8").write(MA.to_sbml(spec))
    rep["stages"]["build"] = topo
    print(f"  {topo['mode']}: {topo['species']} species / {topo['reactions']} reactions; "
          f"calibration drift {topo['calibration_drift']:.3%}; marginal cells {topo['marginal_cells']}")
    for k in ("secreting_cells", "secretion_mods", "cell_flux"):
        if k in topo:
            print(f"  topology [{k:15}] recall {topo[k][0]} precision {topo[k][1]}")

    from pkpd_agent.engines.simbiology import SimBiologyEngine
    sb = SimBiologyEngine()
    try:
        sb.start()
    except Exception as e:                              # noqa: BLE001
        print(f"\n[no matlab.engine here: {str(e)[:80]}] - stage 1 done. Run on a MATLAB host "
              "for stages 2-8."); _write_report(args.out, rep, args); return
    try:
        if not sb.has_simbiology():
            print("\n[MATLAB started but SimBiology not licensed] - stage 1 done.")
            _write_report(args.out, rep, args); return
        sb.eng.addpath(os.path.abspath(args.modeldir), nargout=0)

        # ---- STAGE 2: SIM FILE (transplant + re-wire + save) ----
        print("\n== STAGE 2: SIM FILE (transplant into the clinical shell) ==", flush=True)
        sec_spec, influx = "", ""
        if args.rewire_mtx:
            sec_spec = "*=(1-Anti_CytSec_MTX);IL10=(1+Pro_CytSec_MTX);TGFb=(1+Pro_CytSec_MTX)"
            influx = "(1-Anti_CellInflux_MTX)"
        agent_proj = os.path.abspath("agent_clinical.sbproj")
        sb.eng.sb_agent_clinical(os.path.abspath(args.sbproj), net_xml, agent_proj,
                                 "DAS28_CRP", "", sec_spec, influx, nargout=0)
        rep["stages"]["sim_file"] = {"sbproj": agent_proj, "rewire_mtx": args.rewire_mtx}

        # ---- STAGES 3-6: agent Vpop + first/second line ----
        agent = run_clinical(sb, agent_proj, args, tasks, vt, b_day, r1, r2, mtx, tcz, sub, rc,
                             label="AGENT", modeldir=args.modeldir)
        rep["stages"]["agent_clinical"] = agent

        # ---- STAGE 7: GROUND TRUTH (paper model + paper Vpop) ----
        truth = None
        if args.paper_vpop and os.path.isfile(args.paper_vpop):
            print("\n== STAGE 7: GROUND TRUTH (paper model + paper Vpop) ==", flush=True)
            sb.load_project(os.path.abspath(args.sbproj))
            c1 = (sb.run_vpop(os.path.abspath(args.paper_vpop), dose=mtx, stop_time=r1 + 2,
                              baseline_day=b_day, readout_day=r1, limit=args.limit,
                              states=tasks.get("readout_states")).get("columns") or {})
            c2 = (sb.run_vpop(os.path.abspath(args.paper_vpop),
                              dose=f"{mtx};{tcz}@{int(r1)+1}", stop_time=r2 + 2, baseline_day=b_day,
                              readout_day=r2, limit=args.limit,
                              states=tasks.get("readout_states")).get("columns") or {})
            truth = {"first_line": _acr(c1, rc, "first_line"),
                     "second_line": _acr(c2, rc, "second_line", mask=sub)}
            print(f"  paper model first-line MTX:  {truth['first_line']}")
            print(f"  paper model second-line TCZ: {truth['second_line']}")
        rep["stages"]["ground_truth"] = truth

        # ---- STAGE 8: REPORT ----
        _write_report(args.out, rep, args)
        print(f"\n== STAGE 8: full report -> {os.path.abspath(args.out)} ==")
    finally:
        sb.stop()


def run_clinical(sb, proj, args, tasks, vt, b_day, r1, r2, mtx, tcz, sub, rc, label, modeldir):
    """STAGES 3-6 for one model: sample+qualify a Vpop, then first-line MTX and the TCZ switch."""
    sb.load_project(proj)
    sb.eng.addpath(os.path.abspath(modeldir), nargout=0)
    params = sb.list_parameters()["parameters"]
    marginal = _marginal_cells(args.model)
    spec, chosen = select_severity_params(params, marginal, args.span)
    print(f"\n== STAGE 4: VIRTUAL POP ({label}) - sample {args.n}, qualify to band {vt['band']} ==",
          flush=True)
    res = sb.sample_vpop(spec, n_samples=args.n, baseline_day=b_day, seed=args.seed)
    cols = res.get("columns") or {}
    das = cols.get("DAS28_base") or []
    lo, hi = vt["band"]
    qual = [i for i, d in enumerate(das) if isinstance(d, (int, float)) and d == d and lo <= d <= hi]
    qdas = [das[i] for i in qual]
    print(f"  qualified {len(qual)}/{len(das)}; baseline DAS28 mean "
          f"{statistics.mean(qdas):.2f} (target {vt['mean']}±{vt['sd']})" if qdas else "  none qualified")
    if not qual:
        return {"qualified": 0}
    present = [n for n in chosen if n in cols]
    csvp, xlsx = os.path.abspath("agent_vpop.csv"), os.path.abspath("agent_vpop.xlsx")
    _write_vpop_csv(csvp, present, [{n: cols[n][i] for n in present} for i in qual])
    sb.eng.sb_csv_to_xlsx(csvp, xlsx, nargout=0)

    print(f"== STAGE 5: QUALIFY ({label}) - first-line MTX at day {r1:g} ==", flush=True)
    c1 = (sb.run_vpop(xlsx, dose=mtx, stop_time=r1 + 2, baseline_day=b_day, readout_day=r1,
                      limit=args.limit, states=tasks.get("readout_states")).get("columns") or {})
    first = _acr(c1, rc, "first_line")
    print(f"  {label} first-line MTX: {first}")
    print(f"== STAGE 6: SIMULATE END ({label}) - MTX-IR -> TCZ at day {r2:g} ==", flush=True)
    c2 = (sb.run_vpop(xlsx, dose=f"{mtx};{tcz}@{int(r1)+1}", stop_time=r2 + 2, baseline_day=b_day,
                      readout_day=r2, limit=args.limit,
                      states=tasks.get("readout_states")).get("columns") or {})
    second = _acr(c2, rc, "second_line", mask=sub)
    n_ir = _frac(c2, rc["second_line"]["ACR20"], mask=sub)[1]
    print(f"  {label} second-line TCZ (n_IR={n_ir}): {second}")
    return {"qualified": len(qual), "baseline_das28": round(statistics.mean(qdas), 2),
            "first_line": first, "second_line": second, "n_ir": n_ir}


def _row(name, agent, truth, trial):
    def f(x):
        return "-" if x is None else (f"{x:.1f}" if isinstance(x, (int, float)) else str(x))
    return f"| {name} | {f(agent)} | {f(truth)} | {f(trial)} |"


def _write_report(path, rep, args):
    b = rep["stages"].get("build", {})
    ag = rep["stages"].get("agent_clinical") or {}
    gt = rep["stages"].get("ground_truth") or {}
    ref, radiate = rep["targets"]["first_line_ref"], rep["targets"]["radiate"]
    L = ["# QSP from-scratch pipeline vs ground truth", "",
         f"Model: **{args.model}**  |  build: **{b.get('mode','?')}**  |  "
         f"rewire-mtx: {args.rewire_mtx}  |  force-influx: {args.force_influx or 'none'}", "",
         "## Process metrics (stage 1-2)", "",
         f"- network: {b.get('species','?')} species / {b.get('reactions','?')} reactions; "
         f"steady-state calibration drift **{b.get('calibration_drift',0):.3%}**",
         f"- marginal (birth-death) cells: {b.get('marginal_cells')}"]
    for k, lab in [("secreting_cells", "secreting cells"), ("secretion_mods", "secretion mods"),
                   ("cell_flux", "cell flux")]:
        if k in b:
            L.append(f"- agent topology [{lab}]: recall **{b[k][0]}**, precision **{b[k][1]}**")
    af, asec = ag.get("first_line", {}), ag.get("second_line", {})
    gf, gsec = gt.get("first_line", {}) if gt else {}, gt.get("second_line", {}) if gt else {}
    L += ["", "## Clinical readout: agent vs ground-truth model vs trial", "",
          f"agent baseline DAS28 **{ag.get('baseline_das28','-')}** (target {rep['targets']['vpop']['mean']})",
          "", "| readout | agent | paper model | trial |", "|---|---|---|---|",
          _row("first-line MTX ACR20", af.get("ACR20"), gf.get("ACR20"), ref.get("ACR20")),
          _row("first-line MTX ACR50", af.get("ACR50"), gf.get("ACR50"), ref.get("ACR50")),
          _row("first-line MTX ACR70", af.get("ACR70"), gf.get("ACR70"), ref.get("ACR70")),
          _row("2nd-line TCZ ACR20", asec.get("ACR20"), gsec.get("ACR20"), radiate.get("ACR20")),
          _row("2nd-line TCZ ACR50", asec.get("ACR50"), gsec.get("ACR50"), radiate.get("ACR50")),
          _row("2nd-line TCZ ACR70", asec.get("ACR70"), gsec.get("ACR70"), radiate.get("ACR70")),
          "", "## Evaluation", "",
          "The agent built the immune network from structured data with no project scaffolding, "
          "calibrated it, wore the paper's given clinical shell, and produced a population trial "
          "readout. Distance from the paper model and the trials is the measured cost of the "
          "agent's structural precision and of any data gaps (e.g. a cell with no literature influx "
          "rate is built marginal, blunting an influx-suppressing drug).", "",
          f"_raw stages_: ```{json.dumps(rep['stages'], default=str)[:1500]}```"]
    open(path, "w", encoding="utf-8").write("\n".join(L))


if __name__ == "__main__":
    main()
