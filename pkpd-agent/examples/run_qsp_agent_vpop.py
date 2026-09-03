r"""B) Virtual population over the AGENT model's OWN parameters -> ACR20/50/70 vs the trials.

The paper's Vpop1.xlsx varies the PAPER's parameters (kg_FLS_Baseline, F_TNFa, ...), which the
transplant removed - so it cannot drive the agent's immune mechanism. This builds a Vpop over the
agent's OWN severity knobs instead: the free rates the agent calibrated - the cytokine secretion
scales ksec_<Cyt> and the (non-marginal) cell proliferation rates kprolif_<Cell>. Marginal cells
(no influx: FLS, Macrophages, PlasmaCells) are excluded because perturbing their proliferation has
no restoring force and just destabilises a candidate.

Workflow (each step reuses the SAME MATLAB helpers the paper model uses, so the trial readout is the
model's own event-set ACR flags, never recomputed):

  1. list the agent params in agent_clinical.sbproj; build a log-uniform span around each nominal.
  2. SAMPLE candidates and simulate each to its untreated baseline DAS28 (sb_sample_vpop).
  3. QUALIFY: keep candidates whose baseline DAS28 falls in the target band (the paper's observed
     baseline DAS28 distribution) - the population-qualification step.
  4. write the qualified patients to a Vpop .xlsx (row 1 = agent param names).
  5. TEST first line: run MTX -> ACR20/50/70 at day 284, compare to the calibrated reference.
  6. SIMULATE the switch: MTX-inadequate responders -> second-line TCZ, read ACR at day 600,
     compare to the held-out RADIATE trial.

    python -m examples.run_qsp_agent_vpop --sbproj agent_clinical.sbproj --model ra --n 300 --limit 150

Needs MATLAB + SimBiology + matlab.engine (the run_qsp_paper_pipeline --matlab setup) and openpyxl.
Honest expectation: the agent model over-includes edges, so the ACR fit will be worse than the
paper's - that gap, measured against the trials, is the point.
"""

from __future__ import annotations

import argparse
import json
import os
import re


def _load_targets(model):
    from examples.run_qsp_build_general import _project_dir
    d = os.path.join(_project_dir(model), "data")
    calib = json.load(open(os.path.join(d, "calibration.json")))
    valid = json.load(open(os.path.join(d, "validation.json")))
    tasks = json.load(open(os.path.join(_project_dir(model), "tasks.json")))
    return calib, valid, tasks


def _marginal_cells(model):
    from pkpd_agent.engines import cell_lifecycle as CL
    from examples.run_qsp_build_general import _project_dir
    d = os.path.join(_project_dir(model), "data")
    prov = {p["name"]: p for p in json.load(open(os.path.join(d, "param_provenance.json")))}
    targets = json.load(open(os.path.join(d, "steady_state_targets.json")))
    levels = {t["model_species"]: float(t["target_model_unit"]) for t in targets
              if t.get("target_model_unit") is not None and t.get("model_species")}
    cells = CL.discover_cells(prov, targets, CL.load_cell_aliases(_project_dir(model)))
    return {c for c in cells if CL.fit_base_prolif(cells[c], levels)[1]}


def select_severity_params(params, marginal, span):
    """From the model's parameters, pick the agent severity knobs (ksec_*, kprolif_<non-marginal>)
    and build a ';'-joined 'name,lo,hi,log' spec with a log-uniform span around each nominal."""
    spec, chosen = [], []
    for p in params:
        n, v = p.get("name", ""), p.get("value")
        if v is None or v <= 0:
            continue
        m = re.match(r"(?i)^kprolif_(.+)$", n)
        is_ksec = n.startswith("ksec_")
        if is_ksec or (m and m.group(1) not in marginal):
            spec.append(f"{n},{v / span:.6g},{v * span:.6g},log")
            chosen.append(n)
    return ";".join(spec), chosen


def _frac(cols, flag, mask=None):
    vals = cols.get(flag, [])
    idx = list(range(len(vals)))
    if mask:
        mk = cols.get(mask, [])
        idx = [i for i in idx if i < len(mk) and isinstance(mk[i], (int, float)) and mk[i] >= 0.5]
    sub = [vals[i] for i in idx if i < len(vals)
           and isinstance(vals[i], (int, float)) and vals[i] == vals[i]]
    return (round(100.0 * sum(1 for v in sub if v >= 0.5) / len(sub), 1), len(sub)) if sub \
        else (None, 0)


def _write_vpop_xlsx(path, names, rows):
    import openpyxl
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.append(list(names))
    for r in rows:
        ws.append([r[n] for n in names])
    wb.save(path)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sbproj", required=True, help="the agent-based clinical sbproj")
    ap.add_argument("--model", default="ra")
    ap.add_argument("--n", type=int, default=300, help="candidate patients to sample")
    ap.add_argument("--span", type=float, default=2.0, help="log-uniform fold-range per param")
    ap.add_argument("--limit", type=int, default=0, help="subsample this many qualified patients to run")
    ap.add_argument("--seed", type=int, default=1)
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
    doses = tasks["drugs"]
    mtx = doses["MTX"]["doses"][0]
    tcz = next(d for d in doses["TCZ"]["doses"])
    marginal = _marginal_cells(args.model)

    from pkpd_agent.engines.simbiology import SimBiologyEngine
    sb = SimBiologyEngine()
    print("== starting MATLAB & loading the agent-based sbproj ==", flush=True)
    sb.start()
    try:
        sb.load_project(os.path.abspath(args.sbproj))
        params = sb.list_parameters()["parameters"]
        spec, chosen = select_severity_params(params, marginal, args.span)
        print(f"  agent severity knobs ({len(chosen)}): ksec_* + kprolif_<non-marginal>")
        print(f"  excluded marginal-cell kprolif: {sorted(marginal)}")

        # 2-3. SAMPLE + QUALIFY against the baseline DAS28 band
        print(f"\n== SAMPLE {args.n} candidates -> baseline DAS28, qualify to band {vt['band']} ==",
              flush=True)
        res = sb.sample_vpop(spec, n_samples=args.n, baseline_day=b_day, seed=args.seed)
        cols = res.get("columns") or {}
        das = cols.get("DAS28_base") or cols.get("DAS28_CRP") or []
        lo, hi = vt["band"]
        qual = [i for i, d in enumerate(das)
                if isinstance(d, (int, float)) and d == d and lo <= d <= hi]
        import statistics
        qdas = [das[i] for i in qual]
        print(f"  sampled {len(das)}, qualified {len(qual)} in band; "
              f"baseline DAS28 mean {statistics.mean(qdas):.2f} sd {statistics.pstdev(qdas):.2f} "
              f"(target {vt['mean']}±{vt['sd']})" if qdas else "  no candidate qualified")
        if not qual:
            print("  -> widen --span or --n; the agent baseline may sit outside the band."); return

        # 4. write the qualified patients to a Vpop xlsx
        xlsx = os.path.join(os.path.dirname(os.path.abspath(args.sbproj)), "agent_vpop.xlsx")
        rows = [{n: cols[n][i] for n in chosen if n in cols} for i in qual]
        _write_vpop_xlsx(xlsx, chosen, rows)
        print(f"  wrote {len(rows)} patients -> {xlsx}")

        # 5. TEST first-line MTX
        print(f"\n== TEST first-line MTX -> ACR at day {r1:g} ==", flush=True)
        r_first = sb.run_vpop(xlsx, dose=mtx, stop_time=r1 + 2, baseline_day=b_day,
                              readout_day=r1, limit=args.limit, states=tasks.get("readout_states"))
        c1 = r_first.get("columns") or {}
        first = {k: _frac(c1, col)[0] for k, col in rc["first_line"].items()}
        print(f"  agent MTX first-line: {first}")
        print(f"  paper reference:      {ref}")

        # 6. SIMULATE the switch to second-line TCZ
        print(f"\n== SIMULATE switch: MTX-IR -> TCZ, ACR at day {r2:g} (vs RADIATE) ==", flush=True)
        r_sec = sb.run_vpop(xlsx, dose=f"{mtx};{tcz}@{int(r1)+1}", stop_time=r2 + 2,
                            baseline_day=b_day, readout_day=r2, limit=args.limit,
                            states=tasks.get("readout_states"))
        c2 = r_sec.get("columns") or {}
        sub = rc.get("subgroup_flag", "MTX_NonResp")
        second = {k: _frac(c2, col, mask=sub)[0] for k, col in rc["second_line"].items()
                  if k != "remission"}
        n_ir = _frac(c2, rc["second_line"]["ACR20"], mask=sub)[1]
        print(f"  agent second-line TCZ (n_IR={n_ir}): {second}")
        print(f"  RADIATE held-out:      {{'ACR20': {rt['ACR20']}, 'ACR50': {rt['ACR50']}, "
              f"'ACR70': {rt['ACR70']}}}")
        print("\n  -> the from-scratch agent model, wearing the paper's clinical shell, produces a "
              "population\n     ACR response. Its distance from the trials is the measured cost of "
              "the agent's precision.")
    finally:
        sb.stop()


if __name__ == "__main__":
    main()
