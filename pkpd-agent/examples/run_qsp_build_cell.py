r"""Template A from-scratch build - a CELL life-cycle node, general, no per-cell scaffolding.

The cytokine-secretion builder (run_qsp_build_general) covers template B. This is its template-A
twin: a cell whose count is set by three fluxes - recruitment (influx), proliferation, and death -
each modulated by its own set of cytokines. Everything is discovered from the model's conventions:

  * WHICH cells are buildable    - a steady-state target + a kd_<Cell>_Baseline death rate
  * The THREE regulator sets     - <Cell>Prolif_/Influx_/Apop_Maxby<Cyt> (proliferation / influx /
                                   apoptosis) - so a cell poses THREE topology questions, not one
  * The BASELINE rates           - kd (death) and kIn (influx) looked up; the base proliferation
                                   rate is the FREE parameter, fitted to the steady-state target
  * The HELD-OUT experiment      - knock the strongest proliferation driver, watch the cell respond

    python -m examples.run_qsp_build_cell --model ra --list
    python -m examples.run_qsp_build_cell --model ra --cell Th1
    python -m examples.run_qsp_build_cell --model ra --cell FLS --live   # agent picks each flux
    python -m examples.run_qsp_build_cell --model ra --cell Treg --fit   # provenance split

Pure offline (uses the model's own regulator sets); --live needs a key. Emits a real SBML per cell.
"""

from __future__ import annotations

import argparse
import os
import tempfile

from pkpd_agent.engines import cell_lifecycle as CL, model_assembly as MA
from pkpd_agent.engines.sbml_import import sbml_to_network
from examples.run_qsp_build_general import load_model, _project_dir
from examples.run_qsp_end_to_end import integrate

_FLUX_PROCESS = {"prolif": "proliferation", "influx": "recruitment (influx from blood)",
                 "apop": "apoptosis (death)"}
_MOTIF = {"proliferation_order": "first", "combination": "product", "cap": None}


def _load_cells(model):
    import json
    root = os.path.join(_project_dir(model), "data")
    prov = {p["name"]: p for p in json.load(open(os.path.join(root, "param_provenance.json")))}
    targets = json.load(open(os.path.join(root, "steady_state_targets.json")))
    levels = {t["model_species"]: float(t["target_model_unit"]) for t in targets
              if t.get("target_model_unit") is not None and t.get("model_species")}
    aliases = CL.load_cell_aliases(_project_dir(model))
    cells = CL.discover_cells(prov, targets, aliases)
    return prov, levels, cells


def _agent_flux_regulators(cell, info, levels, call):
    """--live: the agent proposes each flux's regulators independently (three topology questions),
    scored against the model's own set for that flux."""
    cands = sorted(c for c in levels if c != cell)
    chosen = {}
    for flux, process in _FLUX_PROCESS.items():
        if not info[flux]:                                 # model gives this cell no such flux set
            chosen[flux] = set()
            continue
        regs = MA.propose_regulators(cell, cands, process, call)
        picks = {r["cytokine"] for r in regs}
        truth = {c for c in info[flux] if c in levels}
        hit = truth & picks
        rec = len(hit) / len(truth) if truth else 0.0
        prec = len(hit) / len(picks) if picks else 0.0
        print(f"    {flux:7} agent {sorted(picks)}")
        print(f"    {'':7} model {sorted(truth)}  recall {rec:.2f} precision {prec:.2f}")
        chosen[flux] = picks
    return chosen


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="ra")
    ap.add_argument("--cell", default=None)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--live", action="store_true", help="agent proposes each flux's regulators")
    ap.add_argument("--fit", action="store_true", help="print the parameter-origin split")
    args = ap.parse_args()

    prov, levels, cells = _load_cells(args.model)
    if args.list or not cells:
        print(f"== buildable CELL nodes in '{args.model}' (template A, no hardcoding) ==")
        for c, info in sorted(cells.items()):
            kp, marg = CL.fit_base_prolif(info, levels)
            print(f"  {c:12} prolif={sorted(info['prolif'])} influx={sorted(info['influx'])} "
                  f"apop={sorted(info['apop'])}{'  [marginal]' if marg else ''}")
        return

    cell = args.cell or sorted(cells)[0]
    if cell not in cells:
        print(f"'{cell}' not buildable; discovered: {sorted(cells)}"); return
    info = cells[cell]
    print(f"== building CELL '{cell}' (target {info['target']:g}) - template A ==")
    print(f"  model's regulator sets: prolif={sorted(info['prolif'])} "
          f"influx={sorted(info['influx'])} apop={sorted(info['apop'])}")

    # ---- STRUCTURE ----
    chosen = None
    if args.live:
        from pkpd_agent.config import AgentConfig
        from pkpd_agent.engines import llm_tasks as LT
        cfg = AgentConfig(mock=False)
        if not cfg.anthropic_key_present():
            print("  --live but no key; using the model's own regulator sets.")
        else:
            print("  (LIVE) agent proposes each flux's regulators:")
            chosen = _agent_flux_regulators(cell, info, levels, LT.default_call(cfg))

    # ---- ASSEMBLE + FIT + EMIT ----
    kpp = f"kprolif_{cell}"
    per_flux_chosen = None
    if chosen is not None:
        per_flux_chosen = set().union(*chosen.values()) if chosen else set()
    kp, marg = CL.fit_base_prolif(info, levels, chosen=per_flux_chosen)
    rxns, vals = CL.cell_reactions(cell, info, levels, kpp, chosen=per_flux_chosen)
    vals[kpp] = kp
    regcyt = sorted(c for c in CL.all_regulators(info, chosen=per_flux_chosen) if c in levels)
    species = [{"name": cell, "initial": info["target"]}] + \
              [{"name": c, "initial": levels[c], "boundary": True} for c in regcyt]
    spec = {"name": f"{cell}_lifecycle", "species": species,
            "parameters": [{"name": k, "value": v} for k, v in vals.items()],
            "reactions": rxns, "rules": []}
    xml = os.path.join(tempfile.gettempdir(), f"cell_{cell}.xml")
    open(xml, "w", encoding="utf-8").write(MA.to_sbml(spec))
    net = sbml_to_network(xml)
    print(f"\n  assembled: {len(species)} species, {len(rxns)} fluxes "
          f"({'influx+' if any(r['id'].endswith('_influx') for r in rxns) else ''}prolif+death)")
    print(f"  fitted base proliferation rate kprolif_{cell} = {kp:.4g}"
          f"{'  [marginal: births=deaths at the target IC]' if marg else ''}")
    print(f"  emitted {xml}")

    # ---- SIMULATE: hold-at-target + held-out knockdown of the strongest prolif driver ----
    clamp = {c: levels[c] for c in regcyt}
    stay = integrate(net, clamp, t_end=200.0, dt=5e-3)[cell]
    print(f"\n  STEADY STATE (calibration self-consistency): {stay:.4g} "
          f"(target {info['target']:g}, error {abs(stay-info['target'])/info['target']:.1%})")
    up = {c: (m or 1.5) for c, (m, _l) in info["prolif"].items() if c in levels}
    if up:
        drv = max(up, key=lambda c: abs(up[c] - 1))
        held = dict(clamp); held[drv] = levels[drv] * 0.1
        new = integrate(net, held, t_end=200.0, dt=5e-3)[cell]
        arrow = "no restoring force (marginal): drifts" if marg else "settles to a new steady state"
        print(f"  HELD-OUT anti-{drv} (strongest prolif driver -> 10%): {cell} = {new:.4g} "
              f"({(new-info['target'])/info['target']:+.0%}) - {arrow}")

    # ---- PROVENANCE ----
    if args.fit:
        kd_lit = prov[info["kd_param"]].get("from_literature")
        kin_lit = info["kin_lit"]
        print(f"\n== parameter origin for '{cell}' (template A) ==")
        print(f"  LOOK UP  kd_{cell} (death/turnover): from_literature={kd_lit}")
        if info["kin_param"]:
            print(f"  {'LOOK UP ' if kin_lit and info['kin_val'] else 'FIT     '} "
                  f"kIn_{cell} (influx): from_literature={kin_lit}, value="
                  f"{info['kin_val']}")
        print(f"  FIT      kprolif_{cell} (no citation) -> pinned by the one steady-state target")
        print(f"  FIT      every K_<cyt> half-effect -> under-determined by steady state alone "
              "(same data bottleneck as template B)")


if __name__ == "__main__":
    main()
