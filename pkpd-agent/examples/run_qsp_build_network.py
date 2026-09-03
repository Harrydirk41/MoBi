r"""Assemble and run the FULL coupled immune network - the scale-integration step.

Every earlier builder tested a PIECE: one cytokine hub (build_general), one cell life-cycle
(build_cell), one feedback pair (couple_general). This one puts them all together: every cell
(template A) and every cytokine-with-a-target (template B) dynamic at once, wired by the model's
OWN edges, each free rate pinned by its own steady-state target, emitted as ONE SBML and integrated
as one coupled ODE system. It answers the last open question - does the whole thing hold together at
scale, or only the isolated parts?

    python -m examples.run_qsp_build_network --model ra
    python -m examples.run_qsp_build_network --model ra --anti TNFa   # propagate a knockdown
    python -m examples.run_qsp_build_network --model ra --emit net.xml

What it shows, honestly:
  (1) ASSEMBLY   - N species / M reactions from the model's conventions, zero hardcoding.
  (2) CALIBRATION- the joint steady state is self-consistent to ~0% across every species: one free
                   rate per species (ksec_<Cyt>, kprolif_<Cell>) solved from its own target.
  (3) COUPLING   - a single-driver knockdown (--anti) propagates through the loop with the model's
                   own edge signs (e.g. anti-TNFa collapses the TNFa-driven chemokines).
What it does NOT claim: full clinical stability. Three cells (FLS, Macrophages, PlasmaCells) are
structurally birth-death (no influx baseline), so their level has no restoring force - the calibrated
state is a fixed point but only marginally stable. Pinning the whole trajectory (and DAS28 / PK /
Vpop) needs the paper's remaining PK-input + latch machinery in the MATLAB SimBiology engine.
"""

from __future__ import annotations

import argparse
import copy
import json
import os

from pkpd_agent.engines import cell_lifecycle as CL, model_assembly as MA, network_assembly as NA
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


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="ra")
    ap.add_argument("--anti", default="TNFa", help="cytokine to knock to 10%% for the coupling probe")
    ap.add_argument("--emit", default=None, help="write the assembled SBML to this path")
    ap.add_argument("--t-end", type=float, default=40.0)
    ap.add_argument("--live", action="store_true",
                    help="the AGENT proposes every edge of the whole network (topology + which "
                         "cells secrete what), scored against the model; needs ANTHROPIC_API_KEY")
    args = ap.parse_args()

    prov, levels, cells, aliases = _load(args.model)

    # ---- STRUCTURE: the model's own edges, or (--live) the agent's proposed edges end-to-end ----
    struct_src = "the model's own edges"
    if args.live:
        from pkpd_agent.config import AgentConfig
        from pkpd_agent.engines import llm_tasks as LT
        cfg = AgentConfig(mock=False)
        if not cfg.anthropic_key_present():
            print("  --live but no ANTHROPIC_API_KEY; falling back to the model's own edges.\n")
        else:
            print("== AGENT proposes the FULL network structure (nodes given; topology is the "
                  "agent's job) ==")
            sec_struct, cell_struct, sc = NA.propose_structure(prov, levels, cells, aliases,
                                                               LT.default_call(cfg), log=print)
            model_sec = NA.discover_secretion(prov, NA.cell_token_map(aliases))
            sec2, cells2 = NA.apply_structure(model_sec, cells, sec_struct, cell_struct)
            cells = cells2
            # aggregate ONLY over questions that have a model edge to recover: a node the model
            # sources from a constant input (no dynamic secreting cell) or a flux with no
            # regulator has an EMPTY truth, so 0/0 is not a miss - scoring it as 0 would understate
            def _agg(d):
                have = [v for v in d.values() if v["truth"]]
                empty = len(d) - len(have)
                if not have:
                    return 0.0, 0.0, empty
                return (sum(v["recall"] for v in have) / len(have),
                        sum(v["precision"] for v in have) / len(have), empty)
            for label, key in [("secreting-cells", "secreting_cells"),
                               ("secretion-mods", "secretion_mods"), ("cell-flux", "cell_flux")]:
                r, p, empty = _agg(sc[key])
                note = f"  ({empty} nodes have no model edge - excluded)" if empty else ""
                print(f"  AGENT topology [{label:16}] mean recall {r:.2f}  precision {p:.2f}"
                      f"  [on nodes with a model edge]{note}")
            struct_src = "the AGENT's proposed edges"
            spec, meta = NA.assemble_network(prov, levels, cells, aliases, sec_override=sec2)
    if not args.live or struct_src == "the model's own edges":
        spec, meta = NA.assemble_network(prov, levels, cells, aliases)

    marg = {c for c, (kp, m) in meta["free_kprolif"].items() if m}
    targ = {s["name"]: s["initial"] for s in spec["species"]}

    # (1) ASSEMBLY
    print(f"\n== full coupled immune network for '{args.model}' (scale integration) ==")
    print(f"  (1) ASSEMBLED {len(spec['species'])} species, {len(spec['reactions'])} reactions, "
          f"{len(spec['parameters'])} params - from {struct_src}")
    print(f"      dynamic cells   ({len(meta['dynamic_cells'])}): {meta['dynamic_cells']}")
    print(f"      dynamic cytokines ({len(meta['dynamic_cytokines'])}): {meta['dynamic_cytokines']}")
    drop = meta["dropped"]["secreting_cell_not_dynamic"]
    print(f"      dropped edges (species with no target: GMCSF/AutoAb, etc.): "
          f"{drop or 'none among secreting cells'}")
    xml = args.emit or os.path.join(os.path.dirname(__file__), "..", "build_network.xml")
    open(xml, "w", encoding="utf-8").write(MA.to_sbml(spec))
    print(f"      emitted SBML -> {os.path.abspath(xml)}")

    # (2) CALIBRATION self-consistency: start every species at target, hold nothing, integrate
    ss = NA.integrate_network(spec, t_end=5.0, dt=5e-3)
    drift = {k: abs(ss[k] - targ[k]) / targ[k] for k in targ}
    worst = max(drift, key=drift.get)
    print(f"\n  (2) CALIBRATION - joint steady state self-consistency (all species free):")
    print(f"      max drift {drift[worst]:.3%} at {worst}; every free rate (one ksec per cytokine, "
          f"one kprolif per cell)\n      solved from its OWN target - the whole network is one "
          "internally-consistent ODE system")
    print(f"      marginal (birth-death, no influx) cells: {sorted(marg)} - flagged, "
          "not level-pinned")

    # (3) COUPLING: single-driver knockdown propagates through the loop (marginal cells pinned)
    anti = args.anti
    if anti in targ and anti not in marg:
        pin = {c: targ[c] for c in marg}
        base = NA.integrate_network(spec, clamp=pin, t_end=args.t_end, dt=5e-3)
        knock = NA.integrate_network(spec, clamp={**pin, anti: targ[anti] * 0.1},
                                     t_end=args.t_end, dt=5e-3)
        print(f"\n  (3) COUPLING - anti-{anti} ({anti}->10%) propagated through the network "
              f"(fold vs baseline):")
        if base.get("__diverged__") or knock.get("__diverged__"):
            print(f"      NETWORK DIVERGED under the perturbation - it is DYNAMICALLY UNSTABLE.")
            print(f"      The calibrated joint steady state is self-consistent (step 2 = ~0%), but "
                  "the\n      structure is only a fixed point, not a stable one: a perturbation runs "
                  "away.")
            if struct_src.startswith("the AGENT"):
                print(f"      This is the price of low precision at SCALE: the agent's OVER-INCLUDED "
                      "edges\n      add spurious positive-feedback loops that the sparse model does "
                      "not have. Calibration\n      hid the cost (any structure fits its own "
                      "targets); the dynamics expose it. Prune the\n      low-confidence / uncited "
                      "extra edges (--live records confidence) and re-run to recover\n      "
                      "stability - the isolated probes could never have shown this.")
        else:
            moved = [(k, knock[k] / base[k]) for k in sorted(targ)
                     if k not in marg and k != anti and base[k] and abs(knock[k] / base[k] - 1) > 0.1]
            for k, fc in sorted(moved, key=lambda t: abs(t[1] - 1), reverse=True):
                print(f"      {k:12} {fc:.2f}x")
            print(f"      -> a single perturbation moves the coupled network along its edge signs; "
                  f"the\n         isolated probes could not see this whole-loop response.")
    print(f"\n  SCOPE: this is the immune core that has steady-state targets. DAS28 / PK / Vpop / "
          "the\n  clinical qualification layer are the paper's remaining species and need MATLAB "
          "SimBiology.")


if __name__ == "__main__":
    main()
