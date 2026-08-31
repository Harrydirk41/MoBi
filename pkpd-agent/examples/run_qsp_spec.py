r"""Build the layered, provenance-tagged model spec for a node - one decision per layer.

Fills a single `model_spec` by running the pipeline in dependency order: frame (human tier, loaded
from tasks.json + a --node target), then the agent's topology and rate-law form, then constants
tagged by provenance. Prints the spec layer by layer and a provenance rollup (the honest
four-quadrant view), then optionally emits the SBML the spec assembles to.

    python -m examples.run_qsp_spec --model ra --node IL6            # offline all-candidates
    python -m examples.run_qsp_spec --model ra --node IL6 --live     # agent proposes topology+form
    python -m examples.run_qsp_spec --model ra --node IL6 --sbml out.xml

Nothing IL-6-specific; --node picks the target, everything else derives from the model + config.
"""

from __future__ import annotations

import argparse
import json
import os

from pkpd_agent.engines import model_assembly as MA, model_spec as MS
from pkpd_agent.engines.sbml_import import sbml_to_network
from examples.run_qsp_build_general import load_model, discover_nodes, _project_dir


def _frame(model, node):
    """Human tier, loaded from tasks.json where available (objective/acceptance live there);
    scope/scale are a proposed default here (the 'agent proposes, human approves' step)."""
    tj = os.path.join(_project_dir(model), "tasks.json")
    t = json.load(open(tj)) if os.path.isfile(tj) else {}
    return {"objective": t.get("trial_objective") or t.get("readout_desc"),
            "acceptance": t.get("fit_target"),
            "scope": {"target": node, "drugs": list((t.get("drugs") or {}).keys())},
            "scale": "single-compartment; cytokine concentrations; phenomenological Hill"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="ra")
    ap.add_argument("--node", default=None)
    ap.add_argument("--live", action="store_true")
    ap.add_argument("--sbml", default=None, help="also emit the assembled SBML to this path")
    args = ap.parse_args()

    prov, levels, cytokines = load_model(args.model)
    nodes = discover_nodes(prov, levels, cytokines)
    node = args.node or sorted(nodes)[0]
    if node not in nodes:
        print(f"'{node}' not buildable; discovered: {sorted(nodes)}"); return
    truth = set(nodes[node]["regulators"])
    candidates = sorted(c for c in cytokines if c != node)

    # ── FRAME (human tier) ──
    frame = _frame(args.model, node)
    print("== FRAME (human: objective+acceptance; approve: scope+scale) ==")
    print(f"  objective : {frame['objective']}")
    print(f"  acceptance: {frame['acceptance']}")
    print(f"  scale     : {frame['scale']}")

    # ── layer 2+3: the agent's decisions (or offline baseline) ──
    if args.live:
        from pkpd_agent.config import AgentConfig
        from pkpd_agent.engines import llm_tasks as LT
        cfg = AgentConfig(mock=False)
        if cfg.anthropic_key_present():
            call = LT.default_call(cfg)
            chosen = MA.propose_regulators(node, candidates, "secretion", call)
            motif = MA.propose_motif(node, [{"species": r["cytokine"]} for r in chosen], "", call)
            src = "LIVE agent"
        else:
            chosen = [{"cytokine": c} for c in candidates]
            motif = {"proliferation_order": "zeroth", "combination": "product", "cap": None}
            src = "offline (no key)"
    else:
        chosen = [{"cytokine": c} for c in candidates]
        motif = {"proliferation_order": "zeroth", "combination": "product", "cap": None}
        src = "offline all-candidates"
    print(f"  decisions from: {src}")

    spec = MS.build_spec(frame, node, chosen, motif, prov, levels, truth_regulators=truth)

    # ── print the spec layer by layer ──
    print(f"\n== ① NODES ({len(spec['nodes'])}) ==")
    print(f"  dynamic: {node}; inputs: {[n['id'] for n in spec['nodes'] if n['role']=='input']}")
    print(f"\n== ② EDGES / topology ({len(spec['edges'])}) ==")
    for e in spec["edges"]:
        mark = "✓in-model" if e["verify"]["in_model"] else "✗not-in-model"
        print(f"  {e['src']:6} --{e['sign']:4}--> {e['dst']:5} [{e['source']}] {mark}"
              f"  {e['basis'] or ''}")
    print(f"\n== ③ FORM / rate-law ==")
    f = spec["forms"][0]
    print(f"  {f['order']}/{f['combination']} cap={f['cap']}  [{f['source']}]  "
          f"stable={f['verify']['steady_state_stable']} responds={f['verify']['responds_to_single_knockdown']}")
    if f["reason"]:
        print(f"  reason: {f['reason']}")
    print(f"\n== ④ CONSTANTS ({len(spec['constants'])}) ==")
    for c in spec["constants"]:
        tag = c["provenance"]
        extra = []
        if c.get("identifiable") is False:
            extra.append("UNIDENTIFIABLE")
        if c.get("needs_data"):
            extra.append(f"needs {c['needs_data']}")
        print(f"  {c['param']:10} = {str(c['value']):8} [{tag}] {' '.join(extra)}"
              f"{'  ref:'+c['ref'] if c.get('ref') else ''}")

    roll = MS.provenance_rollup(spec)
    print(f"\n== provenance rollup (honest four-quadrant) ==")
    print(f"  constants by source: {roll['constants_by_source']}")
    print(f"  edges: {roll['edges_matching_model']}/{roll['edges_total']} match the model")
    print(f"  must-fit / needs-data: {roll['needs_data']}")

    if args.sbml:
        sub = MS.to_subsystem(spec, prov, levels)
        open(args.sbml, "w", encoding="utf-8").write(MA.to_sbml(sub))
        net = sbml_to_network(args.sbml)
        print(f"\n== assembled SBML -> {args.sbml} "
              f"({len(net['species'])} species, {len(net['reactions'])} reactions) ==")


if __name__ == "__main__":
    main()
