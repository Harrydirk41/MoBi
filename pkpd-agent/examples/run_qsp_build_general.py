r"""General from-scratch build - NO project-specific scaffolding, any node, any model.

Everything the earlier end-to-end driver hardcoded for the IL-6 hub is derived here from the model
and its config instead:

  * WHICH nodes are buildable            - discovered from the model's own naming conventions
                                           (a secretion `<Node>Sec...Maxby<Cyt>` set + a `kcl_<Node>`
                                           clearance + a steady-state target). Nothing says "IL6".
  * WHICH node to build                  - --node, or the first discovered; the candidate list is
                                           printed so an agent/user picks from the model, not a
                                           baked-in choice.
  * The REGULATOR CANDIDATES             - every cytokine node in the target file (general).
  * The STRUCTURE                        - the agent proposes regulators (--live); offline, an
                                           over-inclusive all-candidates baseline (no cherry-picked
                                           recorded answer).
  * The STRENGTHS                        - looked up by the model's `<Node>Sec...Maxby<Cyt>`
                                           convention (config data, node-agnostic).
  * The HELD-OUT EXPERIMENTS             - auto-generated, one knockdown per chosen regulator.

Run it on any node to see there is no IL-6 scaffolding left:

    python -m examples.run_qsp_build_general --model ra --node TNFa
    python -m examples.run_qsp_build_general --model ra --node IL17 --live
    python -m examples.run_qsp_build_general --model ra --list      # show buildable nodes

Pure (offline baseline) / needs a key only for --live. Emits a real SBML file per node.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import tempfile

from pkpd_agent.engines import model_assembly as MA
from pkpd_agent.engines.sbml_import import sbml_to_network
from examples.run_qsp_end_to_end import integrate

_SEC = re.compile(r"(?i)([A-Za-z0-9]+)Sec[A-Za-z0-9]*_Maxby([A-Za-z0-9]+)")


def _project_dir(model):
    """Resolve --model (a name or alias) to its projects/<dir>, via project_aliases in tasks.json."""
    base = os.path.join(os.path.dirname(__file__), "..", "projects")
    if os.path.isdir(os.path.join(base, model, "data")):
        return os.path.join(base, model)
    for d in sorted(os.listdir(base)):
        tj = os.path.join(base, d, "tasks.json")
        if os.path.isfile(tj):
            try:
                aliases = json.load(open(tj)).get("project_aliases") or []
            except Exception:                              # noqa: BLE001
                aliases = []
            if model == d or model in aliases:
                return os.path.join(base, d)
    raise SystemExit(f"no project found for --model '{model}' under {base}")


def load_model(model):
    root = os.path.join(_project_dir(model), "data")
    prov = {p["name"]: p for p in json.load(open(os.path.join(root, "param_provenance.json")))}
    tg = {t["model_species"]: t for t in json.load(open(os.path.join(root,
          "steady_state_targets.json"))) if t.get("model_species")}
    levels = {c: float(tg[c]["target_model_unit"]) for c in tg
              if tg[c].get("target_model_unit") is not None}
    cytokines = {c for c in levels if tg[c].get("kind") == "cytokine"}
    return prov, levels, cytokines


def discover_nodes(prov, levels, cytokines):
    """A node is buildable if the model exposes a secretion-regulator set for it, a clearance
    parameter, and a steady-state target - all read from the model's own conventions."""
    regs = {}
    for n in prov:
        m = _SEC.search(n)
        if m:
            regs.setdefault(m.group(1), set()).add(m.group(2))
    nodes = {}
    for node, rs in regs.items():
        clr = next((n for n in prov if re.match(rf"(?i)(kcl|kd)_{node}\b", n)), None)
        if node in levels and clr and node in cytokines:
            nodes[node] = {"regulators": sorted(r for r in rs if r in levels and r != node),
                           "clearance": clr, "level": levels[node]}
    return nodes


def lookup_max(prov, node, cyt):
    for n, p in prov.items():
        m = _SEC.search(n)
        if m and m.group(1) == node and m.group(2) == cyt and p.get("value_from_reference") is not None:
            return float(p["value_from_reference"])
    return None


def assemble_node(node, chosen, prov, levels, clearance, motif, clamp=None):
    """Node-agnostic assembly: same machinery as the IL-6 driver, parameterised by node name."""
    regs = [c for c in chosen if c in levels and c != node]
    reg_specs, values, from_data, from_prior, excess = [], {}, [], [], []
    for c in regs:
        mx = lookup_max(prov, node, c)
        if mx is not None:
            from_data.append(c)
        else:
            mx = 1.5; from_prior.append(c)                 # generic up-prior for an uncited edge
        reg_specs.append({"species": c, "max_param": f"M_{c}", "k_param": f"K_{c}"})
        values[f"M_{c}"] = mx; values[f"K_{c}"] = levels[c]
        excess.append((mx - 1.0) * 0.5)
    if motif.get("combination") == "capped_sum":
        s = sum(excess); cap = motif.get("cap")
        eff = min(cap, 1.0 + s) if cap not in (None, "", 0) else 1.0 + s
    else:
        eff = 1.0
        for e in excess:
            eff *= 1.0 + e
    kcl = float(prov[clearance]["value_from_reference"])
    kg = levels[node] * kcl / eff
    values[f"kg_{node}"] = kg; values[clearance] = kcl; values[f"{node}_init"] = 0.0
    if clamp is None:
        clamp = {c: levels[c] for c in regs}
    spec = MA.build_subsystem(node, f"kg_{node}", clearance, reg_specs, values, clamp=clamp,
                              motif=motif)
    return spec, from_data, from_prior, kcl


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="ra")
    ap.add_argument("--node", default=None, help="which node to build; default = first discovered")
    ap.add_argument("--list", action="store_true", help="list buildable nodes and exit")
    ap.add_argument("--live", action="store_true", help="agent proposes the structure via the LLM")
    args = ap.parse_args()

    prov, levels, cytokines = load_model(args.model)
    nodes = discover_nodes(prov, levels, cytokines)

    if args.list or not nodes:
        print(f"== buildable nodes discovered in '{args.model}' (no hardcoding) ==")
        for n, info in sorted(nodes.items()):
            print(f"  {n:7} regulators={info['regulators']}  clearance={info['clearance']}")
        return

    node = args.node or sorted(nodes)[0]
    if node not in nodes:
        print(f"'{node}' is not buildable; discovered: {sorted(nodes)}"); return
    info = nodes[node]
    truth = set(info["regulators"])
    candidates = sorted(c for c in cytokines if c != node)     # general candidate list

    print(f"== building node '{node}' (target {info['level']:g}, clearance {info['clearance']}) ==")
    print(f"  candidate regulators (all cytokines): {candidates}")

    # ---- STRUCTURE: agent proposes, or an over-inclusive all-candidates baseline ----
    if args.live:
        from pkpd_agent.config import AgentConfig
        from pkpd_agent.engines import llm_tasks as LT
        cfg = AgentConfig(mock=False)
        if not cfg.anthropic_key_present():
            print("  --live but no key; using the all-candidates baseline.")
            chosen = candidates; motif = {"proliferation_order": "zeroth",
                                          "combination": "product", "cap": None}
        else:
            call = LT.default_call(cfg)
            regs = MA.propose_regulators(node, candidates, "secretion", call)
            chosen = [r["cytokine"] for r in regs]
            motif = MA.propose_motif(node, [{"species": c} for c in chosen], "", call)
            print(f"  (LIVE) agent chose: {chosen}")
            print(f"  (LIVE) rate-law: {motif.get('proliferation_order')}/{motif.get('combination')}")
    else:
        chosen = candidates                                    # no cherry-picked recorded answer
        motif = {"proliferation_order": "zeroth", "combination": "product", "cap": None}
        print("  (no agent: over-inclusive all-candidates baseline; --live for the agent's choice)")

    hit = truth & set(chosen)
    rec = len(hit) / len(truth) if truth else 0
    prec = len(hit) / len(chosen) if chosen else 0
    print(f"  structure vs model's own regulators {sorted(truth)}: recall {rec:.2f}, "
          f"precision {prec:.2f}")

    # ---- ASSEMBLE + FIT + EMIT ----
    spec, fd, fp, kcl = assemble_node(node, chosen, prov, levels, info["clearance"], motif)
    xml = os.path.join(tempfile.gettempdir(), f"build_{node}.xml")
    open(xml, "w", encoding="utf-8").write(MA.to_sbml(spec))
    net = sbml_to_network(xml)
    print(f"\n  assembled: {len(spec['species'])} species, {len(spec['reactions'])} reactions; "
          f"strengths {fd} from data, {fp} from prior")
    print(f"  emitted {xml}")

    # ---- SIMULATE: train + auto-generated held-out knockdowns ----
    clamp = {c: levels[c] for c in chosen if c in levels and c != node}
    ss = integrate(net, clamp)[node]
    print(f"\n  TRAIN steady state: {ss:.4g}  (target {info['level']:g}, "
          f"error {abs(ss-info['level'])/info['level']:.1%})")
    print(f"  HELD-OUT (auto: knock each chosen regulator to 10%):")
    for c in sorted(set(chosen) & set(clamp)):
        held = dict(clamp); held[c] = levels[c] * 0.1
        print(f"    anti-{c:6} -> {node} = {integrate(net, held)[node]:.4g}")
    print(f"\n  -> same code, node '{node}', zero IL-6 scaffolding: discovered the node, its "
          "regulators,\n     its clearance and target from the model; assembled, fit, emitted "
          "SBML, and simulated.")


if __name__ == "__main__":
    main()
