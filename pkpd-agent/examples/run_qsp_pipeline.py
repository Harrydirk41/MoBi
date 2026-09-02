r"""Decomposable build: choose, per step, HOW it is filled - human/read/LLM/data - no scaffolding.

Each modelling layer is a pluggable provider. Flags pick each layer's MODE independently:

    --frame     given | llm          (objective+acceptance human; scope+scale given or agent-proposed)
    --target    given:<node> | llm   (which node to build: human names it, or the agent picks)
    --topology  llm | given:a,b,c | data   (agent proposes edges / human lists / read model structure)
    --form      llm | given          (agent picks the rate-law motif, or a stated default)

Examples:

    # fully human-seeded except topology+form from the agent:
    python -m examples.run_qsp_pipeline --model ra --target given:IL6 --topology llm --form llm

    # no LLM at all (human/data baseline), runs anywhere:
    python -m examples.run_qsp_pipeline --model ra --target given:IL6 --topology data --form given

The model's own structure is just the `data` topology provider - one explicit choice, used as a
baseline or for scoring, never a hidden default. Emits the assembled SBML and a provenance rollup.
"""

from __future__ import annotations

import argparse
import json
import os

from pkpd_agent.engines import model_assembly as MA, model_spec as MS, pipeline as P
from pkpd_agent.engines.sbml_import import sbml_to_network
from examples.run_qsp_build_general import load_model, discover_nodes, _project_dir

_DEFAULT_MOTIF = {"proliferation_order": "zeroth", "combination": "product", "cap": None}


def frame_from_config(ctx):
    t = ctx["tasks"]
    return {"objective": t.get("trial_objective") or t.get("readout_desc"),
            "acceptance": t.get("fit_target"),
            "scope": {"drugs": list((t.get("drugs") or {}).keys())},
            "scale": "single-compartment; concentrations; phenomenological Hill"}


def parse_step(flag, default_mode):
    """'given:IL6' -> (given, 'IL6'); 'llm' -> (llm, None); 'data' -> (data, None)."""
    if flag is None:
        return default_mode, None
    if ":" in flag:
        m, v = flag.split(":", 1)
        return m, v
    return flag, None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="ra")
    ap.add_argument("--frame", default="given")
    ap.add_argument("--target", default=None, help="given:<node> | llm (default: first node)")
    ap.add_argument("--topology", default="data", help="llm | given:a,b,c | data")
    ap.add_argument("--form", default="given", help="llm | given")
    ap.add_argument("--sbml", default=None)
    args = ap.parse_args()

    prov, levels, cytokines = load_model(args.model)
    nodes = discover_nodes(prov, levels, cytokines)
    tasks = json.load(open(os.path.join(_project_dir(args.model), "tasks.json")))

    # LLM boundary (only built if any step needs it)
    call = None
    need_llm = "llm" in (args.frame, args.target, args.topology, args.form)
    if need_llm:
        from pkpd_agent.config import AgentConfig
        from pkpd_agent.engines import llm_tasks as LT
        cfg = AgentConfig(mock=False)
        if cfg.anthropic_key_present():
            call = LT.default_call(cfg)
        else:
            print("[no ANTHROPIC_API_KEY: llm steps fall back to given/data]")

    # ── resolve the target first (topology/candidates depend on it) ──
    tmode, tval = parse_step(args.target, "given")
    if tmode == "llm" and call:
        sysp = "Pick ONE cytokine node to model from this disease's node list. JSON {\"node\": name}."
        target_provider = P.from_llm(sysp, lambda c: "Nodes: " + ", ".join(sorted(nodes)),
                                     lambda s: __import__("json").loads(s).get("node"))
    else:
        target_provider = P.given(tval or sorted(nodes)[0])
    target = target_provider.fn({})
    if target not in nodes:
        print(f"'{target}' not buildable; discovered {sorted(nodes)}"); return
    candidates = sorted(c for c in cytokines if c != target)
    truth = set(nodes[target]["regulators"])

    # ── frame provider ──
    fmode, _ = parse_step(args.frame, "given")
    if fmode == "llm" and call:
        frame_provider = P.from_llm(
            "Propose the scope+scale for this QSP node model. JSON {\"scale\": str, \"scope\": {}}.",
            lambda c: f"Objective: {tasks.get('trial_objective')}. Node: {target}.",
            lambda s: {**frame_from_config({"tasks": tasks}), **__import__("json").loads(s)})
    else:
        frame_provider = P.given(frame_from_config({"tasks": tasks}))

    # ── topology provider ──
    top_mode, top_val = parse_step(args.topology, "data")
    if top_mode == "llm" and call:
        topology_provider = P.from_llm(
            MA._REG_SYS,
            lambda c: (f"Cell: {target}. Process: secretion. Available cytokine nodes: "
                       + ", ".join(candidates) +
                       '\n\nWhich regulate it? JSON {"regulators":[{"cytokine":n,"direction":'
                       '"up"|"down","basis":"one phrase"}]}. Only from the list.'),
            lambda s: [r for r in (MA._parse_json(s).get("regulators") or [])
                       if isinstance(r, dict) and r.get("cytokine") in set(candidates)])
    elif top_mode == "given" and top_val:
        chosen = [c.strip() for c in top_val.split(",") if c.strip() in candidates]
        topology_provider = P.given([{"cytokine": c, "direction": "up"} for c in chosen])
    else:                                                  # data: read the model's own regulators
        topology_provider = P.data(lambda c: [{"cytokine": r, "direction": "up"}
                                              for r in nodes[target]["regulators"]])

    # ── form provider ──
    form_mode, _ = parse_step(args.form, "given")
    if form_mode == "llm" and call:
        form_provider = P.from_llm(
            MA._MOTIF_SYS if hasattr(MA, "_MOTIF_SYS") else "Choose the rate-law form. JSON only.",
            lambda c: f"Cell: {target}. Regulators modulate its secretion. "
                      '\nReturn JSON {"proliferation_order":"zeroth"|"first","combination":'
                      '"product"|"capped_sum","cap":number|null,"per_regulator":"hill","reason":"..."}.',
            lambda s: {**_DEFAULT_MOTIF, **MA._parse_json(s)})
    else:
        form_provider = P.given(dict(_DEFAULT_MOTIF))

    providers = {"frame": frame_provider, "target": target_provider,
                 "topology": topology_provider, "form": form_provider}
    ctx = {"prov": prov, "levels": levels, "truth": truth, "call": call, "tasks": tasks}
    spec = P.run(providers, ctx)

    # ── report: each layer + HOW it was filled ──
    print(f"== build of '{target}': per-layer mode ==")
    for layer, mode in spec["modes"].items():
        print(f"  {layer:9} <- {mode}")
    print(f"\n② topology ({len(spec['edges'])} edges, [{spec['modes']['topology']}]):")
    for e in spec["edges"]:
        print(f"    {e['src']:6} {e['sign']:4} {'✓' if e['verify']['in_model'] else '✗'} "
              f"{e['basis'] or ''}")
    f = spec["forms"][0]
    print(f"\n③ form [{spec['modes']['form']}]: {f['order']}/{f['combination']} "
          f"responds={f['verify']['responds_to_single_knockdown']}")
    roll = MS.provenance_rollup(spec)
    print(f"\n④ constants: {roll['constants_by_source']}")
    print(f"   needs-data params: {len(roll['needs_data'])}")
    print(f"   edges matching model: {roll['edges_matching_model']}/{roll['edges_total']}")

    if args.sbml:
        sub = MS.to_subsystem(spec, prov, levels)
        open(args.sbml, "w", encoding="utf-8").write(MA.to_sbml(sub))
        net = sbml_to_network(args.sbml)
        print(f"\nassembled SBML -> {args.sbml} ({len(net['species'])} species, "
              f"{len(net['reactions'])} reactions)")


if __name__ == "__main__":
    main()
