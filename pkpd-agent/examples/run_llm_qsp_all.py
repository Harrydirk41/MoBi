r"""One-command Stage-1 suite on ANY QSP model: derive keys, run every benchmark, report.

The orchestrator for Problem A/B: load a QSPModel from network.json + a spec, then run all
Stage-1 benchmarks (scope, topology, signs, readout, parameters, and - if the spec carries
a GSA list - sensitivity) each --repeat times, and print one consolidated ladder. Adding a
model = adding a QSPModelSpec to qsp_model.SPECS.

    python -m examples.run_llm_qsp_all --network network.json --model ra --repeat 3
    python -m examples.run_llm_qsp_all --network network.json --model ra --show-key

Needs ANTHROPIC_API_KEY (no live MATLAB). This is LONG: 6 benchmarks x repeat sessions.
"""

from __future__ import annotations

import argparse

from pkpd_agent.config import AgentConfig
from pkpd_agent.engines.qsp_model import QSPModel, get_spec
from pkpd_agent.llm import LLMPolicy
from pkpd_agent.loop import DecisionLoop
from pkpd_agent.state import ModelingSession
from pkpd_agent.tools.registry import ToolRegistry
from pkpd_agent.tools.qsp_topology_loop_tools import register_qsp_topology_loop_tools
from pkpd_agent.tools import qsp_loop_tools as Q


def _mean(xs):
    return sum(xs) / len(xs) if xs else float("nan")


def _bench_defs(model):
    P = model.spec.name
    defs = [
        ("scope", Q.register_qsp_scope_loop_tools, "qsp_scope_final",
         lambda f: ("F1", f["f1"]),
         f"You are scoping the {P} QSP model. Propose the cast (cells + mediators) it "
         "should include; be parsimonious (precision is scored). Use scope_propose then "
         "scope_finalize once."),
        ("topology", register_qsp_topology_loop_tools, "qsp_topo_final",
         lambda f: ("F1", f["topology"]["f1"]),
         f"Reconstruct the {P} regulatory network: propose signed edges among the nodes "
         "with network_propose (batches), then network_finalize once."),
        ("signs", Q.register_qsp_sign_loop_tools, "qsp_sign_final",
         lambda f: ("acc", f["accuracy"]),
         f"For each real edge of the {P} model, decide activate (+1) or inhibit (-1) with "
         "sign_predict, then sign_finalize once. Most edges activate; the inhibitory ones "
         "are anti-inflammatory mediators and negative self-feedback."),
        ("readout", Q.register_qsp_readout_loop_tools, "qsp_readout_final",
         lambda f: ("F1", f["f1"]),
         f"Which {P} nodes is the clinical readout computed from? Propose them with "
         "readout_propose, then readout_finalize once."),
        ("params", Q.register_qsp_params_loop_tools, "qsp_params_final",
         lambda f: ("phys-median", f["physiological"].get("median_log10_err")),
         f"Estimate every {P} parameter's value from its name and units (order-of-"
         "magnitude). Use param_estimate (batches) then param_finalize once."),
    ]
    if model.spec.gsa_top:
        defs.append(
            ("sensitivity", Q.register_qsp_sensitivity_loop_tools, "qsp_sens_final",
             lambda f: ("recall", f["recall"]),
             f"From the pool, rank the {P} parameters that most drive the readout with "
             "sens_rank, then sens_finalize once."))
    return defs


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--network", required=True)
    ap.add_argument("--model", default="ra", help="spec name in SPECS")
    ap.add_argument("--infer", action="store_true",
                    help="heuristically infer the spec from the dump (no hand config); "
                         "skips sensitivity (needs the external GSA list)")
    ap.add_argument("--repeat", type=int, default=1)
    ap.add_argument("--max-steps", type=int, default=20)
    ap.add_argument("--show-key", action="store_true")
    ap.add_argument("--only", help="comma-list of benchmarks to run (default all)")
    args = ap.parse_args()

    model = (QSPModel.inferred(args.network, "auto-inferred") if args.infer
             else QSPModel.from_network_json(args.network, get_spec(args.model)))
    print(f"== model '{model.spec.name}'"
          f"{' [SPEC AUTO-INFERRED]' if args.infer else ''} "
          f"(all keys DERIVED from {args.network}) ==")
    print(f"   {len(model.nodes)} nodes, {len(model.edges)} edges, {len(model.params)} "
          f"params, {len(model.readout_drivers)} readout drivers, "
          f"{len(model.spec.gsa_top)} GSA params\n")
    if args.show_key:
        print("nodes:", model.nodes)
        print("readout drivers:", model.readout_drivers)
        return

    cfg = AgentConfig(mock=False, max_steps=args.max_steps)
    if not cfg.anthropic_key_present():
        print("ANTHROPIC_API_KEY not set.")
        return

    only = set(args.only.split(",")) if args.only else None
    report = []
    for name, reg, key, extract, prompt in _bench_defs(model):
        if only and name not in only:
            continue
        print(f"\n---- {name} ({args.repeat}x) ----", flush=True)
        vals = []
        label = "?"
        for i in range(args.repeat):
            registry = ToolRegistry()
            reg(registry, cfg, {"model": model})
            policy = LLMPolicy(cfg, registry, prompt)
            loop = DecisionLoop(config=cfg, registry=registry, policy=policy)
            session = loop.run(f"Run the {name} benchmark, then finalize.",
                               ModelingSession(goal=name))
            final = session.get(key)
            if final:
                label, v = extract(final)
                if v is not None:
                    vals.append(v)
                    print(f"  run {i + 1}: {label} {v}", flush=True)
        if vals:
            report.append((name, label, _mean(vals), min(vals), max(vals), len(vals)))

    print("\n\n================ STAGE-1 LADDER: " + model.spec.name + " ================")
    print(f"{'benchmark':14s} {'metric':12s} {'mean':>7s} {'min':>7s} {'max':>7s} {'n':>3s}")
    for name, label, mean, lo, hi, n in report:
        print(f"{name:14s} {label:12s} {mean:7.3f} {lo:7.3f} {hi:7.3f} {n:3d}")
    print("(all answer keys were derived from the model dump, not hardcoded.)")


if __name__ == "__main__":
    main()
