r"""Model-AGNOSTIC topology benchmark: point at ANY QSP model's network.json + a spec.

Proves Problem A - the same reconstruction benchmark, with the cast and answer key DERIVED
from the model dump (qsp_model.QSPModel) instead of hardcoded RA vocab. Add a new model by
adding a QSPModelSpec to qsp_model.SPECS; nothing else changes.

    python -m examples.run_llm_qsp_topology --network network.json --model ra
    python -m examples.run_llm_qsp_topology --network network.json --model ra --repeat 5
    python -m examples.run_llm_qsp_topology --network network.json --model ra --show-key
"""

from __future__ import annotations

import argparse

from pkpd_agent.config import AgentConfig
from pkpd_agent.engines.qsp_model import QSPModel, get_spec
from pkpd_agent.llm import LLMPolicy
from pkpd_agent.loop import DecisionLoop
from pkpd_agent.state import Decision, Finish, ModelingSession, Observation
from pkpd_agent.tools.registry import ToolRegistry
from pkpd_agent.tools.qsp_topology_loop_tools import register_qsp_topology_loop_tools


def _sys(model) -> str:
    return (
        f"You are reconstructing the regulatory network of the {model.spec.name} QSP "
        "model. You are given its nodes and must propose the directed, signed edges - "
        "which node up- or down-regulates which. Reason from the biology of the system. "
        "Build the draft in batches with network_propose, then network_finalize once. "
        "Balance recall against precision."
    )


def _var(tag, xs):
    n = len(xs); mean = sum(xs) / n
    sd = (sum((x - mean) ** 2 for x in xs) / (n - 1)) ** 0.5 if n > 1 else 0.0
    print(f"\n  {tag:14s} over {n} runs: mean {mean:.3f}  sd {sd:.3f}  "
          f"min {min(xs):.3f}  max {max(xs):.3f}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--network", required=True)
    ap.add_argument("--model", default="ra", help="spec name in qsp_model.SPECS")
    ap.add_argument("--show-key", action="store_true")
    ap.add_argument("--repeat", type=int, default=1)
    ap.add_argument("--max-steps", type=int, default=20)
    args = ap.parse_args()

    model = QSPModel.from_network_json(args.network, get_spec(args.model))
    print(f"model '{model.spec.name}': {len(model.nodes)} nodes, {len(model.edges)} "
          f"edges, {len(model.params)} params, {len(model.readout_drivers)} readout "
          f"drivers (all DERIVED from {args.network})\n")
    if args.show_key:
        print("nodes:", model.nodes)
        print("edges:", [(e.source, e.sign, e.target) for e in model.edges])
        print("readout drivers:", model.readout_drivers)
        return

    cfg = AgentConfig(mock=False, max_steps=args.max_steps)
    if not cfg.anthropic_key_present():
        print("ANTHROPIC_API_KEY not set.")
        return

    registry = ToolRegistry()
    register_qsp_topology_loop_tools(registry, cfg, {"model": model})
    goal = "Reconstruct the network; propose signed edges, then finalize."
    verbose = args.repeat == 1

    def show(ev):
        if not verbose:
            return
        if isinstance(ev, Decision):
            if ev.text:
                print(f"\n[reason] {ev.text[:900]}")
            for c in ev.calls:
                print(f"  -> {c.name}")
        elif isinstance(ev, Observation):
            print(f"  <- {ev.tool}: {ev.content.get('message', '')}")
        elif isinstance(ev, Finish):
            print(f"\n=== SUMMARY ===\n{ev.text[:1500]}")

    runs = []
    for i in range(args.repeat):
        policy = LLMPolicy(cfg, registry, _sys(model))
        loop = DecisionLoop(config=cfg, registry=registry, policy=policy)
        session = loop.run(goal, ModelingSession(goal=goal), on_event=show)
        final = session.get("qsp_topo_final")
        if final:
            runs.append(final)
            if not verbose:
                t = final["topology"]
                print(f"  run {i + 1}/{args.repeat}: topo F1 {t['f1']} "
                      f"(P {t['precision']} R {t['recall']})")

    if runs and not verbose:
        _var("TOPOLOGY F1", [r["topology"]["f1"] for r in runs])
    elif runs:
        t = runs[0]["topology"]
        print(f"\n== TOPOLOGY vs {model.spec.name} ==\n  F1 {t['f1']}  P {t['precision']}"
              f"  R {t['recall']}  ({t['hit']}/{t['n_truth']})")


if __name__ == "__main__":
    main()
