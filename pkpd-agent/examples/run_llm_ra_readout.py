r"""Stage-3: can the LLM reconstruct the mechanism->endpoint mapping (DAS28-CRP)?

The agent proposes which model nodes the DAS28-CRP readout is built from; scored against
the species the model's readout rule actually depends on (extracted from network.json by
walking the algebraic rule graph). Needs the full network.json (from dump_network.py).

    python -m examples.run_llm_ra_readout --network network.json
    python -m examples.run_llm_ra_readout --network network.json --repeat 5
    python -m examples.run_llm_ra_readout --network network.json --show-key   # just print the answer key
"""

from __future__ import annotations

import argparse

from pkpd_agent.config import AgentConfig
from pkpd_agent.engines import ra_readout as RO
from pkpd_agent.llm import LLMPolicy
from pkpd_agent.loop import DecisionLoop
from pkpd_agent.state import Decision, Finish, ModelingSession, Observation
from pkpd_agent.tools.registry import ToolRegistry
from pkpd_agent.tools.ra_readout_loop_tools import register_ra_readout_loop_tools


def _system_prompt() -> str:
    return (
        "You are a QSP modeler defining a clinical readout. The model computes DAS28-CRP - "
        "a composite disease-severity score - from its physiological state. Propose which "
        "model nodes (cell densities and mediators) that readout is a direct function of. "
        "Reason from what the score physically reflects: synovial cell load (joint "
        "swelling/tenderness) and the systemic acute-phase signal (CRP, driven by a "
        "specific cytokine). It is a readout formula, not the whole network - be focused. "
        "Propose with readout_propose, then readout_finalize once."
    )


def _report_variance(tag: str, xs: list) -> None:
    n = len(xs)
    mean = sum(xs) / n
    sd = (sum((x - mean) ** 2 for x in xs) / (n - 1)) ** 0.5 if n > 1 else 0.0
    print(f"\n  {tag:16s} over {n} runs: mean {mean:.3f}  sd {sd:.3f}  "
          f"min {min(xs):.3f}  max {max(xs):.3f}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--network", required=True, help="network.json (full rule graph)")
    ap.add_argument("--show-key", action="store_true",
                    help="print the extracted readout drivers + raw rules and exit")
    ap.add_argument("--max-steps", type=int, default=10)
    ap.add_argument("--repeat", type=int, default=1)
    ap.add_argument("--model", default=None)
    ap.add_argument("--effort", default=None)
    args = ap.parse_args()

    ext = RO.readout_drivers(args.network)
    drivers = ext["drivers"]
    print(f"readout targets found: {ext['targets_found']}")
    print(f"readout drivers (answer key, {len(drivers)}): {drivers}\n")
    if args.show_key or not drivers:
        for t, r in ext["target_rules"].items():
            print(f"  rule {t} = {r[:400]}")
        if not drivers:
            print("\n  !! no readout drivers extracted - DAS28_CRP is a state variable, "
                  "not an algebraic rule, so the rule-graph walk cannot reach the biology.")
            print("  Find its real definition with:")
            print("    python -c \"import json; d=json.load(open('network.json')); "
                  "[print('RULE',u.get('rule')) for u in d['rules'] if 'DAS28' in "
                  "u.get('rule','')]; [print('RXN',r['name'],'|',r.get('reaction'),'|',"
                  "r.get('rate')) for r in d['reactions'] if 'DAS28' in "
                  "(str(r.get('reaction'))+str(r.get('rate')))]\"")
        return

    cfg = AgentConfig(mock=False, max_steps=args.max_steps)
    if args.model:
        cfg.model = args.model
    if args.effort:
        cfg.effort = args.effort
    if not cfg.anthropic_key_present():
        print("ANTHROPIC_API_KEY not set.")
        return

    registry = ToolRegistry()
    register_ra_readout_loop_tools(registry, cfg, {"drivers": drivers})
    goal = ("Reconstruct the DAS28-CRP readout mapping: propose the model nodes it is "
            "computed from, then finalize to be scored.")
    verbose = args.repeat == 1

    def show(ev):
        if not verbose:
            return
        if isinstance(ev, Decision):
            if ev.text:
                print(f"\n[reason] {ev.text[:1200]}")
            for c in ev.calls:
                n = len((c.arguments or {}).get("nodes", [])) if c.name == "readout_propose" else ""
                print(f"  -> {c.name} {('('+str(n)+' nodes)') if n != '' else ''}")
        elif isinstance(ev, Observation):
            print(f"  <- {ev.tool}: {ev.content.get('message', '')}")
        elif isinstance(ev, Finish):
            print(f"\n=== AGENT SUMMARY ===\n{ev.text}")

    runs = []
    for i in range(args.repeat):
        policy = LLMPolicy(cfg, registry, _system_prompt())
        loop = DecisionLoop(config=cfg, registry=registry, policy=policy)
        session = loop.run(goal, ModelingSession(goal=goal), on_event=show)
        final = session.get("readout_final")
        if not final:
            print(f"  run {i + 1}: agent did not finalize")
            continue
        runs.append(final)
        if not verbose:
            print(f"  run {i + 1}/{args.repeat}: F1 {final['f1']} "
                  f"(P {final['precision']} R {final['recall']}), {final['hit']}/"
                  f"{final['n_truth']}")

    if not runs:
        print("\nno finalized runs.")
        return

    if verbose:
        f = runs[0]
        print("\n== READOUT mapping vs the model ==")
        print(f"  F1 {f['f1']}   precision {f['precision']}   recall {f['recall']}")
        print(f"  {f['hit']}/{f['n_truth']} real drivers recovered")
        print(f"\n  missed (real drivers): {f['missed']}")
        print(f"  extra (not in the readout rule): {f['extra']}")
    else:
        _report_variance("readout F1", [r["f1"] for r in runs])
        _report_variance("recall", [r["recall"] for r in runs])


if __name__ == "__main__":
    main()
