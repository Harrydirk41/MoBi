r"""LLM structure extraction + the RA regression that verifies it.

Reads a model's network.json, has an LLM extract the canonical structure (nodes / signed
edges / readout drivers) from the rule expressions - understanding whatever naming the
author used - and REGRESSES it against the deterministic parser on RA. The extractor must
reproduce the deterministic 88 edges / 9 readout drivers before we would trust it on a
model whose naming convention the regex does not understand.

    python -m examples.run_llm_extract --network network.json --model ra

Needs ANTHROPIC_API_KEY. This is the piece that turns 'derive by regex (fixed convention)'
into 'derive by LLM (any convention), validated against a known convention'.
"""

from __future__ import annotations

import argparse
import json

from pkpd_agent.config import AgentConfig
from pkpd_agent.engines import llm_structure as LS
from pkpd_agent.engines.qsp_model import QSPModel, get_spec


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--network", required=True)
    ap.add_argument("--model", default="vantage_ra", help="spec for the deterministic ground truth")
    ap.add_argument("--llm-model", default=None)
    args = ap.parse_args()

    with open(args.network, encoding="utf-8") as fh:
        data = json.load(fh)

    # deterministic ground truth (the convention we DO understand)
    truth = QSPModel(data, get_spec(args.model))
    print(f"deterministic parser: {len(truth.nodes)} nodes, {len(truth.edges)} edges, "
          f"{len(truth.readout_drivers)} readout drivers\n")

    cfg = AgentConfig(mock=False)
    if args.llm_model:
        cfg.model = args.llm_model
    if not cfg.anthropic_key_present():
        print("ANTHROPIC_API_KEY not set.")
        return

    print("== LLM reading the model structure (any naming) ==", flush=True)
    call = LS.default_call(cfg)
    result = LS.extract_structure(data, call,
                                  readout_hint=truth.spec.readout_targets)
    print(f"LLM extracted: {len(result['nodes'])} nodes, {len(result['edges'])} edges, "
          f"{len(result['readout_drivers'])} readout drivers\n")

    print("== REGRESSION vs the deterministic parser (must reproduce it) ==")
    cmp = LS.compare_to_model(result, truth)
    for k in ("nodes", "edges_topology", "edges_signed", "readout"):
        d = cmp[k]
        print(f"  {k:16s} P {d['precision']}  R {d['recall']}  F1 {d['f1']}  "
              f"({d['hit']}/{d['n_truth']})")
    print("\n  edges the LLM found but the parser didn't:",
          cmp["edge_disagreements"]["llm_only"])
    print("  edges the parser found but the LLM didn't:",
          cmp["edge_disagreements"]["deterministic_only"])
    print("\nRead: high edges/readout recall means the LLM extractor reproduces the known-"
          "convention answer key, so it can be trusted to read a NEW convention. Low recall "
          "means do NOT trust it as a key generator yet.")


if __name__ == "__main__":
    main()
