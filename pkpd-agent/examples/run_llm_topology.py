r"""Stage-1 model building, isolated to its hard core: given the model's NODES and the
paper's references, can an LLM reconstruct the WIRING (which node influences which)?

The mature-modeller workflow is library + assembly, and the creative step is topology -
deciding the edges. This benchmark hands the LLM the answer's node list (the easy part -
species are enumerable) and the literature, asks it to propose the signed influence edges,
and scores that draft against the edges extracted from the model itself (the answer key).

    python -m examples.run_llm_topology --model ra ^
        --sbproj "..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj" ^
        --refs references.txt

--refs is a plain-text file of the paper's references (start with the original authors' own,
then raise the difficulty). Without it the LLM works from the node names alone - a useful
floor showing how far bare biological priors get you.

Needs ANTHROPIC_API_KEY and the MATLAB engine (to dump network.json, the answer key).
"""

from __future__ import annotations

import argparse
import os

from pkpd_agent.config import AgentConfig
from pkpd_agent.engines.simbiology import SimBiologyEngine
from pkpd_agent.engines.qsp_model import QSPModel, get_spec
from pkpd_agent.engines import llm_tasks as LT, llm_topology as TOP


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="ra")
    ap.add_argument("--sbproj", required=True)
    ap.add_argument("--refs", default=None, help="text file of the paper's references")
    ap.add_argument("--out", default="network.json", help="where to dump the answer key")
    ap.add_argument("--llm-model", default=None)
    args = ap.parse_args()

    cfg_llm = AgentConfig(mock=False)
    if args.llm_model:
        cfg_llm.model = args.llm_model
    if not cfg_llm.anthropic_key_present():
        print("ANTHROPIC_API_KEY not set.")
        return

    references = ""
    if args.refs:
        with open(args.refs, encoding="utf-8") as fh:
            references = fh.read().strip()
        print(f"loaded {len(references)} chars of references from {args.refs}")
    else:
        print("no --refs given: LLM works from node names + biological priors alone")

    sb = SimBiologyEngine()
    try:
        print("== starting MATLAB engine =="); sb.start()
        print(f"== loading {os.path.basename(args.sbproj)} =="); sb.load_project(args.sbproj)
        print(f"== dumping network structure -> {args.out} (the answer key) ==", flush=True)
        network = sb.network_json(args.out)
    finally:
        sb.stop()

    spec = get_spec(args.model)
    model = QSPModel(network, spec)
    nodes = list(model.nodes)                         # biological species = the given NODES
    print(f"model has {len(network.get('species', []))} species; "
          f"{len(nodes)} are biological nodes (drug/PK/readout excluded)")

    # answer key: influence edges, restricted to pairs BOTH among the given nodes (a fair
    # comparison - the LLM is only ever given the biological node list).
    all_truth = TOP.ground_truth_edges(network)
    nset = set(nodes)
    truth = {(s, d) for s, d in all_truth if s in nset and d in nset}
    print(f"answer key: {len(truth)} influence edges among the biological nodes "
          f"({len(all_truth)} total including drug/readout)")

    print("== LLM drafting the topology from nodes + references ==", flush=True)
    draft = TOP.draft_topology(nodes, references, LT.default_call(cfg_llm))
    print(f"LLM proposed {len(draft)} edges")

    r = TOP.compare_topology(draft, truth)
    print(f"\n== topology reconstruction ==")
    print(f"  precision {r['precision']}   recall {r['recall']}   f1 {r['f1']}")
    print(f"  {r['hit']} hit / {r['n_draft']} drafted / {r['n_truth']} in answer key")

    if r["missed"]:
        print(f"\n  MISSED ({len(r['missed'])}) - in the model, the LLM did not draw:")
        for s, d in r["missed"][:40]:
            print(f"    {s} -> {d}")
        if len(r["missed"]) > 40:
            print(f"    ... and {len(r['missed']) - 40} more")
    if r["extra"]:
        print(f"\n  EXTRA ({len(r['extra'])}) - the LLM drew, not in the model "
              "(literature-plausible but not wired here):")
        for s, d in r["extra"][:40]:
            print(f"    {s} -> {d}")
        if len(r["extra"]) > 40:
            print(f"    ... and {len(r['extra']) - 40} more")


if __name__ == "__main__":
    main()
