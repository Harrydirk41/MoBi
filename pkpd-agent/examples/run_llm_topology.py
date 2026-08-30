r"""Stage-1 model building, isolated to its hard core: given the model's NODES and the
paper's references, can an LLM reconstruct the WIRING (which node influences which)?

The mature-modeller workflow is library + assembly, and the creative step is topology -
deciding the edges. This benchmark hands the LLM the answer's node list (the easy part -
species are enumerable) plus the literature, asks it to propose the signed influence edges,
and scores that draft against the edges extracted from the model itself (the answer key).

Two ways to supply the literature (do both, compare):
  Route A (cached, reproducible): --paper reads the paper's OWN reference list; optional
      --fetch-abstracts enriches each citation with its PubMed abstract, cached to a file.
      The drafter is a plain tool-less call - it sees only the fixed, inspectable text.
  Route B (self-reading): --web gives the drafter the citation LIST and a web-search tool,
      so it goes and reads the references itself. Non-deterministic, but no pre-assembly.

    python -m examples.run_llm_topology --model ra ^
        --sbproj "..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj" ^
        --paper paper.txt --fetch-abstracts          (Route A, richest material)
    python -m examples.run_llm_topology --model ra ^
        --sbproj "...sbproj" --paper paper.txt --web  (Route B, LLM reads them itself)

--paper may be the full paper text (the References section is sliced out automatically) or
just the reference list. With neither --paper nor --refs the LLM works from node names alone
- a floor showing how far bare biological priors get you.

Needs ANTHROPIC_API_KEY and the MATLAB engine (to dump network.json, the answer key).
"""

from __future__ import annotations

import argparse
import os

from pkpd_agent.config import AgentConfig
from pkpd_agent.engines.simbiology import SimBiologyEngine
from pkpd_agent.engines.qsp_model import QSPModel, get_spec
from pkpd_agent.engines import llm_tasks as LT, llm_topology as TOP, llm_references as REF


def _build_references(args) -> str:
    """Assemble the reading material (Route A). Reads --refs verbatim, or parses --paper's
    reference list (optionally enriching with PubMed abstracts, cached next to the paper)."""
    if args.refs:
        with open(args.refs, encoding="utf-8") as fh:
            return fh.read().strip()
    if not args.paper:
        return ""
    with open(args.paper, encoding="utf-8") as fh:
        raw = fh.read()
    section = REF.extract_references_section(raw) or raw
    cites = REF.parse_references(section)
    print(f"parsed {len(cites)} citations from {os.path.basename(args.paper)}")
    if args.fetch_abstracts and cites:
        cache = os.path.splitext(args.paper)[0] + ".refs_full.txt"
        if os.path.isfile(cache) and not args.refetch:
            print(f"using cached abstracts: {cache}")
            with open(cache, encoding="utf-8") as fh:
                return fh.read().strip()
        print(f"fetching abstracts from PubMed for {len(cites)} citations "
              "(cached for reuse) ...", flush=True)
        REF.fetch_abstracts(cites)
        hits = sum(1 for c in cites if c.get("abstract"))
        print(f"got abstracts for {hits}/{len(cites)} citations")
        text = REF.format_references(cites)
        with open(cache, "w", encoding="utf-8") as fh:
            fh.write(text)
        return text
    return REF.format_references(cites)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="ra")
    ap.add_argument("--sbproj", required=True)
    ap.add_argument("--paper", default=None, help="paper text (its References are sliced out)")
    ap.add_argument("--refs", default=None, help="use this reference text verbatim instead")
    ap.add_argument("--fetch-abstracts", action="store_true", dest="fetch_abstracts",
                    help="Route A: enrich each citation with its PubMed abstract (cached)")
    ap.add_argument("--refetch", action="store_true", help="ignore the abstract cache")
    ap.add_argument("--web", action="store_true",
                    help="Route B: give the drafter web search so it reads references itself")
    ap.add_argument("--out", default="network.json", help="where to dump the answer key")
    ap.add_argument("--llm-model", default=None)
    args = ap.parse_args()

    cfg_llm = AgentConfig(mock=False)
    if args.llm_model:
        cfg_llm.model = args.llm_model
    if not cfg_llm.anthropic_key_present():
        print("ANTHROPIC_API_KEY not set.")
        return

    references = _build_references(args)
    if references:
        print(f"reading material: {len(references)} chars"
              + (" (+ web search)" if args.web else ""))
    elif args.web:
        print("no --paper/--refs: Route B will search the web from the node names alone")
    else:
        print("no references: LLM works from node names + biological priors alone")

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

    # answer key: influence edges, restricted to pairs BOTH among the given nodes (fair -
    # the LLM is only ever given the biological node list).
    all_truth = TOP.ground_truth_edges(network)
    nset = set(nodes)
    truth = {(s, d) for s, d in all_truth if s in nset and d in nset}
    print(f"answer key: {len(truth)} influence edges among the biological nodes "
          f"({len(all_truth)} total including drug/readout)")

    call = LT.default_web_call(cfg_llm) if args.web else LT.default_call(cfg_llm)
    print(f"== LLM drafting the topology ({'Route B: web' if args.web else 'Route A: cached'}) "
          "==", flush=True)
    draft = TOP.draft_topology(nodes, references, call)
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
