r"""Structure-discovery benchmark (B): can the agent REDISCOVER a load-bearing edge that the
literature cannot supply? Ablate one regulatory edge from the calibrated model, read the
resulting disease-steady-state SYMPTOM (which species moved, and which way), then ask the LLM
- given only the symptom and a candidate edge set - which missing edge explains it. Grade by
where the LLM ranks the edge actually removed.

    python -m examples.run_qsp_discover --model ra ^
        --sbproj "..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj" ^
        --edge "IL12->IL6" --distractors 8 --llm

--edge picks the edge to remove (default: the highest-impact regulatory edge found by a quick
knockout scan). --llm ranks candidates (needs ANTHROPIC_API_KEY); without it the symptom + the
mechanical answer are still printed. Needs the MATLAB engine.
"""

from __future__ import annotations

import argparse
import os

from pkpd_agent.engines.simbiology import SimBiologyEngine
from pkpd_agent.engines.qsp_model import QSPModel, get_spec
from pkpd_agent.engines import llm_topology as TOP, llm_discover as D


def _parse_edge(s: str):
    for sep in ("->", "→", ","):
        if sep in s:
            a, b = s.split(sep, 1)
            return (a.strip(), b.strip())
    return None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="ra")
    ap.add_argument("--sbproj", required=True)
    ap.add_argument("--edge", default="", help="edge to remove, e.g. 'IL12->IL6'")
    ap.add_argument("--distractors", type=int, default=8)
    ap.add_argument("--readout-day", type=float, default=199.0)
    ap.add_argument("--top-symptom", type=int, default=6, help="how many moved species to show")
    ap.add_argument("--llm", action="store_true")
    ap.add_argument("--llm-model", default=None)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    sb = SimBiologyEngine()
    try:
        print("== starting MATLAB =="); sb.start()
        print(f"== loading {os.path.basename(args.sbproj)} =="); sb.load_project(args.sbproj)
        net = sb.network_json(_tmp := os.path.join(os.getcwd(), "network.json"))
        model = QSPModel(net, get_spec(args.model))
        nset = set(model.nodes)

        # clean per-edge knobs: '<Dest>...Maxby<Src>' fold-change constants. Setting one to 1.0
        # removes exactly that edge (not the whole combined-effect rule), so the symptom is
        # localized instead of a runaway - the only cleanly ablatable edges in this model.
        pnames = [p["name"] for p in sb.list_parameters().get("parameters", [])]
        medges = D.maxby_edges(pnames, nset)
        if not medges:
            print("no Maxby per-edge knobs found in the model."); return

        true_edge = _parse_edge(args.edge) if args.edge else None
        if true_edge and true_edge not in medges:
            print(f"edge {true_edge} has no clean per-edge knob. Cleanly ablatable edges:")
            for e in list(medges)[:20]:
                print(f"    {e[0]} -> {e[1]}   ({medges[e]})")
            return
        if not true_edge:
            true_edge = next(iter(medges))
        knob = medges[true_edge]
        p0 = next((float(p["value"]) for p in sb.list_parameters().get("parameters", [])
                   if p["name"] == knob), 1.0)
        print(f"removing edge {true_edge[0]} -> {true_edge[1]} "
              f"(set {knob} {p0:g} -> 1.0 = no effect)")

        # symptom: species that move between intact and ablated disease steady state
        base = {k: v[-1] for k, v in sb.simulate(stop_time=args.readout_day + 1.0)
                .get("columns", {}).items() if v}
        sb.set_parameter(knob, 1.0)
        abl = {k: v[-1] for k, v in sb.simulate(stop_time=args.readout_day + 1.0)
               .get("columns", {}).items() if v}
        sb.set_parameter(knob, p0)                      # restore
        moved = []
        for sp in base:
            b, a = base.get(sp), abl.get(sp)
            if b and a is not None and sp in nset and b != 0:
                rel = (a - b) / abs(b)
                if abs(rel) > 0.02:
                    moved.append({"species": sp, "direction": "too low" if rel < 0 else "too high",
                                  "rel_change": f"{rel:+.0%}"})
        moved.sort(key=lambda m: abs(float(m["rel_change"].rstrip("%")) / 100), reverse=True)
        symptom = moved[: args.top_symptom]
        print(f"\n== symptom (ablated vs intact disease steady state) ==")
        for s in symptom:
            print(f"    {s['species']:12} {s['direction']:9} {s['rel_change']}")
        if not symptom:
            print("  (no species moved > 2% - this edge is not observable at steady state)")
            return

        candidates = D.candidate_set(true_edge, list(medges), args.distractors, args.seed)
        print(f"\n== {len(candidates)} candidate missing edges (true one hidden among them) ==")
        for s, d in candidates:
            print(f"    {s} -> {d}")

        if args.llm:
            from pkpd_agent.config import AgentConfig
            from pkpd_agent.engines import llm_tasks as LT
            cfg = AgentConfig(mock=False)
            if args.llm_model:
                cfg.model = args.llm_model
            if not cfg.anthropic_key_present():
                print("ANTHROPIC_API_KEY not set for --llm."); return
            print("\n== LLM ranking candidates from the symptom ==", flush=True)
            ranking = D.rank_candidates(symptom, candidates, LT.default_call(cfg))
            pos = D.rank_of_true(ranking, true_edge)
            print("  LLM ranking (most-likely first):")
            for i, (s, d) in enumerate(ranking, 1):
                mark = "  <-- TRUE" if (s, d) == true_edge else ""
                print(f"    {i}. {s} -> {d}{mark}")
            print(f"\n  the removed edge was ranked #{pos} of {len(candidates)} "
                  + ("(top pick - discovered!)" if pos == 1 else
                     f"(random baseline would be ~{(len(candidates)+1)/2:.1f})"))
    finally:
        sb.stop()


if __name__ == "__main__":
    main()
