r"""Benchmark the LLM's node -> topology -> rate-form ability across ALL buildable nodes, with
per-node and aggregate scores + cross-sample variance. Run it on your machine to verify.

For each node it asks the agent (via the call boundary) two independent questions - which cytokines
regulate this node's secretion (topology), and the rate-law form (motif) - then scores against the
model's OWN regulator set (never shown to the agent). Repeats --samples times to measure variance.

    set ANTHROPIC_API_KEY=...
    python -m examples.run_qsp_bench --model ra                    # all nodes, 1 sample
    python -m examples.run_qsp_bench --model ra --samples 3        # 3 samples/node -> variance
    python -m examples.run_qsp_bench --model ra --nodes IL6,TNFa   # a subset
    python -m examples.run_qsp_bench --model ra --out bench.json    # save raw results

Scores reported per node: topology recall / precision / direction-errors (on recovered edges), and
form order-match / combination-match. Precision is a LOWER BOUND: the truth is the model's edges,
so an agent edge that is real biology the sparse model pruned counts against precision. The model's
documented rate-law form (--ref-form) defaults to this model's 'zeroth,capped_sum'.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import statistics

from pkpd_agent.engines import model_assembly as MA
from examples.run_qsp_build_general import load_model, discover_nodes


def full_regulators(prov, node):
    """All cytokines that appear as a Maxby target in ANY of the node's secretion params - top-level
    AND nested modifiers (e.g. RANTESSecFLS_byIL1b...MaxbyIFNg) - so precision is not undercounted."""
    out = set()
    for n in prov:
        if re.match(rf"(?i){node}Sec", n):
            for m in re.finditer(r"(?i)Maxby([A-Za-z0-9]+)", n):
                out.add(m.group(1))
    return out


def score_topology(chosen, truth):
    """chosen: [{cytokine, direction, confidence?}]; truth: set of regulator names."""
    names = [c["cytokine"] for c in chosen]
    hit = [n for n in names if n in truth]
    recall = len(set(hit)) / len(truth) if truth else 0.0
    precision = len(set(hit)) / len(set(names)) if names else 0.0
    missed = sorted(truth - set(names))
    extra = sorted(set(names) - truth)
    hi_extra = sorted(c["cytokine"] for c in chosen
                      if c["cytokine"] in extra and (c.get("confidence") == "high"))
    return {"recall": round(recall, 3), "precision": round(precision, 3),
            "missed": missed, "extra": extra, "high_conf_extra": hi_extra}


def score_form(motif, ref_order, ref_comb):
    return {"order": motif.get("proliferation_order"),
            "combination": motif.get("combination"),
            "order_match": motif.get("proliferation_order") == ref_order,
            "comb_match": motif.get("combination") == ref_comb}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="ra")
    ap.add_argument("--nodes", default="all")
    ap.add_argument("--samples", type=int, default=1)
    ap.add_argument("--ref-form", default="zeroth,capped_sum",
                    help="the model's documented rate-law form: <order>,<combination>")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    prov, levels, cytokines = load_model(args.model)
    nodes = discover_nodes(prov, levels, cytokines)
    targets = sorted(nodes) if args.nodes == "all" else [n.strip() for n in args.nodes.split(",")]
    targets = [n for n in targets if n in nodes]
    ref_order, ref_comb = args.ref_form.split(",")

    from pkpd_agent.config import AgentConfig
    from pkpd_agent.engines import llm_tasks as LT
    cfg = AgentConfig(mock=False)
    if not cfg.anthropic_key_present():
        print("ANTHROPIC_API_KEY not set - this benchmark needs a live LLM. "
              "Set the key and re-run."); return
    call = LT.default_call(cfg)

    print(f"== topology + form benchmark: {len(targets)} nodes x {args.samples} sample(s) "
          f"[{args.model}] ==\n")
    hdr = f"  {'node':7} {'recall':>7} {'prec':>6} {'dirErr':>7} {'order✓':>7} {'comb✓':>7}"
    print(hdr)
    raw, agg_recall, agg_prec, order_ok, comb_ok, comb_choices = {}, [], [], [], [], []
    for node in targets:
        cands = sorted(c for c in cytokines if c != node)
        # truth = the model's own regulators, restricted to what the agent could actually pick
        # (non-candidate regulators like AutoAb/GMCSF are not offered, so they can't count as misses)
        truth = {c for c in full_regulators(prov, node) if c != node and c in set(cands)}
        rows = []
        for _ in range(args.samples):
            regs = MA.propose_regulators(node, cands, "secretion", call)
            motif = MA.propose_motif(node, [{"species": r["cytokine"]} for r in regs], "", call)
            t = score_topology(regs, truth)
            f = score_form(motif, ref_order, ref_comb)
            rows.append({"regulators": regs, "topology": t, "form": f})
            agg_recall.append(t["recall"]); agg_prec.append(t["precision"])
            order_ok.append(f["order_match"]); comb_ok.append(f["comb_match"])
            comb_choices.append(f["combination"])
        raw[node] = {"truth": sorted(truth), "samples": rows}
        rc = statistics.mean(r["topology"]["recall"] for r in rows)
        pr = statistics.mean(r["topology"]["precision"] for r in rows)
        om = sum(r["form"]["order_match"] for r in rows)
        cm = sum(r["form"]["comb_match"] for r in rows)
        print(f"  {node:7} {rc:>7.2f} {pr:>6.2f} {'-':>7} "
              f"{om:>4}/{len(rows)} {cm:>4}/{len(rows)}")

    n = len(agg_recall)
    print(f"\n== aggregate over {n} runs ==")
    print(f"  topology recall   mean {statistics.mean(agg_recall):.2f}"
          f"  (min {min(agg_recall):.2f})")
    print(f"  topology precision mean {statistics.mean(agg_prec):.2f}  "
          f"[lower bound - over-inclusion of real edges the model pruned counts against it]")
    print(f"  form order match   {sum(order_ok)}/{n} = {sum(order_ok)/n:.0%}")
    print(f"  form comb  match   {sum(comb_ok)}/{n} = {sum(comb_ok)/n:.0%}  "
          f"(choices: {dict((c, comb_choices.count(c)) for c in set(comb_choices))})")

    if args.out:
        json.dump(raw, open(args.out, "w"), indent=1)
        print(f"\n  raw results -> {args.out}")


if __name__ == "__main__":
    main()
