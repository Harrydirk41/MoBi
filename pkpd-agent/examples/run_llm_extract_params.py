r"""Data-acquisition benchmark (A): for each model parameter that cites a reference, have a
web-enabled LLM read the paper and extract the value; grade against MOESM2's "value from
reference". Reports the two rates the modelling bottleneck hinges on: RETRIEVAL (found the
paper) and EXTRACTION accuracy (got the number, within tolerance) - and flags how often the
value is figure-only (needs vision, not text).

    python -m examples.run_llm_extract_params ^
        --provenance projects\vantage_ra\data\param_provenance.json ^
        --limit 20 --web

--web uses the server-side web-search tool (default_web_call). Needs ANTHROPIC_API_KEY.
"""

from __future__ import annotations

import argparse
import json

from pkpd_agent.config import AgentConfig
from pkpd_agent.engines import llm_tasks as LT, llm_extract as EX


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--provenance", required=True)
    ap.add_argument("--limit", type=int, default=0, help="0 = all parameters with a reference")
    ap.add_argument("--tol", type=float, default=0.25)
    ap.add_argument("--web", action="store_true", help="give the LLM web search + fetch")
    ap.add_argument("--figures", default="", help="dir of figure images named <param>.png/.jpg; "
                    "when present, use VISION (read the value off the figure image)")
    ap.add_argument("--llm-model", default=None)
    args = ap.parse_args()

    cfg = AgentConfig(mock=False)
    if args.llm_model:
        cfg.model = args.llm_model
    if not cfg.anthropic_key_present():
        print("ANTHROPIC_API_KEY not set."); return

    with open(args.provenance, encoding="utf-8") as fh:
        prov = json.load(fh)
    items = [p for p in prov if p.get("reference") and p.get("value_from_reference") is not None]
    if args.limit:
        items = items[: args.limit]
    print(f"grading extraction on {len(items)} parameters that cite a reference\n")

    import glob
    web_call = LT.default_web_call(cfg) if args.web else LT.default_call(cfg)
    vis_call = LT.default_vision_call(cfg) if args.figures else None
    found = extracted = hit = figonly = vis_used = 0
    for i, p in enumerate(items, 1):
        imgs = []
        if args.figures:
            imgs = sorted(glob.glob(os.path.join(args.figures, p["name"] + ".*")))
        try:
            if imgs:
                r = EX.extract_value_vision(p, imgs, vis_call); vis_used += 1
            else:
                r = EX.extract_value(p, web_call)
        except Exception as e:
            print(f"[{i}/{len(items)}] {p['name'][:26]:26} ERROR {e}"); continue
        g = EX.grade(r["value"], p["value_from_reference"], tol=args.tol)
        found += r["found_paper"]; extracted += g["extracted"]; hit += g["hit"]
        figonly += 1 if r.get("in_figure_only") else 0
        mark = "HIT " if g["hit"] else ("MISS" if g["extracted"] else "----")
        via = "vision" if imgs else ("fig-only" if r.get("in_figure_only") else "text")
        print(f"[{i}/{len(items)}] {mark} {p['name'][:24]:24} "
              f"got={r['value']} truth={p['value_from_reference']} ({via}) {p['reference']}")
        if r.get("note"):
            print(f"          reason: {r['note']}")

    n = len(items) or 1
    print(f"\n== data-acquisition scorecard ==")
    print(f"  retrieval  (found the paper):        {found}/{n} = {found/n:.0%}")
    print(f"  extraction (returned any value):     {extracted}/{n} = {extracted/n:.0%}")
    print(f"  accuracy   (within {args.tol:.0%} of truth):    {hit}/{n} = {hit/n:.0%}")
    print(f"  value was figure-only (needs vision): {figonly}/{n} = {figonly/n:.0%}")
    if vis_used:
        print(f"  (of these, {vis_used} used a supplied figure image via VISION)")
    print("  -> the gap between retrieval and accuracy is the data-acquisition bottleneck: "
          "reaching\n     a paper is not reading its figures. Supply figures (--figures) to "
          "test the vision ceiling.")


if __name__ == "__main__":
    main()
