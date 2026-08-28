r"""Select the virtual-population parameter set the way the paper did - from the model's
FULL parameter list, by inferred mechanistic category - with no pre-curated list.

Enumerates every parameter in the loaded model (sb_params.m), then an LLM groups them by
the mechanistic category it infers from each name and selects the set to VARY across a
virtual population. This is the paper's category-based selection (it gave no numerical
threshold), done generally: nothing here names a disease. To show the choice is derived,
it reports the overlap with the project's hand-written vpop_drivers.

    python -m examples.run_llm_vpop_select --model ra ^
        --sbproj "..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj"

Needs ANTHROPIC_API_KEY and the MATLAB engine.
"""

from __future__ import annotations

import argparse
import os
from collections import Counter

from pkpd_agent.config import AgentConfig
from pkpd_agent.engines.simbiology import SimBiologyEngine
from pkpd_agent.engines import qsp_config, llm_tasks as LT


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="ra")
    ap.add_argument("--sbproj", required=True)
    ap.add_argument("--llm-model", default=None)
    args = ap.parse_args()

    cfg_llm = AgentConfig(mock=False)
    if args.llm_model:
        cfg_llm.model = args.llm_model
    if not cfg_llm.anthropic_key_present():
        print("ANTHROPIC_API_KEY not set.")
        return

    sb = SimBiologyEngine()
    try:
        print("== starting MATLAB engine =="); sb.start()
        print(f"== loading {os.path.basename(args.sbproj)} =="); sb.load_project(args.sbproj)
        print("== sb_params: enumerating every model parameter ==", flush=True)
        r = sb.list_parameters()
        ml = (r.get("matlab_log") or "").strip()
        if ml:
            print("   [MATLAB] " + ml.replace("\n", "\n   [MATLAB] "))
        params = r.get("parameters") or []
        print(f"model exposes {len(params)} parameters")
        if not params:
            print("no parameters returned - check the MATLAB log above.")
            return

        print("== LLM selecting the Vpop-varied set by inferred category ==", flush=True)
        sel = LT.propose_vpop_set(params, LT.default_call(cfg_llm))
        print(f"\nselected {sel['n_selected']} / {sel['n_candidates']} parameters to vary")
        if sel.get("rationale"):
            print("rationale:", sel["rationale"])
        cats = Counter(sel["categories"].values())
        if cats:
            print("\n== by inferred mechanistic category ==")
            for c, n in cats.most_common():
                print(f"  {n:4d}  {c}")

        # show the choice is DERIVED: overlap with the project's hand-written drivers
        try:
            cur = set(qsp_config.get(args.model).vpop_drivers)
            hit = sorted(n for n in sel["selected"] if n in cur)
            print(f"\noverlap with the config's {len(cur)} hand-written vpop_drivers: {hit}")
        except Exception:
            pass
    finally:
        sb.stop()


if __name__ == "__main__":
    main()
