r"""L1, agent-driven: give the AGENT the real-data inventory (from MOESM1/MOESM2) and let IT
plan the calibration and judge identifiability - which missing parameters the data can pin and
which it cannot - then grade the agent's judgment against the known truth. This is the no-oracle
test: not us proving the identifiability, but whether the agent works it out itself.

    python -m examples.run_qsp_l1_agent           (needs ANTHROPIC_API_KEY)

Truth (for FLS proliferation, from the dynamics): the baseline rate kg IS identifiable from the
steady-state target (it sets the overall scale); the half-effect K's and Hill slopes are NOT
identifiable from a single steady-state target (they trade off, kg re-absorbs them) - they need
per-cytokine dose-responses. The dangerous failure is the agent OVER-claiming it can pin them.
"""

from __future__ import annotations

import json
import os

from pkpd_agent.config import AgentConfig
from pkpd_agent.engines import llm_tasks as LT, llm_calib_plan as CP

_DATA = os.path.join(os.path.dirname(__file__), "..", "projects", "vantage_ra", "data")


def main() -> None:
    cfg = AgentConfig(mock=False)
    if not cfg.anthropic_key_present():
        print("ANTHROPIC_API_KEY not set."); return

    prov = {p["name"]: p for p in json.load(open(os.path.join(_DATA, "param_provenance.json")))}
    cyts = ["TNFa", "IL1b", "IL17", "TGFb", "IL6"]
    known = [f"FLSProlif_Maxby{c}" for c in cyts if prov.get(f"FLSProlif_Maxby{c}", {})
             .get("from_literature")] + ["kd_FLS_Baseline"]
    missing = ["kg_FLS_Baseline"] + [f"HalfEffectConc_FLSProlif_by{c}" for c in cyts] \
        + [f"Slope_FLSProlif_by{c}" for c in cyts]
    available = [
        "FLS steady-state cell density (one disease steady state), from MOESM1",
        "steady-state cytokine levels (TNFa, IL1b, IL17, TGFb, IL6), from MOESM1",
        "clinical trial ACR response rates for MTX/ADA/TCZ (population level), from MOESM1",
        "NO per-cytokine dose-response experiments are available",
    ]

    # ground-truth identifiability: only the baseline rate is pinned by the steady-state target
    truth = {"kg_FLS_Baseline": True}
    for c in cyts:
        truth[f"HalfEffectConc_FLSProlif_by{c}"] = False
        truth[f"Slope_FLSProlif_by{c}"] = False

    print("== giving the agent the real-data inventory; it must judge what is identifiable ==")
    print(f"  known (literature): {len(known)}   missing (to infer): {len(missing)}")
    plan = CP.propose_plan(known, missing, available, LT.default_call(cfg))

    print("\n== the agent's identifiability verdicts ==")
    for param, v in plan["plan"].items():
        mark = "determinable    " if v["determinable"] else "NOT identifiable"
        print(f"  {mark}  {param:36} {v.get('reason') or ''}")

    g = CP.grade_plan(plan, truth)
    print(f"\n== grading the agent's judgment against the truth ==")
    print(f"  accuracy {g['accuracy']}  ({g['correct']}/{g['n']})")
    if g["wrong"]:
        print("  got wrong:")
        for param, why in g["wrong"]:
            print(f"    {param}: {why}")
    print(f"\n  OVER-claimed (said identifiable when it is NOT - the dangerous error): "
          f"{g['overclaimed'] or 'none'}")
    print("\n  -> this tests the no-oracle judgment: without an answer key, does the agent KNOW "
          "which\n     parameters the real data can and cannot pin? Over-claiming is what would "
          "silently\n     corrupt a from-scratch model.")


if __name__ == "__main__":
    main()
