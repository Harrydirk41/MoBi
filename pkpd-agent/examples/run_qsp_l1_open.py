r"""L1, UN-SCAFFOLDED: the agent gets a realistic, distractor-laden data inventory that does NOT
announce what is missing, and must itself work out - for each parameter it has to fit - whether
the data can pin it, what experiment would, and whether that experiment is actually present.
The un-scaffolded win is realising, unprompted, that the identifying data (per-cytokine dose-
responses) is ABSENT - not being told so.

    python -m examples.run_qsp_l1_open           (needs ANTHROPIC_API_KEY)

Grades three things: identifiability accuracy, whether it names the right missing experiment
(dose-response / perturbation / titration), and whether it flags that experiment as absent -
without ever being told the identifying data is missing.
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
    params = [{"name": f"FLSProlif_Maxby{c}",
               "has_literature_value": bool(prov.get(f"FLSProlif_Maxby{c}", {})
                                            .get("from_literature"))} for c in cyts]
    params += [{"name": "kd_FLS_Baseline", "has_literature_value": True},
               {"name": "kg_FLS_Baseline", "has_literature_value": False}]
    params += [{"name": f"HalfEffectConc_FLSProlif_by{c}", "has_literature_value": False}
               for c in cyts]
    params += [{"name": f"Slope_FLSProlif_by{c}", "has_literature_value": False} for c in cyts]

    # realistic inventory WITH distractors, and NO line announcing what is missing
    available = [
        "steady-state synovial FLS cell density in established RA (one disease state)",
        "steady-state synovial cytokine concentrations (TNFa, IL1b, IL17, TGFb, IL6)",
        "clinical trial ACR20/50/70 response rates for MTX, adalimumab, tocilizumab",
        "bulk RNA-seq of RA vs healthy synovium (differential expression)",     # distractor
        "serum autoantibody (RF, anti-CCP) titres across a patient cohort",     # distractor
        "a published DAS28-CRP formula relating joint counts and CRP to the score",  # distractor
    ]

    truth = {"kg_FLS_Baseline": True}
    for c in cyts:
        truth[f"HalfEffectConc_FLSProlif_by{c}"] = False
        truth[f"Slope_FLSProlif_by{c}"] = False

    print("== un-scaffolded: realistic inventory (with distractors), missing data NOT announced ==")
    plan = CP.propose_plan_open(params, available, LT.default_call(cfg))
    for param, v in plan["plan"].items():
        mark = "determinable    " if v["determinable"] else "NOT identifiable"
        av = {True: "have it", False: "MISSING", None: "?"}[v.get("needs_available")]
        print(f"  {mark}  {param:34} needs: {str(v.get('needs'))[:40]:40} [{av}]")

    g = CP.grade_open(plan, truth, needs_kw=["dose", "perturb", "titrat", "concentration-response",
                                             "stimul", "response curve", "gradient", "varying"])
    print(f"\n== grading (un-scaffolded) ==")
    print(f"  identifiability accuracy: {g['id_accuracy']}")
    print(f"  named the right missing experiment (dose-response/perturbation): "
          f"{g['named_missing_experiment']}")
    print(f"  flagged that experiment as ABSENT (unprompted): {g['flagged_data_absent']}")
    print(f"  OVER-claimed identifiable: {g['overclaimed'] or 'none'}")
    print("\n  -> the real test: without being told, does the agent realise the identifying "
          "data\n     (per-cytokine dose-responses) is not in its inventory? That is the "
          "no-oracle,\n     no-scaffold judgment.")


if __name__ == "__main__":
    main()
