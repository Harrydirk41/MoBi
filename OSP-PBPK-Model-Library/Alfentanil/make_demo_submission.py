"""Make a SYNTHETIC demo submission so the grader can be exercised end-to-end.

This does NOT run a PBPK model. It fabricates a plausible-looking submission by
perturbing the observed data (a route-dependent fold factor) so the pipeline
produces a realistic scorecard, and includes a parameter list with one
deliberately implausible value to show the plausibility checker firing.

    python make_demo_submission.py json_input/Alfentanil-Model.input.json \\
        --out demo_submission.json
"""

from __future__ import annotations

import argparse
import json


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("input")
    ap.add_argument("--out", default="demo_submission.json")
    args = ap.parse_args()

    with open(args.input, encoding="utf-8") as fh:
        task = json.load(fh)
    obs = task.get("given_data", {}).get("clinical_observed_data", [])

    # perturb: IV slightly high, PO slightly low -> non-trivial but decent GMFE
    factor = {"IV": 1.25, "PO": 0.8}
    profiles = []
    for o in obs:
        f = factor.get((o.get("route") or "").upper(), 1.15)
        pred = [None if c is None else round(c * f, 8) for c in o.get("conc_mg_L", [])]
        profiles.append({"dataset": o["dataset"], "time_h": o.get("time_h", []),
                         "pred_conc_mg_L": pred})

    submission = {
        "_synthetic": "Fabricated by make_demo_submission.py for pipeline testing; "
                      "not a real simulation.",
        "submission": {
            "structural_model": {
                "distribution_model": "Rodgers & Rowland partition coefficients",
                "absorption_model": "passive transcellular permeability",
                "elimination_pathways": [
                    {"pathway": "CYP3A4 hepatic metabolism", "type": "first-order intrinsic clearance"},
                    {"pathway": "glomerular filtration", "type": "GFR fraction"},
                ],
                "engine": "OSP PK-Sim (illustrative)",
                "notes": "Whole-body PBPK, typical European adult.",
            },
            "parameters": [
                {"parameter": "Lipophilicity", "value": 1.9, "unit": "Log Units",
                 "fixed_or_estimated": "estimated", "rationale": "tuned near literature logD 2.1"},
                {"parameter": "Fraction unbound (plasma)", "value": 0.1, "unit": "",
                 "fixed_or_estimated": "estimated", "rationale": "within literature 8.6-12%"},
                {"parameter": "CYP3A4 intrinsic clearance", "value": 0.5, "unit": "l/min",
                 "fixed_or_estimated": "estimated", "rationale": "fit to IV disposition"},
                {"parameter": "Specific intestinal permeability", "value": 6e-4, "unit": "cm/min",
                 "fixed_or_estimated": "estimated", "rationale": "fit to oral absorption"},
                # deliberately implausible to demonstrate the checker:
                {"parameter": "Fraction unbound (alt)", "value": 1.4, "unit": "",
                 "fixed_or_estimated": "estimated", "rationale": "INTENTIONAL bad value for demo"},
            ],
            "predicted_profiles": profiles,
            "self_assessment": {
                "data_fit_gmfe_overall": None,
                "pct_within_2fold": None,
                "parameter_plausibility_notes": "one alt fu value is out of range (demo)",
                "output_plausibility_notes": "profiles mirror observed shape by construction",
            },
        },
    }
    with open(args.out, "w", encoding="utf-8") as fh:
        json.dump(submission, fh, indent=2, ensure_ascii=False)
    print(f"wrote {args.out}  ({len(profiles)} predicted profiles)")


if __name__ == "__main__":
    main()
