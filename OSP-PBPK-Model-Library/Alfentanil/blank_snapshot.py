"""Create a BLANKED Alfentanil snapshot: the real masked-benchmark start point.

The reference snapshot is a *finished* model - it already contains the fitted
parameter values and the chosen distribution method, i.e. the answers. For a
genuine model-building benchmark the agent must start from a model that is
structurally set up but NOT yet fitted, and recover the parameters (and, if it
wants, the structural choices) from the observed data.

This script reads the reference snapshot, finds every compound parameter that
was ESTIMATED (ValueOrigin = ParameterIdentification), and resets it to a
documented naive prior - a value a modeler would plausibly start from before
fitting. It writes:

  * benchmark/Alfentanil-Model.blanked.json   - the agent's starting model
  * answer_key/Alfentanil-Model.answer_edits.json - the edit spec that reverses
        the blanking (the fitted values), used only to PROVE the edit path can
        reconstruct the reference (reversibility test) - never shown to the agent.

The blanked model still compiles and runs in PK-Sim; it just fits the data
worse, leaving headroom for the agent to close.

    python blank_snapshot.py            # uses json/Alfentanil-Model.json
"""

from __future__ import annotations

import json
import os

REF = "json/Alfentanil-Model.json"
BLANKED = "benchmark/Alfentanil-Model.blanked.json"
ANSWER_EDITS = "answer_key/Alfentanil-Model.answer_edits.json"

# Naive priors for the estimated parameters (documented starting guesses, NOT
# the fitted values). Names must match the compound parameter names exactly.
PRIORS = {
    "Lipophilicity": 2.1,                    # literature logD (fitted ~1.85)
    "Fraction unbound (plasma, reference value)": 0.1,   # literature midpoint
    "Specific intestinal permeability (transcellular)": 1.0e-5,  # low prior
    "Permeability": 1.0e-3,                  # generic organ permeability prior
    "Intrinsic clearance": 0.1,              # under-clearing prior (fitted ~0.53)
}


def _fitted_parameters(comp: dict) -> dict[str, float]:
    """Every compound parameter whose value came from ParameterIdentification."""
    found: dict[str, float] = {}

    def walk(o):
        if isinstance(o, dict):
            if (isinstance(o.get("Value"), (int, float))
                    and (o.get("ValueOrigin") or {}).get("Source")
                    == "ParameterIdentification"):
                found[o["Name"]] = o["Value"]
            for v in o.values():
                walk(v)
        elif isinstance(o, list):
            for v in o:
                walk(v)

    walk(comp)
    return found


def _set_values(comp: dict, values: dict[str, float]) -> list[str]:
    changed = []

    def walk(o):
        if isinstance(o, dict):
            nm = o.get("Name")
            if isinstance(nm, str) and nm in values and "Value" in o:
                o["Value"] = values[nm]
                o.pop("ValueOrigin", None)
                changed.append(nm)
            for v in o.values():
                walk(v)
        elif isinstance(o, list):
            for v in o:
                walk(v)

    walk(comp)
    return changed


def main() -> None:
    with open(REF, encoding="utf-8") as fh:
        snap = json.load(fh)
    comp = snap["Compounds"][0]

    fitted = _fitted_parameters(comp)
    print(f"fitted (estimated) parameters in reference: {len(fitted)}")
    for k, v in fitted.items():
        prior = PRIORS.get(k, "(kept - no prior defined)")
        print(f"   {k:48} fitted={v:<14} -> prior={prior}")

    # blank: set fitted params to naive priors
    priors_to_apply = {k: v for k, v in PRIORS.items() if k in fitted}
    changed = _set_values(comp, priors_to_apply)

    os.makedirs(os.path.dirname(BLANKED), exist_ok=True)
    with open(BLANKED, "w", encoding="utf-8") as fh:
        json.dump(snap, fh, ensure_ascii=False, indent=1)
    print(f"\nwrote {BLANKED}  (blanked {len(changed)} parameters)")

    # answer-edit spec that reverses the blanking (reversibility test only)
    answer_edits = {"parameters": fitted}
    os.makedirs(os.path.dirname(ANSWER_EDITS), exist_ok=True)
    with open(ANSWER_EDITS, "w", encoding="utf-8") as fh:
        json.dump(answer_edits, fh, ensure_ascii=False, indent=2)
    print(f"wrote {ANSWER_EDITS}  (fitted values - GRADER ONLY, do not show agent)")

    print("\nReversibility test (proves the edit path reconstructs the model):")
    print("  # blanked alone -> worse GMFE:")
    print("  python -m examples.osp_run --snapshot "
          "..\\OSP-PBPK-Model-Library\\Alfentanil\\benchmark\\Alfentanil-Model.blanked.json \\")
    print("      --input ..\\OSP-PBPK-Model-Library\\Alfentanil\\json_input\\Alfentanil-Model.input.json")
    print("  # blanked + answer edits -> should return to ~1.45 GMFE:")
    print("      ... --edits ..\\OSP-PBPK-Model-Library\\Alfentanil\\answer_key\\Alfentanil-Model.answer_edits.json")


if __name__ == "__main__":
    main()
