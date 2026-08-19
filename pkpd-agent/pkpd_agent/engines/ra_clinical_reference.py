"""Real-world RA clinical-trial outcomes, transcribed from the paper's ESM1
(41540_2024_454_MOESM1_ESM.xlsx, sheet 'Clinical_trials').

These are the ACTUAL trial numbers the Vantage RA QSP model was calibrated to -
the ground truth for asking "does the model (and an agent driving it) match
reality?", as opposed to "does the agent reproduce the model's own output".

Two conventions are provided per arm:
  * ``raw``   - the reported ACR responder percentages (drug arm).
  * ``pcorr`` - placebo-corrected = drug arm - placebo arm, the drug-attributable
    effect the paper follows (Wang et al. method; see paper Methods Stage 4/5).
The QSP Vpop carries its own disease-progression baseline, so its drug-arm output
tends to line up with ``raw`` rather than ``pcorr`` - but both are kept so the
comparison is explicit rather than assumed.

Endpoints are percentages of patients: ACR20/ACR50/ACR70 and DAS28-CRP<2.6
remission, at the trial's reported week.
"""

from __future__ import annotations

from typing import Any, Optional

# drug arm and placebo arm, by trial, at each reported week (ACR20/50/70/rem %)
RA_TRIALS: dict[str, dict[str, Any]] = {
    "MTX": {
        "trial": "Strand et al. 1999 (MTX monotherapy)",
        "population": "MTX-naive RA",
        "weeks": {
            12: {"drug": {"ACR20": 46.0, "ACR50": 23.0, "ACR70": 9.0},
                 "placebo": {"ACR20": 26.0, "ACR50": 8.0, "ACR70": 4.0}},
        },
    },
    "ADA": {
        "trial": "OPTIMA / Kavanaugh et al. (ADA 40mg Q2W + MTX)",
        "population": "MTX-background, biologic-naive",
        "weeks": {
            24: {"drug": {"ACR20": 70.0, "ACR50": 52.0, "ACR70": 35.0, "rem": 34.0},
                 "placebo": {"ACR20": 57.0, "ACR50": 34.0, "ACR70": 17.0, "rem": 17.0}},
        },
    },
    "TCZ": {
        "trial": "ROSE / Yazici et al. 2012 (TCZ 8mg/kg Q4W + DMARDs)",
        "population": "DMARD-IR (~38% prior anti-TNF)",
        "weeks": {
            12: {"drug": {"ACR20": 50.0, "ACR50": 25.0, "ACR70": 12.0, "rem": 22.0},
                 "placebo": {"ACR20": 29.0, "ACR50": 12.0, "ACR70": 5.0, "rem": 3.0}},
            24: {"drug": {"ACR20": 45.0, "ACR50": 29.0, "ACR70": 13.9, "rem": 38.0},
                 "placebo": {"ACR20": 25.0, "ACR50": 10.0, "ACR70": 1.9, "rem": 2.0}},
        },
    },
}


def trial_target(drug: str, week: int, correction: str = "raw") -> Optional[dict]:
    """Return {ACR20, ACR50, ACR70[, rem]} for a drug/week. ``correction`` is
    'raw' (drug arm) or 'pcorr' (drug minus placebo, floored at 0)."""
    arm = (RA_TRIALS.get(drug) or {}).get("weeks", {}).get(week)
    if not arm:
        return None
    drug_arm = arm["drug"]
    if correction == "raw":
        return {k: v for k, v in drug_arm.items() if k in ("ACR20", "ACR50", "ACR70")}
    pl = arm.get("placebo", {})
    return {k: round(max(0.0, drug_arm[k] - pl.get(k, 0.0)), 1)
            for k in ("ACR20", "ACR50", "ACR70") if k in drug_arm}


def describe(drug: str) -> str:
    t = RA_TRIALS.get(drug) or {}
    return f"{t.get('trial', drug)} [{t.get('population', '')}]"
