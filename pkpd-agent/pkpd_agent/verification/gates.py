"""Verification gates.

Each gate is a function ``(tool_name, result_data, session) -> list[Finding]``.
Gates are cheap, deterministic, and specific. They encode the sanity checks a
careful pharmacometrician applies before believing a result.
"""

from __future__ import annotations

from typing import Any, Callable

from ..state import Finding

Gate = Callable[[str, dict[str, Any], "ModelingSession"], list[Finding]]  # noqa: F821


# --------------------------------------------------------------------------- #
# Estimation gates (pharmpy)
# --------------------------------------------------------------------------- #

def fit_convergence_gate(tool: str, data: dict[str, Any], session) -> list[Finding]:
    if tool != "pharmpy_fit":
        return []
    out: list[Finding] = []
    if data.get("minimization_successful") is False:
        out.append(Finding(
            "block", "fit_convergence",
            "Minimization did NOT succeed - estimates are unreliable. "
            "Reconsider initial estimates, model structure, or estimation method "
            "before trusting this fit.",
        ))
    cond = data.get("condition_number")
    if isinstance(cond, (int, float)) and cond > 1000:
        out.append(Finding(
            "warn", "fit_conditioning",
            f"High condition number ({cond:.0f}) suggests near-collinear "
            "parameters / weak identifiability.",
        ))
    return out


def rse_plausibility_gate(tool: str, data: dict[str, Any], session) -> list[Finding]:
    if tool != "pharmpy_fit":
        return []
    out: list[Finding] = []
    for name, rse in (data.get("relative_standard_errors") or {}).items():
        if isinstance(rse, (int, float)) and rse > 0.50:
            out.append(Finding(
                "warn", "parameter_precision",
                f"Parameter {name} has RSE {rse:.0%} (>50%) - poorly estimated; "
                "may be over-parameterized.",
            ))
    return out


# --------------------------------------------------------------------------- #
# Mechanistic gates (OSP)
# --------------------------------------------------------------------------- #

def physical_sanity_gate(tool: str, data: dict[str, Any], session) -> list[Finding]:
    if tool != "osp_simulate":
        return []
    out: list[Finding] = []
    if data.get("all_values_finite") is False:
        out.append(Finding(
            "block", "numerical_integrity",
            "Simulation produced non-finite values (NaN/Inf) - the ODE system "
            "is misspecified or the solver diverged.",
        ))
    min_c = data.get("min_concentration")
    if isinstance(min_c, (int, float)) and min_c < 0:
        out.append(Finding(
            "block", "physical_sanity",
            f"Negative concentration ({min_c}) is unphysical - check the "
            "reaction/transport definitions.",
        ))
    resid = data.get("mass_balance_residual")
    if isinstance(resid, (int, float)) and abs(resid) > 1e-3:
        out.append(Finding(
            "warn", "mass_balance",
            f"Mass-balance residual {resid:.2e} exceeds tolerance - amount may "
            "not be conserved across compartments.",
        ))
    return out


# --------------------------------------------------------------------------- #
# VPC gate
# --------------------------------------------------------------------------- #

def vpc_coverage_gate(tool: str, data: dict[str, Any], session) -> list[Finding]:
    if tool != "pharmpy_vpc":
        return []
    inside = data.get("pct_observations_within_90_pi")
    if isinstance(inside, (int, float)) and inside < 80:
        return [Finding(
            "warn", "vpc_coverage",
            f"Only {inside}% of observations fall within the 90% PI - model "
            "may misfit the data.",
        )]
    return []


DEFAULT_GATES: list[Gate] = [
    fit_convergence_gate,
    rse_plausibility_gate,
    physical_sanity_gate,
    vpc_coverage_gate,
]


def run_gates(tool: str, data: dict[str, Any], session,
              gates: list[Gate] | None = None) -> list[Finding]:
    gates = gates if gates is not None else DEFAULT_GATES
    findings: list[Finding] = []
    for gate in gates:
        findings.extend(gate(tool, data, session))
    return findings
