"""Authoritative catalog of the OSP PBPK action space.

The agent should not need to *know* OSP a priori. The parameter names, current
structure, and expressed molecules come authoritatively from the snapshot; this
module supplies the descriptions, plausible ranges, and the menu of legal
calculation methods and addable process types, so an LLM with only general
pharmacology knowledge can decide what to do.

Parameter descriptions/ranges are general PBPK knowledge (not drug-specific).
Process structures are the OSP snapshot forms verified against PK-Sim 12.
"""

from __future__ import annotations

from typing import Any

# --------------------------------------------------------------------------- #
# compound parameters: description, plausible range, typical role
# --------------------------------------------------------------------------- #
# role: "measured"  -> normally fixed at the literature value
#       "estimate"  -> commonly identified (IVIVE unreliable)
PARAM_CATALOG: dict[str, dict[str, Any]] = {
    "Molecular weight": {
        "description": "molecular weight", "role": "measured"},
    "Solubility at reference pH": {
        "description": "aqueous solubility at the reference pH", "role": "measured"},
    "Reference pH": {
        "description": "pH at which the solubility is defined", "role": "measured"},
    "Lipophilicity": {
        "description": "effective lipophilicity / membrane affinity (logP-like) "
        "driving tissue partitioning; the model value often differs from measured "
        "logD, so it is commonly refined",
        "range": [-2.0, 7.0], "role": "estimate"},
    "Fraction unbound (plasma, reference value)": {
        "description": "unbound fraction in plasma; scales distribution and "
        "clearance. Measured, but small errors have large effect - sometimes refined",
        "range": [0.001, 1.0], "role": "measured"},
    "Specific intestinal permeability (transcellular)": {
        "description": "transcellular intestinal permeability governing the rate "
        "and extent of oral absorption; in-vitro->in-vivo transfer is weak, so "
        "usually estimated",
        "range": [1e-8, 1e-2], "unit": "cm/min", "role": "estimate"},
    "Permeability": {
        "description": "specific organ (cellular) permeability governing the "
        "kinetics of tissue distribution; usually estimated",
        "range": [1e-6, 1.0], "unit": "cm/min", "role": "estimate"},
    "Intrinsic clearance": {
        "description": "first-order metabolic intrinsic clearance by the enzyme "
        "(per the process it belongs to); IVIVE unreliable, so usually estimated",
        "range": [1e-3, 5.0], "unit": "l/min", "role": "estimate"},
    "GFR fraction": {
        "description": "fraction of glomerular filtration contributing to renal "
        "clearance of unbound drug (0 disables renal clearance)",
        "range": [0.0, 1.0], "role": "measured"},
    "Vmax": {"description": "maximum rate of a Michaelis-Menten process",
             "range": [1e-3, 1e3], "role": "estimate"},
    "Km": {"description": "Michaelis constant (half-saturation concentration)",
           "range": [1e-4, 1e3], "role": "estimate"},
}


def describe_parameter(name: str) -> dict[str, Any]:
    return dict(PARAM_CATALOG.get(name, {}))


def param_role(name: str) -> str:
    return PARAM_CATALOG.get(name, {}).get("role", "unknown")


# --------------------------------------------------------------------------- #
# calculation methods (with descriptions)
# --------------------------------------------------------------------------- #
PARTITION_METHOD_INFO = {
    "PK-Sim Standard": "Willmann et al.: lipids/proteins/water; membrane affinity "
                       "as the lipophilicity measure.",
    "Rodgers and Rowland": "adds drug ionization + electrostatic binding to acidic "
                           "phospholipids (good for moderate-strong bases).",
    "Schmitt": "water / neutral lipids / phospholipids / protein fractions with "
               "electrostatic interactions.",
    "Poulin and Theil": "foundational lipid/water partitioning method.",
    "Berezhkovskiy": "Poulin & Theil with a mass-balance correction.",
}
PERMEABILITY_METHOD_INFO = {
    "PK-Sim Standard": "default cellular permeability model.",
    "Charge dependent Schmitt": "charge-dependent permeability (Schmitt).",
}


# --------------------------------------------------------------------------- #
# addable process types (verified snapshot structure)
# --------------------------------------------------------------------------- #
PROCESS_TYPES: dict[str, dict[str, Any]] = {
    "metabolization_first_order": {
        "internal_name": "MetabolizationIntrinsic_FirstOrder",
        "data_source": "1st order CL",
        "applies_to": "enzyme",
        "parameters": [{"name": "Intrinsic clearance", "unit": "l/min",
                        "default": 0.1}],
        "description": "first-order metabolic clearance by an enzyme; adds a "
                       "-CLint*C_unbound term wherever that enzyme is expressed.",
    },
    "glomerular_filtration": {
        "internal_name": "GlomerularFiltration",
        "data_source": "GFR",
        "applies_to": "system",
        "parameters": [{"name": "GFR fraction", "unit": "", "default": 1.0}],
        "description": "renal clearance by glomerular filtration of unbound drug.",
    },
}


def addable_process_types(expressed_molecules: list[dict]) -> list[dict[str, Any]]:
    """The process types that can be legally added, with which molecules each
    can attach to (an enzyme process needs an expressed enzyme)."""
    enzymes = [m["molecule"] for m in expressed_molecules
               if (m.get("type") or "").lower() == "enzyme"]
    out = []
    for key, spec in PROCESS_TYPES.items():
        entry = {"type": key, "description": spec["description"],
                 "parameters": spec["parameters"]}
        if spec["applies_to"] == "enzyme":
            entry["can_attach_to"] = enzymes
        else:
            entry["can_attach_to"] = ["(system - no molecule needed)"]
        out.append(entry)
    return out
