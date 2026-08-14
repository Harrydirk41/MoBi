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
    "Plasma clearance": {"description": "lumped plasma clearance for a whole-organ "
                         "(hepatic/renal) clearance process - use when clearance is "
                         "not attributed to a specific enzyme",
                         "range": [1e-3, 100.0], "unit": "ml/min/kg", "role": "estimate"},
    # --- large molecules (proteins / monoclonal antibodies) --------------
    "Radius (solute)": {"description": "hydrodynamic radius of a large molecule "
                        "(protein/mAb); governs size-limited tissue permeation",
                        "range": [1e-3, 2e-2], "unit": "µm", "role": "measured"},
    "Kd (FcRn) in endosomal space": {"description": "FcRn binding affinity in the "
                        "endosome; drives antibody recycling and hence half-life "
                        "(large molecules only)",
                        "range": [1e-2, 1e2], "unit": "µmol/l", "role": "estimate"},
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
# Structures are the exact OSP snapshot forms taken from published library models
# (Alfentanil, Digoxin, Midazolam, Verapamil). ``validated`` = confirmed to run
# through PKSim.CLI here (first-order metabolism, GFR); the rest are structurally
# correct from published snapshots but should be validated once on your machine.
# ``applies_to``: enzyme | transporter | target | system.
PROCESS_TYPES: dict[str, dict[str, Any]] = {
    "metabolization_first_order": {
        "internal_name": "MetabolizationIntrinsic_FirstOrder",
        "data_source": "1st order CL", "applies_to": "enzyme", "validated": True,
        "parameters": [{"name": "Intrinsic clearance", "unit": "l/min",
                        "default": 0.1}],
        "description": "first-order (linear) metabolic clearance by an enzyme; "
                       "adds -CLint*C_unbound where the enzyme is expressed.",
    },
    "metabolization_mm": {
        "internal_name": "MetabolizationLiverMicrosomes_MM",
        "data_source": "liver microsomes", "applies_to": "enzyme",
        "validated": False,
        "parameters": [
            {"name": "In vitro Vmax for liver microsomes",
             "unit": "pmol/min/mg mic. protein", "default": 100.0},
            {"name": "Km", "unit": "µmol/l", "default": 1.0},
            {"name": "kcat", "unit": "1/min", "default": 1.0}],
        "description": "saturable (Michaelis-Menten) metabolism by an enzyme - "
                       "use for nonlinear/dose-dependent clearance.",
    },
    "active_transport_mm": {
        "internal_name": "ActiveTransportSpecific_MM",
        "data_source": "active transport", "applies_to": "transporter",
        "validated": False,
        "parameters": [
            {"name": "Transporter concentration", "unit": "µmol/l", "default": 1.0},
            {"name": "Vmax", "unit": "µmol/l/min", "default": 1.0},
            {"name": "Km", "unit": "µmol/l", "default": 1.0},
            {"name": "kcat", "unit": "1/min", "default": 1.0}],
        "description": "saturable active transport (efflux/uptake) by a "
                       "transporter, e.g. P-gp/ABCB1 efflux or OATP uptake.",
    },
    "specific_binding": {
        "internal_name": "SpecificBinding", "data_source": "binding",
        "applies_to": "target", "validated": False,
        "parameters": [{"name": "koff", "unit": "1/min", "default": 1.0},
                       {"name": "Kd", "unit": "nmol/l", "default": 1.0}],
        "description": "specific (target) binding - target-mediated disposition "
                       "(saturable binding to a tissue target).",
    },
    "glomerular_filtration": {
        "internal_name": "GlomerularFiltration", "data_source": "GFR",
        "applies_to": "system", "validated": True,
        "parameters": [{"name": "GFR fraction", "unit": "", "default": 1.0}],
        "description": "renal clearance by glomerular filtration of unbound drug.",
    },
    "liver_clearance": {
        "internal_name": "LiverClearance", "data_source": "plasma clearance",
        "applies_to": "system", "validated": False,
        "parameters": [{"name": "Plasma clearance", "unit": "ml/min/kg",
                        "default": 1.0}],
        "description": "lumped whole-liver (hepatic) plasma clearance - a simpler "
                       "alternative to enzyme-specific metabolism when you only "
                       "have total clearance.",
    },
    "kidney_clearance": {
        "internal_name": "KidneyClearance", "data_source": "plasma clearance",
        "applies_to": "system", "validated": False,
        "parameters": [{"name": "Plasma clearance", "unit": "ml/min/kg",
                        "default": 1.0}],
        "description": "lumped renal plasma clearance (beyond passive GFR).",
    },
    "metabolization_specific_first_order": {
        "internal_name": "MetabolizationSpecific_FirstOrder",
        "data_source": "specific CL", "applies_to": "enzyme", "validated": False,
        "parameters": [{"name": "Specific clearance", "unit": "1/min",
                        "default": 0.1}],
        "description": "first-order metabolism scaled by the enzyme's tissue "
                       "concentration (specific clearance).",
    },
}


def addable_process_types(expressed_molecules: list[dict]) -> list[dict[str, Any]]:
    """The process types that can be legally added, each with the molecules it
    can attach to (typed: metabolism->enzymes, transport->transporters,
    binding->any expressed molecule, GFR->system)."""
    by_type: dict[str, list[str]] = {"enzyme": [], "transporter": [], "any": []}
    for m in expressed_molecules:
        t = (m.get("type") or "").lower()
        by_type.setdefault(t, []).append(m["molecule"])
        by_type["any"].append(m["molecule"])
    out = []
    for key, spec in PROCESS_TYPES.items():
        at = spec["applies_to"]
        if at == "enzyme":
            attach = by_type.get("enzyme", [])
        elif at == "transporter":
            attach = by_type.get("transporter", [])
        elif at == "target":
            attach = by_type.get("any", [])
        else:
            attach = ["(system - no molecule needed)"]
        out.append({"type": key, "description": spec["description"],
                    "validated": spec.get("validated", False),
                    "parameters": spec["parameters"], "can_attach_to": attach})
    return out
