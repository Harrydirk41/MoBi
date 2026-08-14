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
# tier: how the identification strategy MUST treat the parameter -
#   "constant"      -> a measured physical constant (MW, pKa, reference pH). It has
#                      one true value; fitting it means modeling a different
#                      molecule. NEVER estimated (hard rule).
#   "measured_soft" -> measured but genuinely refinable (fraction unbound,
#                      solubility). May be estimated ONLY when justified, and ONLY
#                      within its measured uncertainty range.
#   "estimate"      -> not reliably knowable a priori (effective lipophilicity,
#                      intrinsic clearance, permeabilities). Freely estimated.
PARAM_CATALOG: dict[str, dict[str, Any]] = {
    "Molecular weight": {
        "description": "molecular weight", "role": "measured", "tier": "constant"},
    "Solubility at reference pH": {
        "description": "aqueous solubility at the reference pH", "role": "measured",
        "tier": "measured_soft"},
    "Reference pH": {
        "description": "pH at which the solubility is defined", "role": "measured",
        "tier": "constant"},
    "Lipophilicity": {
        "description": "effective lipophilicity / membrane affinity (logP-like) "
        "driving tissue partitioning; the model value often differs from measured "
        "logD, so it is commonly refined",
        "range": [-2.0, 7.0], "role": "estimate", "tier": "estimate"},
    "Fraction unbound (plasma, reference value)": {
        "description": "unbound fraction in plasma; scales distribution and "
        "clearance. Measured, but small errors have large effect - sometimes refined",
        "range": [0.001, 1.0], "role": "measured", "tier": "measured_soft"},
    "Specific intestinal permeability (transcellular)": {
        "description": "transcellular intestinal permeability governing the rate "
        "and extent of oral absorption; in-vitro->in-vivo transfer is weak, so "
        "usually estimated",
        "range": [1e-8, 1e-2], "unit": "cm/min", "role": "estimate", "tier": "estimate"},
    "Permeability": {
        "description": "specific organ (cellular) permeability governing the "
        "kinetics of tissue distribution; usually estimated",
        "range": [1e-6, 1.0], "unit": "cm/min", "role": "estimate", "tier": "estimate"},
    "Intrinsic clearance": {
        "description": "first-order metabolic intrinsic clearance by the enzyme "
        "(per the process it belongs to); IVIVE unreliable, so usually estimated",
        "range": [1e-3, 5.0], "unit": "l/min", "role": "estimate", "tier": "estimate"},
    "GFR fraction": {
        "description": "fraction of glomerular filtration contributing to renal "
        "clearance of unbound drug (0 disables renal clearance)",
        "range": [0.0, 1.0], "role": "measured", "tier": "measured_soft"},
    "Vmax": {"description": "maximum rate of a Michaelis-Menten process",
             "range": [1e-3, 1e3], "role": "estimate", "tier": "estimate"},
    "Km": {"description": "Michaelis constant (half-saturation concentration)",
           "range": [1e-4, 1e3], "role": "estimate", "tier": "estimate"},
    "Plasma clearance": {"description": "lumped plasma clearance for a whole-organ "
                         "(hepatic/renal) clearance process - use when clearance is "
                         "not attributed to a specific enzyme",
                         "range": [1e-3, 100.0], "unit": "ml/min/kg", "role": "estimate",
                         "tier": "estimate"},
    # --- large molecules (proteins / monoclonal antibodies) --------------
    "Radius (solute)": {"description": "hydrodynamic radius of a large molecule "
                        "(protein/mAb); governs size-limited tissue permeation",
                        "range": [1e-3, 2e-2], "unit": "µm", "role": "measured",
                        "tier": "measured_soft"},
    "Kd (FcRn) in endosomal space": {"description": "FcRn binding affinity in the "
                        "endosome; drives antibody recycling and hence half-life "
                        "(large molecules only)",
                        "range": [1e-2, 1e2], "unit": "µmol/l", "role": "estimate",
                        "tier": "estimate"},
}


def describe_parameter(name: str) -> dict[str, Any]:
    return dict(PARAM_CATALOG.get(name, {}))


def param_role(name: str) -> str:
    return PARAM_CATALOG.get(name, {}).get("role", "unknown")


def param_tier(name: str) -> str:
    """Identification tier: 'constant' | 'measured_soft' | 'estimate'.

    Unknown parameters default to 'estimate'; anything that looks like a measured
    physical constant by name (molecular weight, pKa, reference pH) is a constant
    even if not catalogued, so it can never be fitted by accident."""
    entry = PARAM_CATALOG.get(name)
    if entry and entry.get("tier"):
        return entry["tier"]
    n = name.lower()
    if "molecular weight" in n or "pka" in n or "reference ph" in n:
        return "constant"
    return "estimate"


def measured_range(name: str, literature: list) -> tuple | None:
    """The measured-uncertainty range [lo, hi] for a 'measured_soft' parameter,
    taken from the study's literature values, so a refinement can be constrained
    to what the measurement actually supports. Returns None if no range is known."""
    key = name.lower()
    for e in literature or []:
        pn = str(e.get("parameter", "")).lower()
        if "unbound" in key and "unbound" in pn:
            rp = e.get("reported_range_percent")
            if rp:
                return (min(rp) / 100.0, max(rp) / 100.0)
            vals = ([float(e["value"])] if e.get("value") is not None else []) + \
                   [float(x) for x in e.get("reported_values", [])]
            if vals:
                return (min(vals), max(vals))
        if "solub" in key and "solub" in pn and e.get("value") is not None:
            v = float(e["value"])
            return (v / 2.0, v * 2.0)   # solubility rarely refined beyond ~2-fold
    return None


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
        # simulation-level mirror: {"Name": "<label>-<data_source>",
        #                           "SystemicProcessType": <systemic_type>}
        "systemic_label": "Glomerular Filtration", "systemic_type": "GFR",
        "parameters": [{"name": "GFR fraction", "unit": "", "default": 1.0}],
        "description": "renal clearance by glomerular filtration of unbound drug.",
    },
    "liver_clearance": {
        "internal_name": "LiverClearance", "data_source": "plasma clearance",
        "applies_to": "system", "validated": False,
        "systemic_label": "Liver Plasma Clearance", "systemic_type": "LiverClearance",
        "parameters": [{"name": "Plasma clearance", "unit": "ml/min/kg",
                        "default": 1.0}],
        "description": "lumped whole-liver (hepatic) plasma clearance - a simpler "
                       "alternative to enzyme-specific metabolism when you only "
                       "have total clearance.",
    },
    "kidney_clearance": {
        "internal_name": "KidneyClearance", "data_source": "plasma clearance",
        "applies_to": "system", "validated": False,
        "systemic_label": "Kidney Plasma Clearance", "systemic_type": "KidneyClearance",
        "parameters": [{"name": "Plasma clearance", "unit": "ml/min/kg",
                        "default": 1.0}],
        "description": "lumped renal plasma clearance (beyond passive GFR).",
    },
    "metabolization_specific_first_order": {
        "internal_name": "MetabolizationSpecific_FirstOrder",
        "data_source": "specific CL", "applies_to": "enzyme", "validated": False,
        "provenance": "OSP library snapshots",
        "parameters": [{"name": "Specific clearance", "unit": "1/min",
                        "default": 0.1}],
        "description": "first-order metabolism scaled by the enzyme's tissue "
                       "concentration (specific clearance).",
    },
    # --- additional kinetic variants (PK-Sim naming convention; VERIFY the exact
    #     parameter names/units once on your PK-Sim before benchmark use) --------
    "metabolization_hepatocytes_mm": {
        "internal_name": "MetabolizationHepatocytes_MM",
        "data_source": "hepatocytes", "applies_to": "enzyme", "validated": False,
        "provenance": "PK-Sim convention (parallels liver-microsomes MM); verify",
        "parameters": [
            {"name": "In vitro Vmax for hepatocytes",
             "unit": "pmol/min/10^6cells", "default": 100.0},
            {"name": "Km", "unit": "µmol/l", "default": 1.0},
            {"name": "kcat", "unit": "1/min", "default": 1.0}],
        "description": "saturable (Michaelis-Menten) metabolism scaled from an "
                       "in-vitro hepatocyte assay.",
    },
    "metabolization_specific_mm": {
        "internal_name": "MetabolizationSpecific_MM",
        "data_source": "specific MM", "applies_to": "enzyme", "validated": False,
        "provenance": "PK-Sim convention (specific MM); verify",
        "parameters": [
            {"name": "Vmax", "unit": "µmol/l/min", "default": 1.0},
            {"name": "Km", "unit": "µmol/l", "default": 1.0},
            {"name": "kcat", "unit": "1/min", "default": 1.0}],
        "description": "saturable metabolism scaled by the enzyme's tissue "
                       "concentration (specific Michaelis-Menten).",
    },
    "active_transport_first_order": {
        "internal_name": "ActiveTransportSpecific_FirstOrder",
        "data_source": "active transport 1st order", "applies_to": "transporter",
        "validated": False,
        "provenance": "PK-Sim convention (parallels active-transport MM); verify",
        "parameters": [{"name": "Transporter concentration", "unit": "µmol/l",
                        "default": 1.0},
                       {"name": "kcat", "unit": "1/min", "default": 1.0}],
        "description": "linear (non-saturable) active transport by a transporter - "
                       "use when the transporter is far from saturation.",
    },
    "biliary_clearance": {
        "internal_name": "BiliaryClearance", "data_source": "plasma clearance",
        "applies_to": "system", "validated": False,
        "systemic_label": "Biliary Clearance", "systemic_type": "BiliaryClearance",
        "provenance": "PK-Sim systemic process; verify InternalName/params",
        "parameters": [{"name": "Plasma clearance", "unit": "ml/min/kg",
                        "default": 1.0}],
        "description": "lumped biliary (hepatobiliary) plasma clearance into bile.",
    },
}


# --------------------------------------------------------------------------- #
# DDI / interaction mechanisms (perpetrator acting on a victim's enzyme).
# --------------------------------------------------------------------------- #
# These are part of PK-Sim's full process library, but they are NOT single-
# compound processes: an interaction links a PERPETRATOR compound to the enzyme
# of a VICTIM via the snapshot's ``Interactions`` block (a multi-compound DDI
# simulation). They are therefore DOCUMENTED here for completeness but are NOT
# offered through ``add_processes`` (which builds single-compound mechanisms);
# adding them requires a DDI setup with >1 Compound. ``validated`` marks the
# InternalNames confirmed against the Rifampicin-Digoxin-DDI library snapshot.
INTERACTION_PROCESS_TYPES: dict[str, dict[str, Any]] = {
    "competitive_inhibition": {
        "internal_name": "CompetitiveInhibition", "validated": True,
        "provenance": "Rifampicin-Digoxin-DDI snapshot",
        "parameters": [{"name": "Ki", "unit": "µmol/l"}],
        "description": "reversible competitive enzyme/transporter inhibition (Ki).",
    },
    "uncompetitive_inhibition": {
        "internal_name": "UncompetitiveInhibition", "validated": False,
        "provenance": "PK-Sim standard; verify InternalName",
        "parameters": [{"name": "Ki", "unit": "µmol/l"}],
        "description": "reversible uncompetitive inhibition (Ki).",
    },
    "noncompetitive_inhibition": {
        "internal_name": "NonCompetitiveInhibition", "validated": False,
        "provenance": "PK-Sim standard; verify InternalName",
        "parameters": [{"name": "Ki", "unit": "µmol/l"}],
        "description": "reversible noncompetitive inhibition (Ki).",
    },
    "mixed_inhibition": {
        "internal_name": "MixedInhibition", "validated": False,
        "provenance": "PK-Sim standard; verify InternalName/params",
        "parameters": [{"name": "Ki_c", "unit": "µmol/l"},
                       {"name": "Ki_u", "unit": "µmol/l"}],
        "description": "mixed competitive/uncompetitive inhibition (Ki_c, Ki_u).",
    },
    "mechanism_based_inhibition": {
        "internal_name": "IrreversibleInhibition", "validated": True,
        "provenance": "Rifampicin-Digoxin-DDI snapshot",
        "parameters": [{"name": "kinact", "unit": "1/min"},
                       {"name": "Ki", "unit": "µmol/l"}],
        "description": "irreversible / mechanism-based (time-dependent) inhibition "
                       "(kinact, Ki) - enzyme is inactivated.",
    },
    "induction": {
        "internal_name": "Induction", "validated": True,
        "provenance": "Rifampicin-Digoxin-DDI snapshot",
        "parameters": [{"name": "EC50", "unit": "µmol/l"},
                       {"name": "Emax", "unit": ""}],
        "description": "enzyme induction (EC50, Emax) - perpetrator up-regulates "
                       "the victim's enzyme, increasing its clearance.",
    },
}


def interaction_process_types() -> list[dict[str, Any]]:
    """The DDI/interaction mechanisms in PK-Sim's library. Documented for the
    full action space, but NOT addable to a single-compound model (they need a
    multi-compound DDI setup with an Interactions block)."""
    return [{"type": k, "internal_name": s["internal_name"],
             "description": s["description"], "validated": s.get("validated", False),
             "provenance": s.get("provenance"), "parameters": s["parameters"]}
            for k, s in INTERACTION_PROCESS_TYPES.items()]


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
                    "provenance": spec.get("provenance", "OSP library snapshots"),
                    "parameters": spec["parameters"], "can_attach_to": attach})
    return out
