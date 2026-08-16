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
    # --- ionisation: compound type + pKa drive the charge-based partition methods
    # (Rodgers & Rowland, Schmitt) and the charged-fraction absorption penalty. A
    # compound is Neutral / Base / Acid with up to three pKa values. Measured
    # constants - never fitted - but the agent MUST see them to justify a
    # charge-aware method choice.
    "pKa": {
        "description": "acid/base dissociation constant. With the compound type "
        "(Neutral/Base/Acid; up to 3 pKa) it sets the ionised fraction at each pH - "
        "which drives the electrostatic terms in Rodgers & Rowland / Schmitt "
        "partitioning and the charged-fraction absorption penalty. Choose a "
        "charge-aware partition method when the drug is appreciably ionised",
        "range": [0.0, 14.0], "role": "measured", "tier": "constant"},
    "Compound type": {
        "description": "protonation class: Neutral / Base / Acid (may be zwitterion "
        "with multiple pKa). Determines the sign of the electrostatic tissue-binding "
        "term - bases bind acidic phospholipids (favours Rodgers & Rowland)",
        "role": "measured", "tier": "constant"},
    "Lipophilicity": {
        "description": "effective lipophilicity / membrane affinity (logP-like) "
        "driving tissue partitioning; the model value often differs from measured "
        "logD, so it is commonly refined",
        "range": [-2.0, 7.0], "role": "estimate", "tier": "estimate"},
    "Fraction unbound (plasma, reference value)": {
        "description": "unbound fraction in plasma; scales distribution and "
        "clearance. Measured, but small errors have large effect - sometimes "
        "refined. The binding partner (albumin / alpha-1-acid glycoprotein / "
        "unknown) sets how fu scales with age (ontogeny)",
        "range": [0.001, 1.0], "role": "measured", "tier": "measured_soft"},
    "Specific intestinal permeability (paracellular)": {
        "description": "paracellular (between-cell) intestinal permeability - a "
        "second absorption route beside the transcellular one, relevant for small "
        "hydrophilic/charged molecules",
        "range": [1e-9, 1e-3], "unit": "cm/min", "role": "estimate", "tier": "estimate"},
    "Hill coefficient": {
        "description": "cooperativity exponent of a Hill (sigmoidal saturable) "
        "process; n=1 reduces to Michaelis-Menten, n>1 = positive cooperativity",
        "range": [0.5, 4.0], "role": "estimate", "tier": "estimate"},
    "In vitro half-life (hepatocytes)": {
        "description": "in-vitro depletion half-life in a hepatocyte assay; scaled "
        "to whole-liver clearance. IVIVE unreliable, so usually refined",
        "range": [0.1, 1e4], "unit": "min", "role": "estimate", "tier": "estimate"},
    "In vitro half-life (microsomes)": {
        "description": "in-vitro depletion half-life in a liver-microsome assay; "
        "scaled to whole-liver clearance. Usually refined",
        "range": [0.1, 1e4], "unit": "min", "role": "estimate", "tier": "estimate"},
    "Cell concentration": {
        "description": "hepatocyte density used to scale an in-vitro half-life to "
        "the whole liver; a measured assay condition, not normally fitted",
        "range": [1e-3, 1e3], "unit": "10^6cells/ml", "role": "measured",
        "tier": "measured_soft"},
    "Microsomal protein content": {
        "description": "microsomal protein concentration used to scale in-vitro "
        "kinetics to the liver; a measured scaling input, not normally fitted",
        "range": [1e-3, 1e3], "unit": "mg/ml", "role": "measured", "tier": "measured_soft"},
    "In vitro Vmax": {
        "description": "in-vitro maximal transport/metabolic rate (e.g. from a "
        "vesicular assay); scaled in vivo, usually refined",
        "range": [1e-4, 1e4], "role": "estimate", "tier": "estimate"},
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
    # --- in-vitro enzyme-kinetic inputs (Michaelis-Menten / recombinant) -----
    # measured in vitro but IVIVE (in-vitro->in-vivo) transfer is unreliable, so
    # the effective in-vivo value is usually refined against the clinical data.
    "In vitro Vmax for liver microsomes": {
        "description": "in-vitro maximal metabolic rate measured in liver "
        "microsomes (Vmax) for a Michaelis-Menten enzyme; scaled to the whole "
        "liver in vivo. IVIVE is unreliable, so often refined",
        "range": [1e-4, 1e4], "role": "estimate", "tier": "estimate"},
    "In vitro Vmax for hepatocytes": {
        "description": "in-vitro maximal metabolic rate measured in a hepatocyte "
        "assay (Vmax); scaled to the whole liver in vivo, usually refined",
        "range": [1e-4, 1e4], "role": "estimate", "tier": "estimate"},
    "In vitro Vmax/recombinant enzyme": {
        "description": "in-vitro maximal metabolic rate per recombinant enzyme "
        "(Vmax) for a Michaelis-Menten pathway; scaled by enzyme abundance in "
        "vivo. Often refined against the clinical data",
        "range": [1e-4, 1e4], "role": "estimate", "tier": "estimate"},
    "In vitro CL/recombinant enzyme": {
        "description": "in-vitro intrinsic clearance per recombinant enzyme "
        "(first-order CLint) for a specific enzyme; scaled by enzyme abundance in "
        "vivo. IVIVE unreliable, so usually refined",
        "range": [1e-4, 1e4], "role": "estimate", "tier": "estimate"},
    "In vitro CL for liver microsomes": {
        "description": "in-vitro intrinsic clearance measured in liver microsomes "
        "(first-order CLint); scaled to the whole liver in vivo, usually refined",
        "range": [1e-4, 1e4], "unit": "µl/min/mg mic. protein",
        "role": "estimate", "tier": "estimate"},
    "Vmax (liver tissue)": {
        "description": "maximal metabolic rate expressed per whole-liver tissue "
        "mass (intrinsic Michaelis-Menten Vmax); usually refined",
        "range": [1e-4, 1e4], "unit": "µmol/min/kg tissue",
        "role": "estimate", "tier": "estimate"},
    "kcat": {"description": "catalytic turnover number of an enzyme (Vmax per unit "
             "enzyme) for a Michaelis-Menten pathway; the enzyme-normalised rate, "
             "scaled by enzyme abundance in vivo. Usually refined",
             "range": [1e-4, 1e6], "unit": "1/min", "role": "estimate", "tier": "estimate"},
    "Content of CYP proteins in liver microsomes": {
        "description": "measured enzyme abundance (pmol enzyme per mg microsomal "
        "protein) used to scale in-vitro kinetics to the whole liver; a measured "
        "scaling input, not normally fitted",
        "range": [1e-2, 1e3], "role": "measured", "tier": "measured_soft"},
    "Transporter concentration": {
        "description": "relative abundance of a membrane transporter driving active "
        "uptake/efflux; scales the transport rate. A structural/measured input",
        "range": [1e-3, 1e3], "role": "measured", "tier": "measured_soft"},
    # --- binding kinetics (specific/target binding) --------------------------
    "koff": {"description": "dissociation rate constant of specific (target/tissue) "
             "binding; with Kd sets the on-rate. Governs slow-binding kinetics "
             "(e.g. digoxin tissue binding). Often refined",
             "range": [1e-6, 1e3], "unit": "1/min", "role": "estimate", "tier": "estimate"},
    "Kd": {"description": "equilibrium dissociation constant of specific (target/"
           "tissue) binding - binding affinity (lower = tighter). Often refined "
           "against the data",
           "range": [1e-6, 1e3], "unit": "µmol/l", "role": "estimate", "tier": "estimate"},
    # --- experimental physicochemical variants -------------------------------
    "Lipophilicity (experiment)": {
        "description": "experimentally measured lipophilicity (logP/logD) as an "
        "input; the effective in-vivo value used for partitioning may differ and "
        "is sometimes refined",
        "range": [-2.0, 7.0], "role": "measured", "tier": "measured_soft"},
    "Fraction unbound (experiment)": {
        "description": "experimentally measured unbound fraction (plasma or the "
        "stated matrix); scales distribution and clearance. Measured, small errors "
        "have large effect",
        "range": [1e-4, 1.0], "role": "measured", "tier": "measured_soft"},
    "F": {"description": "bioavailable fraction of an oral dose that reaches the "
          "systemic circulation; usually an OUTCOME of absorption + first-pass "
          "rather than a free knob - set only if the model uses it directly",
          "range": [0.0, 1.0], "role": "estimate", "tier": "estimate"},
    "Cl": {"description": "generic first-order clearance parameter on a process; "
           "prefer the process-specific clearance (Intrinsic clearance, "
           "CLspec/[Enzyme], Plasma clearance) when present",
           "range": [1e-4, 1e3], "role": "estimate", "tier": "estimate"},
    "Solubility table": {
        "description": "tabulated aqueous solubility versus pH (a measured "
        "solubility profile) rather than a single reference value; a given input",
        "range": None, "role": "measured", "tier": "measured_soft"},
    "Plasma clearance": {"description": "lumped plasma clearance for a whole-organ "
                         "(hepatic/renal) clearance process - use when clearance is "
                         "not attributed to a specific enzyme",
                         "range": [1e-3, 100.0], "unit": "ml/min/kg", "role": "estimate",
                         "tier": "estimate"},
    # --- specific-clearance (enzyme-kinetic) family ----------------------
    # On a MetabolizationSpecific process the INPUT you fit is CLspec/[Enzyme];
    # 'Specific clearance' and 'Enzyme concentration' are its derived/structural
    # siblings, not independent knobs - fit CLspec/[Enzyme], not those.
    "CLspec/[Enzyme]": {
        "description": "specific intrinsic clearance normalised to enzyme "
        "concentration (CLint per unit enzyme) for a specific-clearance "
        "metabolization process; THIS is the fittable clearance for such a "
        "process. With several enzymes on one compound each is targeted by the "
        "qualified name 'CLspec/[Enzyme]@<Molecule>'. Plasma parent data pin only "
        "total clearance, so the split across enzymes is weakly identifiable "
        "without in-vitro CLint",
        "range": [1e-3, 100.0], "unit": "1/min", "role": "estimate", "tier": "estimate"},
    "Specific clearance": {
        "description": "DERIVED product CLspec/[Enzyme] x enzyme concentration - a "
        "structural output of the specific-clearance process, NOT an independent "
        "knob; fit CLspec/[Enzyme] instead (this usually reads 0 until built)",
        "range": [0.0, 1e3], "unit": "1/min", "role": "derived", "tier": "estimate"},
    "Enzyme concentration": {
        "description": "enzyme abundance the specific clearance is scaled by; a "
        "structural input set by the expression profile, not a free fitting knob - "
        "fit CLspec/[Enzyme] instead",
        "range": [1e-6, 1e3], "unit": "µmol/l", "role": "derived", "tier": "estimate"},
    # --- DDI interaction-mechanism parameters (perpetrator on victim enzyme) ---
    "Ki": {"description": "reversible inhibition constant - perpetrator "
           "concentration for half-maximal enzyme/transporter inhibition (lower = "
           "stronger inhibitor)",
           "range": [1e-4, 1e3], "unit": "µmol/l", "role": "estimate", "tier": "estimate"},
    "Ki_c": {"description": "competitive component of mixed inhibition (Ki for the "
             "competitive term); trades off with Ki_u from a single ratio",
             "range": [1e-4, 1e3], "unit": "µmol/l", "role": "estimate", "tier": "estimate"},
    "Ki_u": {"description": "uncompetitive component of mixed inhibition; trades off "
             "with Ki_c - fix one, fit the other from a single interaction ratio",
             "range": [1e-4, 1e3], "unit": "µmol/l", "role": "estimate", "tier": "estimate"},
    "kinact": {"description": "maximal inactivation rate of mechanism-based "
               "(irreversible) inhibition; with K_kinact_half only their ratio "
               "(kinact/K_kinact_half, the inactivation efficiency) is identified "
               "from a single interaction ratio",
               "range": [1e-4, 10.0], "unit": "1/min", "role": "estimate", "tier": "estimate"},
    "K_kinact_half": {"description": "perpetrator concentration for half-maximal "
                      "inactivation in mechanism-based inhibition; trades off with "
                      "kinact - fix to its in-vitro value and fit kinact",
                      "range": [1e-3, 1e3], "unit": "µmol/l", "role": "estimate", "tier": "estimate"},
    "EC50": {"description": "perpetrator concentration for half-maximal enzyme "
             "induction; trades off with Emax unless arms span [I] far below and "
             "above it - usually fixed to its in-vitro value",
             "range": [1e-3, 1e3], "unit": "µmol/l", "role": "estimate", "tier": "estimate"},
    "Emax": {"description": "maximal fold-induction of the target enzyme by the "
             "perpetrator (dimensionless, added to baseline 1); trades off with "
             "EC50 from a single perpetrator dose",
             "range": [0.0, 50.0], "unit": "", "role": "estimate", "tier": "estimate"},
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
    "Charge dependent Schmitt": "charge-dependent permeability (Schmitt) - "
                                "penalises permeation of the ionised fraction.",
    "Charge dependent Schmitt normalized to PK-Sim": "charge-dependent Schmitt "
                                "rescaled to the PK-Sim Standard baseline.",
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
    # --- recombinant-CYP and intrinsic/microsomal variants VERIFIED against the
    #     OSP library snapshots (exact InternalName + parameter names read from
    #     real models) -------------------------------------------------------- #
    "metabolization_recombinant_mm": {
        "internal_name": "rCYP450_MM", "data_source": "recombinant enzyme",
        "applies_to": "enzyme", "validated": False, "internal_name_verified": True,
        "provenance": "verified: Sildenafil/Efavirenz/Felodipine/Fluvoxamine/"
                      "Itraconazole/Erythromycin snapshots",
        "parameters": [
            {"name": "In vitro Vmax/recombinant enzyme",
             "unit": "pmol/min/pmol rec. enzyme", "default": 1.0},
            {"name": "Km", "unit": "µmol/l", "default": 1.0},
            {"name": "kcat", "unit": "1/min", "default": 1.0}],
        "description": "saturable (Michaelis-Menten) metabolism from a RECOMBINANT-"
                       "enzyme (rCYP) assay, scaled by enzyme abundance - the "
                       "standard CYP metabolism form in the library (e.g. "
                       "sildenafil's CYP3A4/2C9/2C19).",
    },
    "metabolization_recombinant_first_order": {
        "internal_name": "rCYP450_FirstOrder", "data_source": "recombinant enzyme",
        "applies_to": "enzyme", "validated": False, "internal_name_verified": True,
        "provenance": "verified: Montelukast/Ethinylestradiol snapshots",
        "parameters": [
            {"name": "In vitro CL/recombinant enzyme",
             "unit": "µl/min/pmol rec. enzyme", "default": 1.0}],
        "description": "first-order (linear) metabolism from a RECOMBINANT-enzyme "
                       "(rCYP) intrinsic-clearance assay, scaled by enzyme abundance.",
    },
    "metabolization_intrinsic_mm": {
        "internal_name": "MetabolizationIntrinsic_MM", "data_source": "intrinsic MM",
        "applies_to": "enzyme", "validated": False, "internal_name_verified": True,
        "provenance": "verified: Moclobemide snapshot",
        "parameters": [
            {"name": "Vmax (liver tissue)", "unit": "µmol/min/kg tissue", "default": 1.0},
            {"name": "Km", "unit": "µmol/l", "default": 1.0}],
        "description": "saturable (Michaelis-Menten) metabolism expressed as a "
                       "whole-liver-tissue Vmax (intrinsic), not scaled per enzyme.",
    },
    "metabolization_microsomes_first_order": {
        "internal_name": "MetabolizationLiverMicrosomes_FirstOrder",
        "data_source": "liver microsomes", "applies_to": "enzyme",
        "validated": False, "internal_name_verified": True,
        "provenance": "verified: S-Mephenytoin snapshot",
        "parameters": [
            {"name": "In vitro CL for liver microsomes",
             "unit": "µl/min/mg mic. protein", "default": 1.0},
            {"name": "Content of CYP proteins in liver microsomes",
             "unit": "pmol/mg mic. protein", "default": 1.0}],
        "description": "first-order (linear) metabolism from a liver-microsome "
                       "intrinsic-clearance assay, scaled to the whole liver.",
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
    # NOTE: PK-Sim has NO first-order/linear transporter process - transporters
    # are Michaelis-Menten or Hill only (verified against the PK-Sim v12 docs), so
    # a linear active-transport type was removed as spurious.
    "active_transport_hill": {
        "internal_name": "ActiveTransportSpecific_Hill",
        "data_source": "active transport Hill", "applies_to": "transporter",
        "validated": False, "internal_name_verified": False,
        "provenance": "PK-Sim v12 docs (Specific active transport - Hill); "
                      "InternalName inferred - run validate_processes.py to confirm",
        "parameters": [
            {"name": "Transporter concentration", "unit": "µmol/l", "default": 1.0},
            {"name": "Vmax", "unit": "µmol/l/min", "default": 1.0},
            {"name": "Km", "unit": "µmol/l", "default": 1.0},
            {"name": "Hill coefficient", "unit": "", "default": 1.0}],
        "description": "cooperative (Hill) saturable active transport - use when "
                       "transport shows sigmoidal, cooperative saturation.",
    },
    "active_transport_vesicular_mm": {
        "internal_name": "ActiveTransportVesicular_MM",
        "data_source": "vesicular assay", "applies_to": "transporter",
        "validated": False, "internal_name_verified": False,
        "provenance": "PK-Sim v12 docs (In vitro active transport, vesicular assay - "
                      "MM); InternalName inferred - validate before use",
        "parameters": [
            {"name": "In vitro Vmax", "unit": "pmol/min/mg protein", "default": 1.0},
            {"name": "Km", "unit": "µmol/l", "default": 1.0}],
        "description": "saturable active transport scaled from an in-vitro vesicular "
                       "(membrane-vesicle) transport assay.",
    },
    "metabolization_hill": {
        "internal_name": "MetabolizationSpecific_Hill",
        "data_source": "specific Hill", "applies_to": "enzyme",
        "validated": False, "internal_name_verified": False,
        "provenance": "PK-Sim v12 docs (In vitro clearance - Hill); InternalName "
                      "inferred - run validate_processes.py to confirm",
        "parameters": [
            {"name": "Vmax", "unit": "µmol/l/min", "default": 1.0},
            {"name": "Km", "unit": "µmol/l", "default": 1.0},
            {"name": "Hill coefficient", "unit": "", "default": 1.0}],
        "description": "cooperative (Hill) saturable metabolism - sigmoidal "
                       "concentration dependence, use over MM when cooperativity "
                       "is evident.",
    },
    "tubular_secretion_first_order": {
        "internal_name": "TubularSecretion_FirstOrder",
        "data_source": "tubular secretion", "applies_to": "transporter",
        "validated": False, "internal_name_verified": False,
        "provenance": "PK-Sim v12 docs (Tubular Secretion - First Order); "
                      "InternalName inferred - validate before use",
        "parameters": [{"name": "Transporter concentration", "unit": "µmol/l",
                        "default": 1.0},
                       {"name": "Specific clearance", "unit": "1/min", "default": 0.1}],
        "description": "active renal TUBULAR SECRETION (first order) by a kidney "
                       "transporter - renal clearance BEYOND passive GFR. Use when "
                       "CL_renal exceeds fu x GFR.",
    },
    "tubular_secretion_mm": {
        "internal_name": "TubularSecretion_MM",
        "data_source": "tubular secretion MM", "applies_to": "transporter",
        "validated": False, "internal_name_verified": False,
        "provenance": "PK-Sim v12 docs (Tubular Secretion - Michaelis-Menten); "
                      "InternalName inferred - validate before use",
        "parameters": [
            {"name": "Transporter concentration", "unit": "µmol/l", "default": 1.0},
            {"name": "Vmax", "unit": "µmol/l/min", "default": 1.0},
            {"name": "Km", "unit": "µmol/l", "default": 1.0}],
        "description": "saturable (Michaelis-Menten) active renal tubular secretion "
                       "beyond passive GFR.",
    },
    "hepatic_clearance_hepatocytes_thalf": {
        "internal_name": "MetabolizationHepatocytes_tHalf",
        "data_source": "hepatocytes t1/2", "applies_to": "system",
        "validated": False, "internal_name_verified": False,
        "systemic_label": "Total Hepatic Clearance", "systemic_type": "LiverClearance",
        "provenance": "PK-Sim v12 docs (In vitro hepatocytes - t1/2, Total Hepatic "
                      "Clearance); InternalName inferred - validate before use",
        "parameters": [{"name": "In vitro half-life (hepatocytes)", "unit": "min",
                        "default": 30.0},
                       {"name": "Cell concentration", "unit": "10^6cells/ml",
                        "default": 1.0}],
        "description": "whole-liver clearance scaled from an in-vitro hepatocyte "
                       "depletion HALF-LIFE - a systemic hepatic clearance you use "
                       "when you have a t1/2, not enzyme kinetics.",
    },
    "hepatic_clearance_microsomes_thalf": {
        "internal_name": "MetabolizationLiverMicrosomes_tHalf",
        "data_source": "microsomes t1/2", "applies_to": "system",
        "validated": False, "internal_name_verified": False,
        "systemic_label": "Total Hepatic Clearance", "systemic_type": "LiverClearance",
        "provenance": "PK-Sim v12 docs (In vitro liver microsomes - t1/2, Total "
                      "Hepatic Clearance); InternalName inferred - validate first",
        "parameters": [{"name": "In vitro half-life (microsomes)", "unit": "min",
                        "default": 30.0},
                       {"name": "Microsomal protein content", "unit": "mg/ml",
                        "default": 1.0}],
        "description": "whole-liver clearance scaled from an in-vitro liver-"
                       "microsome depletion HALF-LIFE.",
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
# InternalNames and parameter names below are VERIFIED against the OSP library
# DDI snapshots (Rifampicin, Erythromycin, Fluconazole, Atazanavir, ...).
# ``validated`` still means "confirmed to run end-to-end through PKSim.CLI here",
# which requires the multi-compound DDI build (not yet exercised), so it stays
# False; ``internal_name_verified`` marks the name/params as ground-truth-correct.
INTERACTION_PROCESS_TYPES: dict[str, dict[str, Any]] = {
    "competitive_inhibition": {
        "internal_name": "CompetitiveInhibition", "validated": False,
        "internal_name_verified": True,
        "provenance": "verified: Rifampicin/Itraconazole/Cimetidine snapshots",
        "parameters": [{"name": "Ki", "unit": "µmol/l"}],
        "description": "reversible competitive enzyme/transporter inhibition (Ki).",
    },
    "uncompetitive_inhibition": {
        "internal_name": "UncompetitiveInhibition", "validated": False,
        "internal_name_verified": False,
        "provenance": "PK-Sim standard; not seen in the library DDI set - verify",
        "parameters": [{"name": "Ki", "unit": "µmol/l"}],
        "description": "reversible uncompetitive inhibition (Ki).",
    },
    "noncompetitive_inhibition": {
        "internal_name": "NoncompetitiveInhibition", "validated": False,
        "internal_name_verified": True,
        "provenance": "verified: Fluconazole snapshot",
        "parameters": [{"name": "Ki", "unit": "µmol/l"}],
        "description": "reversible noncompetitive inhibition (Ki).",
    },
    "mixed_inhibition": {
        "internal_name": "MixedInhibition", "validated": False,
        "internal_name_verified": True,
        "provenance": "verified: Atazanavir/Mefenamic_acid snapshots",
        "parameters": [{"name": "Ki_c", "unit": "µmol/l"},
                       {"name": "Ki_u", "unit": "µmol/l"}],
        "description": "mixed competitive/uncompetitive inhibition (Ki_c, Ki_u).",
        # a single AUCR at one victim/inhibitor level cannot separate the
        # competitive (Ki_c) from the uncompetitive (Ki_u) component.
        "identifiability": {
            "tradeoff_pair": ["Ki_c", "Ki_u"],
            "identifiable_combination": "the net inhibition at the studied "
                "victim/inhibitor concentration",
            "fix_first": "Ki_u",
            "note": "Ki_c and Ki_u trade off; from one AUCR fit Ki_c and fix "
                    "Ki_u to its in-vitro value (or set it high to reduce to "
                    "competitive inhibition).",
        },
    },
    "mechanism_based_inhibition": {
        "internal_name": "IrreversibleInhibition", "validated": False,
        "internal_name_verified": True,
        "provenance": "verified: Erythromycin/Clarithromycin/Moclobemide snapshots",
        "parameters": [{"name": "kinact", "unit": "1/min"},
                       {"name": "K_kinact_half", "unit": "µmol/l"}],
        "description": "irreversible / mechanism-based (time-dependent) inhibition "
                       "(kinact, K_kinact_half) - the enzyme is inactivated.",
        # from a single AUCR only the inactivation EFFICIENCY kinact/K_kinact_half
        # (the first-order inactivation rate when [I] << K_half) is identified;
        # the two parameters individually trade off.
        "identifiability": {
            "tradeoff_pair": ["kinact", "K_kinact_half"],
            "identifiable_combination": "kinact / K_kinact_half (inactivation "
                "efficiency)",
            "fix_first": "K_kinact_half",
            "note": "kinact and K_kinact_half trade off; a single AUCR pins only "
                    "their ratio. Fix K_kinact_half to its in-vitro value and fit "
                    "kinact (or report the ratio, not the split).",
        },
    },
    "induction": {
        "internal_name": "Induction", "validated": False,
        "internal_name_verified": True,
        "provenance": "verified: Rifampicin/Carbamazepine/Efavirenz snapshots",
        "parameters": [{"name": "EC50", "unit": "µmol/l"},
                       {"name": "Emax", "unit": ""}],
        "description": "enzyme induction (EC50, Emax) - perpetrator up-regulates "
                       "the target enzyme, increasing its clearance.",
        # from one AUCR only the induction magnitude at the studied [I],
        # Emax*[I]/(EC50+[I]), is identified; EC50 and Emax individually trade off
        # unless arms span [I] both well below and well above EC50.
        "identifiability": {
            "tradeoff_pair": ["EC50", "Emax"],
            "identifiable_combination": "Emax*[I]/(EC50+[I]) at the studied "
                "perpetrator exposure",
            "fix_first": "EC50",
            "note": "EC50 and Emax trade off from a single perpetrator dose. Fix "
                    "EC50 to its in-vitro value and fit Emax; only arms spanning "
                    "[I] far below and far above EC50 identify both.",
        },
    },
}


def interaction_by_internal_name(name: str) -> dict[str, Any]:
    """Look up a DDI mechanism's catalog entry by its PK-Sim InternalName
    (e.g. 'IrreversibleInhibition'). Returns {} if unknown."""
    for entry in INTERACTION_PROCESS_TYPES.values():
        if entry.get("internal_name") == name:
            return entry
    return {}


def interaction_process_types() -> list[dict[str, Any]]:
    """The DDI/interaction mechanisms in PK-Sim's library. Documented for the
    full action space, but NOT addable to a single-compound model (they need a
    multi-compound DDI setup with an Interactions block)."""
    return [{"type": k, "internal_name": s["internal_name"],
             "description": s["description"], "validated": s.get("validated", False),
             "internal_name_verified": s.get("internal_name_verified", False),
             "provenance": s.get("provenance"), "parameters": s["parameters"],
             # a 2-parameter mechanism trades off from a single interaction ratio;
             # tell the agent up front which to fix and what is identifiable
             **({"identifiability": s["identifiability"]}
                if s.get("identifiability") else {})}
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
