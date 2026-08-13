"""Build a clean, agent-ready INPUT file from an OSP snapshot JSON.

An OSP (PK-Sim / MoBi) snapshot mixes three very different things in one file:

  1. the INPUTS you are handed  -> observed clinical data, literature
     physicochemical/in-vitro properties, study designs, demographics
  2. the MODELING STRATEGY       -> distribution/permeability methods, which
     parameters were estimated, the parameter-identification setup
  3. the RESULTS                 -> the fitted parameter values and the
     simulated concentration-time output

For a modeling benchmark you want to hand an agent only (1), keep (2)/(3) as a
hidden answer key, and ask "can you rebuild a model that reproduces the observed
data?". This script extracts (1) into a tidy JSON under ``json_input/``.

What it keeps (the clean input):
  * objective + compound identity
  * given physicochemical / in-vitro data  (ValueOrigin = Publication/measured)
  * the list of parameters that WERE estimated -- names + units only, values
    hidden -- so the agent knows what is unknown without seeing the answer
  * clinical observed concentration-time profiles (time in h, conc in mg/L,
    plus the native values and the arithmetic SD when reported) with full
    dosing metadata (study, route, dose, matrix)
  * study designs (protocols: route, dose, formulation)
  * demographics (individuals / populations: species, age, weight, gender)

What it drops (goes to the answer key, not the input):
  * fitted parameter VALUES
  * distribution / cellular-permeability calculation methods
  * parameter-identification configuration and simulation results

Usage:
    python build_clean_input.py json/Alfentanil-Model.json
    python build_clean_input.py json/Alfentanil-Pediatrics.json --objective "..."
    python build_clean_input.py <snapshot.json> --outdir json_input

Writes ``json_input/<snapshot-stem>.input.json`` (stdlib only, runs anywhere).
"""

from __future__ import annotations

import argparse
import json
import os
from typing import Any


# --------------------------------------------------------------------------- #
# unit helpers
# --------------------------------------------------------------------------- #

def _num(v: Any) -> float | None:
    """Coerce to float; missing markers ('', 'NaN', None, lists) -> None."""
    if isinstance(v, bool) or isinstance(v, (list, dict)):
        return None
    if v is None:
        return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return None if f != f else f                          # drop NaN


def _to_mg_L(values: list, unit: str, mol_weight: float | None) -> list:
    """Convert a concentration column to mg/L (non-numeric entries -> None)."""
    nums = [_num(v) for v in values]
    u = (unit or "").strip().lower()
    mass = {"mg/l": 1.0, "µg/l": 1e-3, "ug/l": 1e-3, "ng/ml": 1e-3,
            "ng/l": 1e-6, "g/l": 1e3}
    if u in mass:
        f = mass[u]
    elif "mol" in u and mol_weight:                      # molar -> mass
        f = mol_weight / 1000.0                           # µmol/L * g/mol /1000 = mg/L
        if u.startswith("nmol"):
            f /= 1000.0
    else:
        return nums                                       # unknown unit: pass through
    return [None if v is None else v * f for v in nums]


def _time_to_h(values: list, unit: str) -> list:
    nums = [_num(v) for v in values]
    u = (unit or "").strip().lower()
    if u.startswith("min"):
        return [None if v is None else v / 60.0 for v in nums]
    if u.startswith("day"):
        return [None if v is None else v * 24.0 for v in nums]
    return nums


# --------------------------------------------------------------------------- #
# extractors
# --------------------------------------------------------------------------- #

def _ext_props(od: dict[str, Any]) -> dict[str, Any]:
    return {e.get("Name"): e.get("Value") for e in od.get("ExtendedProperties") or []}


def clinical_observed(data: dict[str, Any]) -> list[dict[str, Any]]:
    """Tidy every observed profile: metadata + time(h) + conc(mg/L) + SD."""
    out = []
    for od in data.get("ObservedData") or []:
        base = od.get("BaseGrid") or {}
        time_h = _time_to_h(list(base.get("Values") or []), base.get("Unit", ""))
        cols = od.get("Columns") or []
        if not cols or not time_h:
            continue
        col = cols[0]
        mw = (col.get("DataInfo") or {}).get("MolWeight")
        conc_native = list(col.get("Values") or [])
        conc_mg_L = _to_mg_L(conc_native, col.get("Unit", ""), mw)

        # arithmetic SD, if a related column carries it (already reported in mg/L)
        sd_mg_L = None
        for rc in col.get("RelatedColumns") or []:
            if (rc.get("DataInfo") or {}).get("AuxiliaryType") == "ArithmeticStdDev":
                sd_mg_L = _to_mg_L(list(rc.get("Values") or []), rc.get("Unit", ""), mw)
                break

        ep = _ext_props(od)
        n = min(len(time_h), len(conc_native))
        rec = {
            "dataset": od.get("Name"),
            "study": ep.get("Study Id"),
            "molecule": ep.get("Molecule"),
            "route": ep.get("Route"),
            "dose": ep.get("Dose"),
            "matrix": "/".join(x for x in (ep.get("Organ"), ep.get("Compartment")) if x),
            "patient_group": ep.get("Patient Id"),
            "n_points": n,
            "time_unit": "h",
            "conc_unit": "mg/L",
            "mol_weight_g_per_mol": mw,
            "time_h": [None if t is None else round(t, 5) for t in time_h[:n]],
            "conc_mg_L": [None if c is None else round(c, 8) for c in conc_mg_L[:n]],
            "conc_native": conc_native[:n],
            "conc_native_unit": col.get("Unit"),
        }
        if sd_mg_L is not None:
            rec["sd_mg_L"] = [None if s is None else round(s, 8) for s in sd_mg_L[:n]]
        out.append(rec)
    return out


def given_physchem_and_unknowns(data: dict[str, Any]):
    """Split compound parameters into GIVEN (measured/published) inputs and the
    list of parameters that were ESTIMATED (names/units only -- values hidden)."""
    comp = (data.get("Compounds") or [{}])[0]
    given: list[dict[str, Any]] = []
    to_identify: list[dict[str, Any]] = []
    seen: set[tuple] = set()

    def walk(obj: Any, alt: str = "") -> None:
        if isinstance(obj, dict):
            if "Name" in obj and isinstance(obj.get("Value"), (int, float)):
                vo = obj.get("ValueOrigin") or {}
                src = vo.get("Source")
                key = (obj["Name"], obj.get("Unit"))
                if key not in seen:
                    seen.add(key)
                    if src == "ParameterIdentification":
                        to_identify.append({
                            "parameter": obj["Name"],
                            "unit": obj.get("Unit", ""),
                            "note": "estimated during model building (value withheld)",
                        })
                    else:
                        given.append({
                            "parameter": obj["Name"],
                            "value": obj["Value"],
                            "unit": obj.get("Unit", ""),
                            "source": vo.get("Description") or src or "snapshot default",
                        })
            for k, v in obj.items():
                walk(v, obj.get("Name", alt) if k == "Parameters" else alt)
        elif isinstance(obj, list):
            for v in obj:
                walk(v, alt)

    walk(comp)

    # pKa (structural physchem, always a given input)
    pka = []
    for p in comp.get("PkaTypes") or []:
        pka.append({"type": p.get("Type"), "pka": p.get("Pka"),
                    "source": (p.get("ValueOrigin") or {}).get("Description")})

    meta = {
        "name": comp.get("Name"),
        "is_small_molecule": comp.get("IsSmallMolecule"),
        "plasma_protein_binding_partner": comp.get("PlasmaProteinBindingPartner"),
        "pka": pka,
    }
    return meta, given, to_identify


def study_designs(data: dict[str, Any]) -> list[dict[str, Any]]:
    out = []
    for pr in data.get("Protocols") or []:
        params = {p.get("Name"): (p.get("Value"), p.get("Unit"))
                  for p in pr.get("Parameters") or []}
        dose, dunit = params.get("InputDose", (None, None))
        start, sunit = params.get("Start time", (None, None))
        out.append({
            "name": pr.get("Name"),
            "application_type": pr.get("ApplicationType"),
            "dosing_interval": pr.get("DosingInterval"),
            "dose": dose, "dose_unit": dunit,
            "start_time": start, "start_time_unit": sunit,
        })
    return out


def demographics(data: dict[str, Any]) -> list[dict[str, Any]]:
    out = []
    for ind in data.get("Individuals") or []:
        od = ind.get("OriginData") or {}
        age = od.get("Age") or {}
        wt = od.get("Weight") or {}
        out.append({
            "name": ind.get("Name"),
            "kind": "individual",
            "species": od.get("Species"),
            "population": od.get("Population"),
            "gender": od.get("Gender"),
            "age": age.get("Value"), "age_unit": age.get("Unit"),
            "weight": wt.get("Value"), "weight_unit": wt.get("Unit"),
            "source": (od.get("ValueOrigin") or {}).get("Description"),
        })
    for pop in data.get("Populations") or []:
        settings = pop.get("Settings") or {}
        out.append({
            "name": pop.get("Name"),
            "kind": "population",
            "number_of_individuals": settings.get("NumberOfIndividuals")
            or pop.get("NumberOfIndividuals"),
            "details": "see snapshot Populations block for full demographic ranges",
        })
    return out


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #

def build_clean_input(path: str, objective: str | None = None) -> dict[str, Any]:
    with open(path, encoding="utf-8") as fh:
        data = json.load(fh)

    meta, given, to_identify = given_physchem_and_unknowns(data)
    obs = clinical_observed(data)

    routes = sorted({o["route"] for o in obs if o.get("route")})
    studies = sorted({o["study"] for o in obs if o.get("study")})

    return {
        "schema": "osp-clean-input/v1",
        "source_snapshot": os.path.basename(path),
        "compound": meta["name"],
        "objective": objective or (
            f"Build a PBPK model for {meta['name']} that reproduces the observed "
            f"plasma concentration-time data below, using the given "
            f"physicochemical/in-vitro inputs and estimating the listed unknown "
            f"parameters. Modeling strategy and fitted values are withheld."),
        "data_overview": {
            "n_observed_datasets": len(obs),
            "routes": routes,
            "studies": studies,
            "n_given_physchem_params": len(given),
            "n_parameters_to_identify": len(to_identify),
            "n_study_designs": len(data.get("Protocols") or []),
            "n_demographics": len(data.get("Individuals") or [])
            + len(data.get("Populations") or []),
        },
        "compound_identity": meta,
        "given_physicochemical": given,
        "parameters_to_identify": to_identify,
        "clinical_observed_data": obs,
        "study_designs": study_designs(data),
        "demographics": demographics(data),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("snapshot", help="path to the OSP snapshot .json")
    ap.add_argument("--outdir", default="json_input",
                    help="output folder (default: json_input)")
    ap.add_argument("--objective", default=None,
                    help="override the default objective string")
    args = ap.parse_args()

    clean = build_clean_input(args.snapshot, args.objective)
    os.makedirs(args.outdir, exist_ok=True)
    stem = os.path.splitext(os.path.basename(args.snapshot))[0]
    out = os.path.join(args.outdir, f"{stem}.input.json")
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(clean, fh, indent=2, ensure_ascii=False)

    ov = clean["data_overview"]
    print(f"wrote {out}")
    print(f"  compound            : {clean['compound']}")
    print(f"  observed datasets   : {ov['n_observed_datasets']}  "
          f"routes={ov['routes']}")
    print(f"  given physchem      : {ov['n_given_physchem_params']}")
    print(f"  params to identify  : {ov['n_parameters_to_identify']}")
    print(f"  study designs       : {ov['n_study_designs']}")
    print(f"  demographics        : {ov['n_demographics']}")


if __name__ == "__main__":
    main()
