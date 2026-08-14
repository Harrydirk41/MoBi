"""Apply structured model edits to an OSP snapshot (the agent's action space).

The snapshot is the editable *source* of a PBPK model. This module applies a
declarative edit spec so the agent can change both parameters and structure and
then have PK-Sim recompile+run the result:

    edits = {
      "parameters": {                      # compound parameter values
          "Lipophilicity": 2.0,
          "Intrinsic clearance": 0.4,
      },
      "calculation_methods": {             # distribution / permeability choice
          "partition": "Rodgers and Rowland",     # one of PARTITION_METHODS
          "permeability": "PK-Sim Standard",      # one of PERMEABILITY_METHODS
      },
      "processes": {                       # enable/disable a process
          "CYP3A4": True,                  # False -> remove the enzyme process
          "GFR": True,
      },
    }

Every edit is reported back (applied / not-found) so a wrong name never fails
silently. Method-name validity is ultimately decided by PK-Sim; the lists below
are guidance and are not hard-enforced.
"""

from __future__ import annotations

import copy
from typing import Any

from .osp_catalog import PROCESS_TYPES

# valid PK-Sim distribution / permeability methods (snapshot short names)
PARTITION_METHODS = [
    "PK-Sim Standard", "Rodgers and Rowland", "Schmitt",
    "Poulin and Theil", "Berezhkovskiy",
]
PERMEABILITY_METHODS = ["PK-Sim Standard", "Charge dependent Schmitt"]

_PARTITION_PREFIX = "Cellular partition coefficient method - "
_PERMEABILITY_PREFIX = "Cellular permeability - "


def apply_edits(snapshot: dict, edits: dict | None) -> tuple[dict, dict]:
    """Return (new_snapshot, report). Does not mutate the input."""
    snap = copy.deepcopy(snapshot)
    report: dict[str, Any] = {"parameters": {}, "calculation_methods": {},
                              "processes": {}, "not_found": []}
    if not edits:
        return snap, report
    comp = (snap.get("Compounds") or [{}])[0]
    comp_name = comp.get("Name")
    sims = snap.get("Simulations") or []
    expressed = {ep.get("Molecule"): (ep.get("Type") or "").lower()
                 for ep in (snap.get("ExpressionProfiles") or [])
                 if ep.get("Molecule")}

    _apply_parameters(comp, edits.get("parameters") or {}, report)
    _apply_calc_methods(comp, edits.get("calculation_methods") or {}, report)
    _apply_processes(comp, sims, comp_name, edits.get("processes") or {}, report)
    _apply_add_processes(comp, sims, comp_name,
                         edits.get("add_processes") or [], expressed, report)
    return snap, report


# --------------------------------------------------------------------------- #
# simulation-level process mirroring
# --------------------------------------------------------------------------- #
# A process lives in TWO places in a PK-Sim snapshot: the compound
# (Compounds[i].Processes, the definition) and *every* simulation that uses the
# compound (Simulations[n].Compounds[j].Processes, a reference carrying that
# simulation's parameter values). The simulation reference links back by
# MoleculeName (enzymes/transporters) or SystemicProcessType (GFR & lumped
# clearances). If we remove/add a process on the compound but not on the
# simulations, PK-Sim's snapshot->project mapper hits a dangling (or missing)
# reference and MapToModel throws ("no .pksim5 produced by snap"). These helpers
# keep the two in sync.

def _sim_compounds(sims: list, comp_name: str | None):
    """Yield each simulation compound-block that is the compound we edited."""
    for sim in sims:
        if not isinstance(sim, dict):
            continue
        for sc in sim.get("Compounds") or []:
            if not isinstance(sc, dict):
                continue
            if comp_name is None or sc.get("Name") is None \
                    or sc.get("Name") == comp_name:
                yield sc


def _sim_proc_matches(sim_proc: dict, comp_proc: dict) -> bool:
    """Does this simulation-level process reference the given compound process?"""
    mol = comp_proc.get("Molecule")
    if mol:                                   # enzyme / transporter
        return sim_proc.get("MoleculeName") == mol
    ds = (comp_proc.get("DataSource") or "").lower()          # systemic (GFR, ...)
    spt = (sim_proc.get("SystemicProcessType") or "")
    if spt and ds and spt.lower() == ds:
        return True
    nm = (sim_proc.get("Name") or "").lower()
    return bool(ds) and nm.endswith("-" + ds)


def _remove_sim_processes(sims: list, comp_name: str | None,
                          removed: list[dict]) -> int:
    """Drop simulation-level entries that reference any removed compound process."""
    n = 0
    for sc in _sim_compounds(sims, comp_name):
        procs = sc.get("Processes")
        if not isinstance(procs, list):
            continue
        keep = [sp for sp in procs
                if not any(_sim_proc_matches(sp, cp) for cp in removed)]
        n += len(procs) - len(keep)
        sc["Processes"] = keep
    return n


def _add_sim_process(sims: list, comp_name: str | None,
                     spec: dict, mol: str | None) -> int:
    """Mirror a newly-added compound process onto every simulation's compound."""
    if mol:                                   # enzyme / transporter reference
        entry = {"Name": f"{mol}-{spec['data_source']}", "MoleculeName": mol}
        def dup(sp): return sp.get("MoleculeName") == mol \
            and (sp.get("Name") or "").endswith("-" + spec["data_source"])
    else:                                     # systemic reference
        label = spec.get("systemic_label") or spec["internal_name"]
        styp = spec.get("systemic_type") or spec["data_source"]
        entry = {"Name": f"{label}-{spec['data_source']}",
                 "SystemicProcessType": styp}
        def dup(sp): return (sp.get("SystemicProcessType") or "") == styp
    n = 0
    for sc in _sim_compounds(sims, comp_name):
        procs = sc.setdefault("Processes", [])
        if any(dup(sp) for sp in procs if isinstance(sp, dict)):
            continue
        procs.append(dict(entry))
        n += 1
    return n


# --------------------------------------------------------------------------- #
# parameters
# --------------------------------------------------------------------------- #

def _apply_parameters(comp: dict, params: dict, report: dict) -> None:
    want = {k.lower(): (k, v) for k, v in params.items()}
    hit: set[str] = set()

    def walk(obj: Any) -> None:
        if isinstance(obj, dict):
            nm = obj.get("Name")
            if isinstance(nm, str) and nm.lower() in want and "Value" in obj:
                orig, val = want[nm.lower()]
                obj["Value"] = val
                # a user-set value is no longer a fit result / default
                obj.pop("ValueOrigin", None)
                report["parameters"][nm] = val
                hit.add(nm.lower())
            for v in obj.values():
                walk(v)
        elif isinstance(obj, list):
            for v in obj:
                walk(v)

    walk(comp)
    for low, (orig, _) in want.items():
        if low not in hit:
            report["not_found"].append(f"parameter:{orig}")


# --------------------------------------------------------------------------- #
# calculation methods (distribution / permeability)
# --------------------------------------------------------------------------- #

def _apply_calc_methods(comp: dict, methods: dict, report: dict) -> None:
    cms = comp.get("CalculationMethods")
    if cms is None:
        cms = comp["CalculationMethods"] = []

    def set_method(prefix: str, name: str) -> None:
        full = name if " - " in name else prefix + name
        # replace an existing entry with the same prefix, else append
        for i, m in enumerate(cms):
            if isinstance(m, str) and m.startswith(prefix):
                cms[i] = full
                return
        cms.append(full)

    for key, name in methods.items():
        k = key.lower()
        if k in ("partition", "partition_coefficient", "distribution"):
            set_method(_PARTITION_PREFIX, name)
            report["calculation_methods"]["partition"] = name
        elif k in ("permeability", "cellular_permeability"):
            set_method(_PERMEABILITY_PREFIX, name)
            report["calculation_methods"]["permeability"] = name
        else:
            report["not_found"].append(f"calculation_method:{key}")


# --------------------------------------------------------------------------- #
# processes (enable / disable an enzyme / transporter / clearance)
# --------------------------------------------------------------------------- #

def _process_label(p: dict) -> str:
    return (p.get("Molecule") or p.get("InternalName") or p.get("DataSource")
            or p.get("Name") or "process")


def _apply_processes(comp: dict, sims: list, comp_name: str | None,
                     processes: dict, report: dict) -> None:
    procs = comp.get("Processes") or []

    def matches(p: dict, key: str) -> bool:
        k = key.lower()
        for field in ("Molecule", "InternalName", "DataSource", "Name"):
            v = p.get(field)
            if isinstance(v, str) and k in v.lower():
                return True
        # allow "GFR" / "glomerular" to hit the renal process
        if k in ("gfr", "renal", "glomerular"):
            return any(isinstance(p.get(f), str) and "glomerul" in p[f].lower()
                       for f in ("InternalName", "DataSource"))
        return False

    keep = []
    removed_procs = []
    disabled = []
    keys_that_removed: set[str] = set()
    for p in procs:
        hit_keys = [key for key, enable in processes.items()
                    if (not enable) and matches(p, key)]
        if hit_keys:
            removed_procs.append(p)
            disabled.append(_process_label(p))
            keys_that_removed.update(hit_keys)
        else:
            keep.append(p)
    comp["Processes"] = keep
    # keep the simulations in sync so PK-Sim's snapshot mapper does not choke on
    # a dangling process reference.
    n_sim = _remove_sim_processes(sims, comp_name, removed_procs) \
        if removed_procs else 0

    for key, enable in processes.items():
        present = any(matches(p, key) for p in procs)
        if enable and not present:
            report["not_found"].append(
                f"process:{key} (use add_processes to add a new mechanism)")
        elif not enable:
            report["processes"][key] = (
                "disabled" if key in keys_that_removed else "not_found")
    for d in disabled:
        report["processes"].setdefault(d, "disabled")
    if n_sim:
        report["processes"]["_simulation_refs_removed"] = n_sim


def _apply_add_processes(comp: dict, sims: list, comp_name: str | None,
                         additions: list, expressed: set, report: dict) -> None:
    """Attach a NEW process (mechanism) to the compound. An enzyme process needs
    an already-expressed enzyme; the process structure follows the catalog."""
    procs = comp.setdefault("Processes", [])
    for a in additions:
        if not isinstance(a, dict):
            continue
        typ = a.get("type")
        mol = a.get("molecule")
        spec = PROCESS_TYPES.get(typ)
        if not spec:
            report["not_found"].append(f"process_type:{typ}")
            continue
        at = spec["applies_to"]
        if at in ("enzyme", "transporter", "target"):
            if not mol:
                report["not_found"].append(f"add_process:{typ} needs a molecule")
                continue
            moltype = expressed.get(mol)
            if moltype is None:
                report["not_found"].append(
                    f"add_process:{mol} is not expressed in the model")
                continue
            if at == "enzyme" and moltype != "enzyme":
                report["not_found"].append(
                    f"add_process:{mol} is a {moltype}, not an enzyme")
                continue
            if at == "transporter" and moltype != "transporter":
                report["not_found"].append(
                    f"add_process:{mol} is a {moltype}, not a transporter")
                continue
        # already there? skip
        if any(p.get("InternalName") == spec["internal_name"]
               and p.get("Molecule") == mol for p in procs):
            report["processes"][f"add:{mol or typ}"] = "already_present"
            continue
        block: dict[str, Any] = {"InternalName": spec["internal_name"],
                                 "DataSource": spec["data_source"],
                                 "Species": "Human"}
        if mol:
            block["Molecule"] = mol
        given = a.get("parameters") or {}
        block["Parameters"] = [
            {"Name": p["name"], "Value": float(given.get(p["name"], p["default"])),
             "Unit": p["unit"]}
            for p in spec["parameters"]]
        procs.append(block)
        # mirror onto the simulations, or the added process is inert (a compound
        # process that no simulation references is never built into the model).
        n_sim = _add_sim_process(sims, comp_name, spec, mol)
        report["processes"][f"add:{mol or typ}"] = (
            "added (no matching simulation compound)"
            if sims and not n_sim else "added")
