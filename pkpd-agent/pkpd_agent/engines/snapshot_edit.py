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

    _apply_parameters(comp, edits.get("parameters") or {}, report)
    _apply_calc_methods(comp, edits.get("calculation_methods") or {}, report)
    _apply_processes(comp, edits.get("processes") or {}, report)
    return snap, report


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


def _apply_processes(comp: dict, processes: dict, report: dict) -> None:
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
    disabled = []
    for p in procs:
        remove = any((not enable) and matches(p, key)
                     for key, enable in processes.items())
        if remove:
            disabled.append(_process_label(p))
        else:
            keep.append(p)
    comp["Processes"] = keep

    for key, enable in processes.items():
        present = any(matches(p, key) for p in procs)
        if enable and not present:
            report["not_found"].append(
                f"process:{key} (enable not supported - add it in the snapshot)")
        elif not enable:
            report["processes"][key] = "disabled" if any(
                key.lower() in d.lower() or d.lower() in key.lower()
                for d in disabled) else "not_found"
    for d in disabled:
        report["processes"].setdefault(d, "disabled")
