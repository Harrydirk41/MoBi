"""Reference library for LIBRARY-ASSISTED modeling.

The realistic pharma workflow: a modeler building a new compound reads the
finished models of ANALOGOUS compounds and reuses their structural choices. The
OSP library models are public and qualified, so giving the agent the OTHER
compounds' finished models is fair - it is the same library a human has. The one
inviolable rule is LEAVE-ONE-OUT: the target compound's own reference (values,
methods, held-out data) is never in the readable set, so held-out grading stays
honest.

Two digest levels:
  * structural (default): each other model's calculation methods, its process
    types (which enzymes/transporters), and the NAMES of its estimated
    parameters - teaches "how a CYP3A4 substrate is usually built" without
    handing over the fitted numbers (the agent must still fit them to the
    target's own data).
  * full: additionally includes the fitted parameter VALUES.
"""
from __future__ import annotations

import glob
import json
import os

_CLEAR_KW = ("kcat", "clint", "intrinsic clearance", "specific clearance",
             "clspec", "plasma clearance", "vmax")


def target_clearing_molecules(snapshot_path: str) -> list[str]:
    """The molecules (enzymes/transporters) the target's model attaches a
    mechanism to, read from the snapshot's processes. In normal mode the process
    STRUCTURE is given (only values are blanked), so this is given information,
    not an answer - used only to pick 'same-type' analogues. Empty in hard mode
    (structure stripped)."""
    try:
        with open(snapshot_path, encoding="utf-8") as fh:
            snap = json.load(fh)
    except (OSError, ValueError):
        return []
    mols: list[str] = []
    for comp in snap.get("Compounds") or []:
        for p in comp.get("Processes") or []:
            m = p.get("Molecule")
            if m and m not in mols:
                mols.append(m)
    return mols


def _digest(ak: dict, comp_dir: str, full: bool) -> tuple[dict, set]:
    """One reference model's digest + the set of molecules it involves."""
    sc = ak.get("structural_choices", {}) or {}
    methods = sc.get("calculation_methods", []) or []
    procs = sc.get("metabolizing_processes", []) or sc.get("processes", []) or []
    est = ak.get("estimated_parameters", []) or []
    proc_list = [{"molecule": p.get("molecule"), "type": p.get("type")}
                 for p in procs if isinstance(p, dict)]
    mols = {p["molecule"] for p in proc_list if p.get("molecule")}
    for e in est:
        pn = e.get("parameter", "")
        if "@" in pn:
            mols.add(pn.split("@", 1)[1])
    digest = {
        "compound": ak.get("compound", comp_dir),
        "calculation_methods": methods,
        "processes": proc_list,
        "estimated_parameters": (
            [{"parameter": e.get("parameter"), "value": e.get("value")} for e in est]
            if full else [e.get("parameter") for e in est]),
    }
    return digest, mols


def build_reference_library(target_dir: str, lib_dir: str, *, same_type: bool = False,
                            full: bool = False,
                            target_molecules: list[str] | None = None,
                            only: list[str] | None = None) -> list[dict]:
    """Digests of every OTHER compound's finished model (strict leave-one-out).

    `target_dir` is the target compound's directory NAME (e.g. 'Triazolam'); ALL
    of its variants (adult/pediatric) are excluded. `same_type` keeps only models
    sharing >=1 molecule with `target_molecules`. DDI answer keys are skipped."""
    tset = set(target_molecules or [])
    best: dict[str, dict] = {}                          # compound -> richest digest
    for akf in sorted(glob.glob(os.path.join(lib_dir, "*", "answer_key",
                                             "*.answer_key.json"))):
        if "ddi" in os.path.basename(akf):
            continue
        comp_dir = os.path.basename(os.path.dirname(os.path.dirname(akf)))
        if comp_dir == target_dir:                     # leave-one-out (all variants)
            continue
        if only is not None and comp_dir not in only:  # modeler restricted the set
            continue
        try:
            ak = json.load(open(akf, encoding="utf-8"))
        except (OSError, ValueError):
            continue
        digest, mols = _digest(ak, comp_dir, full)
        if same_type and tset and not (mols & tset):
            continue
        # collapse a compound's variants (adult/pediatric) to its richest digest
        key = digest["compound"]
        if key not in best or len(digest["estimated_parameters"]) > \
                len(best[key]["estimated_parameters"]):
            best[key] = digest
    return [best[k] for k in sorted(best)]


def library_for_snapshot(snapshot_path: str, *, same_type: bool = False,
                         full: bool = False, only: list[str] | None = None) -> dict:
    """Assemble the reference_library block to inject into a task input, deriving
    the library root and target compound from the snapshot path
    (.../OSP-PBPK-Model-Library/<Compound>/benchmark/<stem>.blanked.json)."""
    bench_dir = os.path.dirname(snapshot_path)             # <Compound>/benchmark
    comp_dir_path = os.path.dirname(bench_dir)             # <Compound>
    lib_dir = os.path.dirname(comp_dir_path)               # library root
    target_dir = os.path.basename(comp_dir_path)
    tmols = target_clearing_molecules(snapshot_path) if same_type else None
    models = build_reference_library(target_dir, lib_dir, same_type=same_type,
                                     full=full, target_molecules=tmols, only=only)
    note = (
        "REFERENCE LIBRARY (leave-one-out): finished OSP models for OTHER compounds - "
        "public, qualified models you may reuse as a human modeler would. The target "
        "compound's own reference is NOT here. Use these to choose the STRUCTURE "
        "(distribution/permeability method, which enzyme/transporter processes) by "
        "analogy to compounds with similar chemistry/clearance, then FIT the parameters "
        "to THIS compound's own data - do not copy their numbers blindly."
        + ("" if full else " Only their structural choices and parameter NAMES are shown; "
           "you must fit the values yourself."))
    return {"mode": ("same-type" if same_type else "all") + ("/full" if full else "/structural"),
            "note": note, "models": models}
