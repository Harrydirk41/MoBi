"""Generic benchmark-file generator for the OSP model library.

Follows the Alfentanil pattern but reads each snapshot's OWN metadata instead of
a hand-curated card, so it generalises across compounds without inventing text:

  * fitted parameters are those with ValueOrigin.Source == 'ParameterIdentification'
    (the values a modeler identified) - these become the ANSWER KEY and are
    blanked to a documented, no-leak type default (never a function of the true
    value) to make the agent's starting model;
  * given physicochemical values are those with ValueOrigin.Source == 'Publication'
    (or 'In Vitro'); their ValueOrigin.Description carries the literature source,
    so literature_physicochemical is extracted from the snapshot, not fabricated;
  * known biology is derived from the model's own structure (active processes and
    calculation methods) - factual, from the snapshot.

For each single-compound small-molecule snapshot it writes:
    json_input/<stem>.input.json          (agent-facing; no answers)
    benchmark/<stem>.blanked.json         (agent start point; fitted params reset)
    answer_key/<stem>.answer_edits.json   (fitted values; grading only)
    answer_key/<stem>.answer_key.json     (fitted params + structure; grading only)

Audit mode classifies every snapshot and flags the ones that need category-
specific handling (DDI, metabolite, biologic, multi-compound) rather than
silently producing a wrong single-compound benchmark.

    python build_benchmark.py --audit                 # classify everything, write nothing
    python build_benchmark.py --audit <dir-or-file>
    python build_benchmark.py <Compound>/json/<Stem>.json      # generate one
"""

from __future__ import annotations

import argparse
import glob
import json
import os
from typing import Any

# --------------------------------------------------------------------------- #
# no-leak blanking defaults: a documented naive start value per parameter,
# INDEPENDENT of the fitted value (chosen from general PBPK priors / catalog
# ranges), so the blanked model encodes no information about the answer.
# --------------------------------------------------------------------------- #
BLANK_DEFAULTS: dict[str, float] = {
    "Lipophilicity": 2.0,
    "Fraction unbound (plasma, reference value)": 0.5,
    "Specific intestinal permeability (transcellular)": 1.0e-5,
    "Permeability": 1.0e-4,
    "Intrinsic clearance": 0.1,
    "Specific clearance": 0.1,
    "Plasma clearance": 1.0,
    "GFR fraction": 1.0,
    "Vmax": 1.0, "Km": 1.0, "kcat": 1.0,
    "Kd": 1.0, "koff": 1.0,
    "Transporter concentration": 1.0,
    "Ki": 1.0, "K_kinact_half": 1.0, "kinact": 0.1,
    "EC50": 1.0, "Emax": 1.0,
}

INTERACTION_NAMES = {"CompetitiveInhibition", "UncompetitiveInhibition",
                     "NoncompetitiveInhibition", "MixedInhibition",
                     "IrreversibleInhibition", "Induction"}


def blank_default(name: str) -> float:
    """A documented, no-leak starting value for a fitted parameter - a value a
    modeler would plausibly start from, chosen from the parameter TYPE only
    (never from the fitted value). Falls back by keyword so an unusual/templated
    name (e.g. 'CLspec/[Enzyme]', 'In vitro Vmax for liver microsomes') still
    gets a sensible prior rather than blocking generation."""
    if name in BLANK_DEFAULTS:
        return BLANK_DEFAULTS[name]
    n = name.lower()
    if "vmax" in n and "microsom" in n:
        return 100.0
    if "vmax" in n:
        return 1.0
    if "km" in n:
        return 1.0
    if "clearance" in n or "clspec" in n or n.startswith("cl"):
        return 0.1
    if "permeab" in n:
        return 1.0e-4
    if "solub" in n:
        return 100.0
    if "lipophil" in n:
        return 2.0
    if "unbound" in n or n.startswith("fu"):
        return 0.5
    if any(k in n for k in ("kd", "ki", "koff", "kcat", "kinact", "ec50", "emax")):
        return 1.0
    return 1.0                       # generic positive prior


def _num(v: Any) -> float | None:
    if isinstance(v, bool) or isinstance(v, (list, dict)) or v is None:
        return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return None if f != f else f


def _to_mg_L(values: list, unit: str, mw: float | None) -> list:
    nums = [_num(v) for v in values]
    u = (unit or "").strip().lower()
    mass = {"mg/l": 1.0, "µg/l": 1e-3, "ug/l": 1e-3, "ng/ml": 1e-3,
            "ng/l": 1e-6, "g/l": 1e3}
    if u in mass:
        f = mass[u]
    elif "mol" in u and mw:
        f = mw / 1000.0
        if u.startswith("nmol"):
            f /= 1000.0
    else:
        return nums
    return [None if v is None else v * f for v in nums]


def _time_to_h(values: list, unit: str) -> list:
    nums = [_num(v) for v in values]
    u = (unit or "").strip().lower()
    if u.startswith("min"):
        return [None if v is None else v / 60.0 for v in nums]
    if u.startswith("day"):
        return [None if v is None else v * 24.0 for v in nums]
    return nums


def _ext_props(od: dict) -> dict:
    return {e.get("Name"): e.get("Value") for e in od.get("ExtendedProperties") or []}


# --------------------------------------------------------------------------- #
# extractors (agent-facing, from the snapshot)
# --------------------------------------------------------------------------- #

def clinical_observed(data: dict) -> list[dict]:
    out = []
    for od in data.get("ObservedData") or []:
        base = od.get("BaseGrid") or {}
        time_h = _time_to_h(list(base.get("Values") or []), base.get("Unit", ""))
        cols = od.get("Columns") or []
        if not cols or not time_h:
            continue
        col = cols[0]
        mw = (col.get("DataInfo") or {}).get("MolWeight")
        conc_mg_L = _to_mg_L(list(col.get("Values") or []), col.get("Unit", ""), mw)
        sd_mg_L = None
        for rc in col.get("RelatedColumns") or []:
            if (rc.get("DataInfo") or {}).get("AuxiliaryType") == "ArithmeticStdDev":
                sd_mg_L = _to_mg_L(list(rc.get("Values") or []), rc.get("Unit", ""), mw)
                break
        ep = _ext_props(od)
        n = min(len(time_h), len(conc_mg_L))
        rec = {"dataset": od.get("Name"), "study": ep.get("Study Id"),
               "molecule": ep.get("Molecule"), "route": ep.get("Route"),
               "dose": ep.get("Dose"), "n_points": n,
               "time_unit": "h", "conc_unit": "mg/L",
               "time_h": [None if t is None else round(t, 5) for t in time_h[:n]],
               "conc_mg_L": [None if c is None else round(c, 8) for c in conc_mg_L[:n]]}
        if sd_mg_L is not None:
            rec["sd_mg_L"] = [None if s is None else round(s, 8) for s in sd_mg_L[:n]]
        out.append(rec)
    return out


def literature_physicochemical(comp: dict) -> list[dict]:
    """Given (measured/published) physchem values, from the compound's own
    ValueOrigin - the source citation is in ValueOrigin.Description."""
    out, seen = [], set()

    def walk(o):
        if isinstance(o, dict):
            nm = o.get("Name")
            vo = o.get("ValueOrigin") or {}
            if isinstance(nm, str) and isinstance(o.get("Value"), (int, float)) \
                    and vo.get("Source") in ("Publication", "In Vitro", "Database", "Internet") \
                    and nm not in seen:
                seen.add(nm)
                out.append({"parameter": nm, "value": o["Value"],
                            "unit": o.get("Unit", ""),
                            "source": vo.get("Description")})
            for v in o.values():
                walk(v)
        elif isinstance(o, list):
            for v in o:
                walk(v)

    walk(comp)
    return out


def _walk_fitted(node):
    """(name, value, unit) for every ParameterIdentification value under node."""
    out = []

    def w(o):
        if isinstance(o, dict):
            nm = o.get("Name")
            if isinstance(nm, str) and isinstance(o.get("Value"), (int, float)) \
                    and (o.get("ValueOrigin") or {}).get("Source") == "ParameterIdentification":
                out.append((nm, o["Value"], o.get("Unit", "")))
            for v in o.values():
                w(v)
        elif isinstance(o, list):
            for v in o:
                w(v)

    w(node)
    return out


def fitted_answer(comp: dict) -> dict:
    """{key: value} of the parameters to recover. A key is QUALIFIED as
    ``<Name>@<Molecule>`` when the same parameter name is fitted on more than one
    process (e.g. a per-enzyme CLspec on UGT1A9, UGT2B7, a CYP) - so the distinct
    per-enzyme values are kept, not collapsed. Compound-level params and
    single-occurrence process params keep their plain name."""
    from collections import defaultdict
    comp_no_proc = {k: v for k, v in comp.items() if k != "Processes"}
    answer = {n: v for (n, v, _u) in _walk_fitted(comp_no_proc)}   # compound-level
    proc_occ = defaultdict(list)
    for p in comp.get("Processes") or []:
        mol = p.get("Molecule") or p.get("InternalName")
        for (n, v, _u) in _walk_fitted({"Parameters": p.get("Parameters") or []}):
            proc_occ[n].append((mol, v))
    for n, occs in proc_occ.items():
        if len(occs) == 1:
            answer[n] = occs[0][1]
        else:                                   # collision -> qualify by molecule
            for mol, v in occs:
                answer[f"{n}@{mol}"] = v
    return answer


def fitted_parameters(comp: dict) -> list[dict]:
    """Flat list of fitted parameters (qualified where names collide)."""
    return [{"parameter": k, "value": v} for k, v in fitted_answer(comp).items()]


def known_biology(comp: dict) -> list[str]:
    """Factual structure notes derived from the model itself."""
    facts = []
    procs = comp.get("Processes") or []
    enz = sorted({p.get("Molecule") for p in procs
                  if p.get("Molecule") and "Metaboli" in (p.get("InternalName") or "")})
    trans = sorted({p.get("Molecule") for p in procs
                    if p.get("Molecule") and "Transport" in (p.get("InternalName") or "")})
    if enz:
        facts.append(f"Metabolized by: {', '.join(enz)}.")
    if trans:
        facts.append(f"Transported by: {', '.join(trans)}.")
    if any("Glomerul" in (p.get("InternalName") or "") for p in procs):
        facts.append("Renal clearance via glomerular filtration is present.")
    methods = [m for m in comp.get("CalculationMethods") or []]
    if methods:
        facts.append("Reference distribution/permeability methods exist but are "
                     "withheld; choose your own.")
    return facts


# --------------------------------------------------------------------------- #
# classification / audit
# --------------------------------------------------------------------------- #

def classify(data: dict) -> dict:
    comps = data.get("Compounds") or []
    small = all(c.get("IsSmallMolecule", True) for c in comps)
    has_inter = any(s.get("Interactions") for s in data.get("Simulations") or [])
    perps = [c.get("Name") for c in comps
             if any(p.get("InternalName") in INTERACTION_NAMES
                    for p in c.get("Processes") or [])]
    names = [c.get("Name") for c in comps]
    is_metab = len(comps) > 1 and any(
        any(m in (n or "").lower() for m in ("hydroxy", "keto", "desalkyl",
            "[15n]", "n-des", "eso", "r-omeprazole"))
        for n in names)
    if not small:
        kind = "biologic"
    elif is_metab:
        kind = "metabolite"
    elif has_inter and len(comps) > 1:
        kind = "ddi-victim"
    elif has_inter:
        kind = "ddi-auto"
    else:
        kind = "single"
    comp0 = comps[0] if comps else {}
    fitted = fitted_parameters(comp0)
    return {"kind": kind, "n_compounds": len(comps), "compounds": names,
            "perpetrators": perps, "n_fitted": len(fitted),
            "fitted": [f["parameter"] for f in fitted],
            "n_observed": len(clinical_observed(data)),
            "missing_defaults": [f["parameter"] for f in fitted
                                 if f["parameter"] not in BLANK_DEFAULTS]}


# --------------------------------------------------------------------------- #
# generation (single-compound small molecules)
# --------------------------------------------------------------------------- #

def _blank(data: dict) -> tuple[dict, list[dict], list[str]]:
    """Return (blanked_snapshot, answer_edits, warnings). Resets every fitted
    parameter (across ALL compounds) to its no-leak default."""
    import copy
    snap = copy.deepcopy(data)
    answer, warns = [], []

    def walk(o):
        if isinstance(o, dict):
            nm = o.get("Name")
            vo = o.get("ValueOrigin") or {}
            if isinstance(nm, str) and isinstance(o.get("Value"), (int, float)) \
                    and vo.get("Source") == "ParameterIdentification":
                answer.append({"parameter": nm, "value": o["Value"],
                               "unit": o.get("Unit", "")})
                if nm not in BLANK_DEFAULTS:
                    warns.append(nm)          # used a keyword fallback, not an exact default
                o["Value"] = blank_default(nm)
                o["ValueOrigin"] = {"Source": "Unknown",
                                    "Description": "benchmark naive prior (blanked)"}
            for v in o.values():
                walk(v)
        elif isinstance(o, list):
            for v in o:
                walk(v)

    # blank only the COMPOUND (drug) parameters, matching the Alfentanil design -
    # NOT formulation dissolution or protocol parameters, which are a different,
    # study-specific kind of fit and are described in the given study designs.
    for comp in snap.get("Compounds") or []:
        walk(comp)
    return snap, answer, warns


def build_files(path: str) -> dict:
    with open(path, encoding="utf-8") as fh:
        data = json.load(fh)
    stem = os.path.splitext(os.path.basename(path))[0]
    info = classify(data)
    if info["kind"] != "single":
        return {"skip": True, "stem": stem, "reason": info["kind"], "info": info}
    if info["n_fitted"] == 0:
        return {"skip": True, "stem": stem,
                "reason": "no fitted parameters (scaling/extrapolation, not a "
                          "parameter-recovery task)", "info": info}
    if info["n_observed"] == 0:
        return {"skip": True, "stem": stem,
                "reason": "no observed data to fit against", "info": info}

    comp = (data.get("Compounds") or [{}])[0]
    obs = clinical_observed(data)
    lit = literature_physicochemical(comp)
    fitted = fitted_parameters(comp)

    agent_input = {
        "schema": "osp-agent-task/v2",
        "task_id": stem,
        "source_snapshot": os.path.basename(path),
        "compound": comp.get("Name"),
        "objective": (f"Build a whole-body PBPK model for {comp.get('Name')} that "
                      "reproduces the observed plasma concentration-time data. "
                      "Decide the model structure and any parameters you cannot "
                      "obtain from the given data, and justify them. No modeling "
                      "strategy and no reference results are provided."),
        "background": {"description": f"{comp.get('Name')} PBPK modeling task.",
                       "literature_facts": known_biology(comp)},
        "given_data": {
            "compound_identity": {"name": comp.get("Name"),
                                  "is_small_molecule": comp.get("IsSmallMolecule")},
            "literature_physicochemical": lit,
            "clinical_observed_data": obs,
            "study_designs": [{"name": pr.get("Name"),
                               "application_type": pr.get("ApplicationType")}
                              for pr in data.get("Protocols") or []],
            "demographics": [{"name": ind.get("Name"),
                              "species": (ind.get("OriginData") or {}).get("Species")}
                             for ind in data.get("Individuals") or []],
        },
        "data_overview": {"n_observed_datasets": len(obs),
                          "routes": sorted({o["route"] for o in obs if o.get("route")}),
                          "studies": sorted({o["study"] for o in obs if o.get("study")})},
        "unknowns_guidance": (
            "Some compound parameters cannot be read off the given data and must "
            "be estimated so the model reproduces the observations. Decide which "
            "to fix at literature values and which to estimate, choose a "
            "distribution/permeability approach, and justify each choice."),
    }

    blanked, answer_edits, warns = _blank(data)
    answer_key = {
        "schema": "osp-answer-key/v1",
        "warning": "REFERENCE ANSWERS - grading only, never shown to the agent.",
        "compound": comp.get("Name"),
        "structural_choices": {
            "calculation_methods": comp.get("CalculationMethods") or [],
            "metabolizing_processes": [{"molecule": p.get("Molecule"),
                                        "type": p.get("InternalName")}
                                       for p in comp.get("Processes") or []
                                       if p.get("Molecule")]},
        "estimated_parameters": fitted,
    }
    return {"skip": False, "stem": stem, "input": agent_input,
            "blanked": blanked,
            # qualified answer (<Name>@<Molecule> where a param name is fitted on
            # more than one process) so per-enzyme values are not collapsed.
            "answer_edits": {"parameters": fitted_answer(comp)},
            "answer_key": answer_key, "warnings": warns, "info": info}


def _write(base_dir: str, stem: str, res: dict) -> None:
    for sub in ("json_input", "benchmark", "answer_key"):
        os.makedirs(os.path.join(base_dir, sub), exist_ok=True)
    with open(os.path.join(base_dir, "json_input", f"{stem}.input.json"), "w",
              encoding="utf-8") as fh:
        json.dump(res["input"], fh, indent=2, ensure_ascii=False)
    with open(os.path.join(base_dir, "benchmark", f"{stem}.blanked.json"), "w",
              encoding="utf-8") as fh:
        json.dump(res["blanked"], fh, indent=2, ensure_ascii=False)
    with open(os.path.join(base_dir, "answer_key", f"{stem}.answer_edits.json"), "w",
              encoding="utf-8") as fh:
        json.dump(res["answer_edits"], fh, indent=2, ensure_ascii=False)
    with open(os.path.join(base_dir, "answer_key", f"{stem}.answer_key.json"), "w",
              encoding="utf-8") as fh:
        json.dump(res["answer_key"], fh, indent=2, ensure_ascii=False)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("target", nargs="?", default=".",
                    help="a snapshot .json, or a directory to scan")
    ap.add_argument("--audit", action="store_true",
                    help="classify snapshots and flag exceptions; write nothing")
    ap.add_argument("--write", action="store_true",
                    help="actually write the files (default is a dry preview)")
    ap.add_argument("--selftest", action="store_true",
                    help="verify the generator reproduces the Alfentanil gold files "
                         "and leaks nothing")
    args = ap.parse_args()

    if args.selftest:
        _selftest()
        return

    if os.path.isfile(args.target):
        files = [args.target]
    else:
        files = sorted(glob.glob(os.path.join(args.target, "*/json/*.json")))

    if args.audit:
        buckets: dict[str, list[str]] = {}
        print(f"{'snapshot':40} {'kind':12} {'#cmp':4} {'#fit':4} {'#obs':5} exceptions")
        for f in files:
            with open(f, encoding="utf-8") as fh:
                info = classify(json.load(fh))
            stem = os.path.basename(f)
            buckets.setdefault(info["kind"], []).append(stem)
            exc = []
            if info["kind"] != "single":
                exc.append(info["kind"].upper())
            if info["missing_defaults"]:
                exc.append("no-default:" + ",".join(info["missing_defaults"]))
            if info["n_observed"] == 0:
                exc.append("NO-OBSERVED-DATA")
            print(f"{stem:40} {info['kind']:12} {info['n_compounds']:<4} "
                  f"{info['n_fitted']:<4} {info['n_observed']:<5} {'; '.join(exc)}")
        print("\nby kind:", {k: len(v) for k, v in sorted(buckets.items())})
        return

    for f in files:
        res = build_files(f)
        stem = res["stem"]
        if res["skip"]:
            print(f"SKIP {stem}: {res['reason']} (needs category-specific handling)")
            continue
        base = os.path.dirname(os.path.dirname(f))   # <Compound>/
        n_obs = res["info"]["n_observed"]
        n_fit = res["info"]["n_fitted"]
        warn = (" WARN no-default: " + ",".join(res["warnings"])) if res["warnings"] else ""
        print(f"{'WROTE' if args.write else 'DRY  '} {stem}: "
              f"{n_fit} fitted, {n_obs} observed datasets{warn}")
        if args.write:
            _write(base, stem, res)


def _selftest() -> None:
    """Regression: reproduce Alfentanil's gold fitted set and blank with no leak."""
    ref = "Alfentanil/json/Alfentanil-Model.json"
    if not os.path.exists(ref):
        print("selftest skipped: Alfentanil snapshot not found"); return
    res = build_files(ref)
    gold = json.load(open("Alfentanil/answer_key/Alfentanil-Model.answer_edits.json"))
    gold_p = gold.get("parameters", gold)
    gen_p = res["answer_edits"]["parameters"]
    assert sorted(gen_p) == sorted(gold_p), (sorted(gen_p), sorted(gold_p))
    assert all(abs(gen_p[k] - gold_p[k]) < 1e-9 for k in gen_p)
    # blanked values must differ from the fitted answer (no leak)
    def find(o, nm):
        if isinstance(o, dict):
            if o.get("Name") == nm and "Value" in o:
                return o["Value"]
            for v in o.values():
                r = find(v, nm)
                if r is not None:
                    return r
        elif isinstance(o, list):
            for v in o:
                r = find(v, nm)
                if r is not None:
                    return r
    comp = res["blanked"]["Compounds"][0]
    for k, v in gen_p.items():
        assert abs((find(comp, k) or 0) - v) > 1e-9, f"{k} not blanked (leak!)"
    # literature (given) must not include any fitted parameter
    lit = {l["parameter"] for l in res["input"]["given_data"]["literature_physicochemical"]}
    assert not (lit & set(gen_p)), lit & set(gen_p)
    print("selftest OK: Alfentanil reproduced (5 fitted, blanked no-leak, no input leak)")


if __name__ == "__main__":
    main()
