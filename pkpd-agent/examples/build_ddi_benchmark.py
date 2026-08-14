r"""Generate a DDI benchmark from a perpetrator+victim OSP snapshot.

The DDI analogue of build_benchmark.py. A victim-DDI model's ANSWER is not a set
of single-compound drug parameters but the perpetrator's INTERACTION parameters
(e.g. erythromycin's kinact / K_kinact_half on CYP3A4) that were identified. The
benchmark quantity is the interaction RATIO (victim AUC treatment/control).

This reads each snapshot's own metadata (ValueOrigin) to find which interaction
parameters were fitted, blanks them to no-leak defaults, and writes:

    <Compound>/json_input/<stem>.ddi_input.json     agent-facing DDI task
    <Compound>/benchmark/<stem>.ddi_blanked.json    agent start point
    <Compound>/answer_key/<stem>.ddi_answer_edits.json   fitted interaction params
    <Compound>/answer_key/<stem>.ddi_answer_key.json     + observed ratios

    python -m examples.build_ddi_benchmark ^
        ..\OSP-PBPK-Model-Library\Erythromycin\json\Erythromycin-Model.json --write
    python -m examples.build_ddi_benchmark --audit ..\OSP-PBPK-Model-Library
"""

from __future__ import annotations

import argparse
import copy
import glob
import json
import os

from pkpd_agent.engines import osp_ddi
from pkpd_agent.engines.osp_ddi import _auc

INTERACTION_NAMES = {"CompetitiveInhibition", "UncompetitiveInhibition",
                     "NoncompetitiveInhibition", "MixedInhibition",
                     "IrreversibleInhibition", "Induction"}

# no-leak naive priors for the interaction parameters (independent of the fit)
INTERACTION_DEFAULTS = {"kinact": 0.1, "K_kinact_half": 1.0, "Ki": 1.0,
                        "Ki_c": 1.0, "Ki_u": 1.0, "EC50": 1.0, "Emax": 1.0}


def fitted_interaction_params(snapshot: dict) -> list[dict]:
    """The interaction parameters that were IDENTIFIED (the DDI answer)."""
    specs = []
    for c in snapshot.get("Compounds") or []:
        for p in c.get("Processes") or []:
            if p.get("InternalName") not in INTERACTION_NAMES:
                continue
            fitted = {par["Name"]: par["Value"] for par in p.get("Parameters") or []
                      if (par.get("ValueOrigin") or {}).get("Source")
                      == "ParameterIdentification"
                      and isinstance(par.get("Value"), (int, float))}
            if fitted:
                specs.append({"perpetrator": c.get("Name"),
                              "internal_name": p.get("InternalName"),
                              "target": p.get("Molecule"),
                              "parameters": fitted})
    return specs


def _blank(snapshot: dict, specs: list) -> dict:
    snap = copy.deepcopy(snapshot)
    want = {(s["perpetrator"], s["internal_name"], s["target"]): set(s["parameters"])
            for s in specs}
    for c in snap.get("Compounds") or []:
        for p in c.get("Processes") or []:
            key = (c.get("Name"), p.get("InternalName"), p.get("Molecule"))
            names = want.get(key)
            if not names:
                continue
            for par in p.get("Parameters") or []:
                if par.get("Name") in names:
                    par["Value"] = INTERACTION_DEFAULTS.get(par["Name"], 1.0)
                    par["ValueOrigin"] = {"Source": "Unknown",
                                          "Description": "benchmark naive prior (blanked)"}
    return snap


_MASS = {"mg/l": 1.0, "µg/l": 1e-3, "ug/l": 1e-3, "ng/ml": 1e-3, "ng/l": 1e-6, "g/l": 1e3}


def _to_mg_L(vals, unit, mw):
    u = (unit or "").strip().lower()
    out = []
    for v in vals:
        if not isinstance(v, (int, float)):
            out.append(None); continue
        if u in _MASS:
            out.append(v * _MASS[u])
        elif "mol" in u and mw:
            f = mw / 1000.0 * (1e-3 if u.startswith("nmol") else 1.0)
            out.append(v * f)
        else:
            out.append(v)
    return out


def _observed_ratios_and_profiles(snapshot: dict, victim: str,
                                  perpetrators: list | None = None) -> dict:
    """Extract the victim's control/treatment plasma profiles (time, conc) and
    compute the OBSERVED interaction ratio AUCR = AUC_treatment / AUC_control per
    study+route. Among the VICTIM's profiles, a TREATMENT arm is one whose name
    mentions the perpetrator (victim co-administered); a CONTROL arm is the victim
    alone. The named 'AUCratio_*' datasets are a theoretical-AUC plotting artifact
    (not the ratio) and are ignored."""
    perp_l = [p.lower() for p in (perpetrators or [])]
    profiles = []
    for od in snapshot.get("ObservedData") or []:
        name = od.get("Name") or ""
        if "aucratio" in name.lower():
            continue
        ep = {e.get("Name"): e.get("Value") for e in od.get("ExtendedProperties") or []}
        if ep.get("Molecule") != victim:
            continue
        base = od.get("BaseGrid") or {}
        col = (od.get("Columns") or [{}])[0]
        mw = (col.get("DataInfo") or {}).get("MolWeight")
        time_h = [v for v in (base.get("Values") or []) if isinstance(v, (int, float))]
        conc = _to_mg_L(list(col.get("Values") or []), col.get("Unit", ""), mw)
        n = min(len(time_h), len(conc))
        if n < 2:
            continue
        nm = name.lower()
        mentions_perp = any(p in nm for p in perp_l)
        explicit_control = ("control" in nm or "placebo" in nm or "alone" in nm)
        # among the victim's own profiles: with the perpetrator named -> treatment;
        # otherwise (or explicitly labelled control) -> the victim alone.
        is_control = explicit_control or not mentions_perp
        is_treat = mentions_perp and not explicit_control
        profiles.append({"dataset": name, "study": ep.get("Study Id"),
                         "route": ep.get("Route"),
                         "arm": "control" if is_control else ("treatment" if is_treat else "?"),
                         "time_h": [round(t, 5) for t in time_h[:n]],
                         "conc_mg_L": [None if c is None else round(c, 8) for c in conc[:n]]})

    # pair control<->treatment on (study, route) and compute observed AUCR
    def key(p):
        r = (p.get("route") or "").lower()
        r = "PO" if r in ("po", "oral") else ("IV" if r == "iv" else r)
        return ((p.get("study") or "").lower(), r)
    ctrl = {key(p): p for p in profiles if p["arm"] == "control"}
    ratios = []
    for p in profiles:
        if p["arm"] != "treatment":
            continue
        c = ctrl.get(key(p))
        if not c:
            continue
        ac = _auc(c["time_h"], [x for x in c["conc_mg_L"]])
        at = _auc(p["time_h"], [x for x in p["conc_mg_L"]])
        if ac > 0:
            ratios.append({"study": p.get("study"), "route": key(p)[1],
                           "treatment_dataset": p["dataset"], "control_dataset": c["dataset"],
                           "observed_aucr": round(at / ac, 3)})
    return {"observed_ratios": ratios, "victim_profiles": profiles}


def build(path: str) -> dict:
    with open(path, encoding="utf-8") as fh:
        snap = json.load(fh)
    stem = os.path.splitext(os.path.basename(path))[0]
    ddi = osp_ddi.analyze_ddi(snap)
    if not ddi:
        return {"skip": True, "stem": stem, "reason": "not a DDI snapshot"}
    if not ddi.get("victims") or not ddi.get("pairs"):
        return {"skip": True, "stem": stem,
                "reason": "no victim / no control-treatment pairs (auto-TDI only)"}

    specs = fitted_interaction_params(snap)
    if not specs:
        return {"skip": True, "stem": stem,
                "reason": "no fitted interaction parameters to recover"}

    victim = ddi["victims"][0]
    perp_names = [p["name"] for p in ddi["perpetrators"]]
    obs = _observed_ratios_and_profiles(snap, victim, perp_names)
    blanked = _blank(snap, specs)

    perp = ddi["perpetrators"][0]["name"]
    agent_input = {
        "schema": "osp-ddi-task/v1",
        "task_id": stem, "source_snapshot": os.path.basename(path),
        "objective": (
            f"Predict the drug-drug interaction: how {perp} (perpetrator) changes "
            f"the exposure of {victim} (victim). The victim and perpetrator PBPK "
            "models are given; the interaction MECHANISM is present but its "
            "parameters are not set. Choose/estimate the interaction parameters so "
            "the model reproduces the observed victim exposure with and without the "
            "perpetrator, and report the predicted interaction ratio "
            "(AUC_treatment / AUC_control)."),
        "ddi_structure": {
            "perpetrators": [{"name": p["name"],
                              "mechanisms": [{"internal_name": m["internal_name"],
                                              "target": m["target"]}
                                             for m in p["mechanisms"]]}
                             for p in ddi["perpetrators"]],
            "victim": victim,
            "control_treatment_pairs": ddi["pairs"],
        },
        "given_data": {
            "victim_observed_profiles": obs["victim_profiles"],
            "observed_interaction_ratios": obs["observed_ratios"],
        },
        "unknowns_guidance": (
            "The interaction parameters (e.g. kinact & K_kinact_half for "
            "mechanism-based inhibition, or Ki for competitive inhibition, or "
            "EC50 & Emax for induction) are NOT given. Estimate them from the "
            "observed control vs treatment victim data. Do not change the victim's "
            "own disposition parameters; the interaction is the unknown."),
        "how_scored": ("Predicted vs observed interaction ratio (AUCR = "
                       "AUC_treatment/AUC_control): GMFE on the ratio, and "
                       "fraction within 2-fold."),
    }
    answer_edits = {"interaction_parameters": specs}
    answer_key = {
        "schema": "osp-ddi-answer-key/v1",
        "warning": "REFERENCE ANSWERS - grading only, never shown to the agent.",
        "perpetrator": perp, "victim": victim,
        "fitted_interaction_parameters": specs,
        "observed_interaction_ratios": obs["observed_ratios"],
    }
    return {"skip": False, "stem": stem, "input": agent_input, "blanked": blanked,
            "answer_edits": answer_edits, "answer_key": answer_key,
            "n_pairs": len(ddi["pairs"]), "specs": specs, "victim": victim,
            "n_ratios": len(obs["observed_ratios"])}


def _write(base_dir: str, stem: str, res: dict) -> None:
    for sub in ("json_input", "benchmark", "answer_key"):
        os.makedirs(os.path.join(base_dir, sub), exist_ok=True)
    files = {
        os.path.join(base_dir, "json_input", f"{stem}.ddi_input.json"): res["input"],
        os.path.join(base_dir, "benchmark", f"{stem}.ddi_blanked.json"): res["blanked"],
        os.path.join(base_dir, "answer_key", f"{stem}.ddi_answer_edits.json"): res["answer_edits"],
        os.path.join(base_dir, "answer_key", f"{stem}.ddi_answer_key.json"): res["answer_key"],
    }
    for path, obj in files.items():
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(obj, fh, indent=2, ensure_ascii=False)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("target", help="a DDI snapshot .json, or a directory to scan")
    ap.add_argument("--write", action="store_true")
    ap.add_argument("--audit", action="store_true")
    args = ap.parse_args()

    files = ([args.target] if os.path.isfile(args.target)
             else sorted(glob.glob(os.path.join(args.target, "*/json/*.json"))))
    for f in files:
        res = build(f)
        stem = res["stem"]
        if res["skip"]:
            if args.audit or os.path.isfile(args.target):
                print(f"SKIP {stem}: {res['reason']}")
            continue
        print(f"{'WROTE' if (args.write and not args.audit) else 'DDI '} {stem}: "
              f"perpetrator->{res['victim']}, {len(res['specs'])} interaction "
              f"mechanism(s) blanked, {res['n_pairs']} pairs, "
              f"{res['n_ratios']} observed ratio(s)")
        for s in res["specs"]:
            print(f"      {s['perpetrator']}: {s['internal_name']} on {s['target']} "
                  f"-> recover {list(s['parameters'])}")
        if args.write and not args.audit:
            _write(os.path.dirname(os.path.dirname(f)), stem, res)


if __name__ == "__main__":
    main()
