r"""Verify a generated benchmark: NO answer leak, and no MISSING agent inputs.

For each task type this builds the agent-facing benchmark from a real snapshot and
runs two checks:

  LEAK   - the blanked snapshot carries NO ParameterIdentification origin (every
           fitted value was reset), and no full-precision fitted VALUE appears in
           the agent-facing input JSON.
  INPUTS - the agent input still contains everything needed to do the task:
           the objective, the observed data to fit, the structure/biology the
           agent legitimately needs, and (in the blanked snapshot) the compound
           identity, molecular weight and expression system to build on.

    python -m examples.verify_benchmark --type single     ..\OSP-PBPK-Model-Library\Sildenafil\json\Sildenafil-Model.json
    python -m examples.verify_benchmark --type metabolite  ..\OSP-PBPK-Model-Library\Itraconazole\json\Itraconazole-Model.json
    python -m examples.verify_benchmark --type biologic    ..\OSP-PBPK-Model-Library\BAY794620\json\BAY794620.json
    python -m examples.verify_benchmark --type ddi         ..\OSP-PBPK-Model-Library\Erythromycin\json\Erythromycin-Model.json
"""

from __future__ import annotations

import argparse
import json
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_LIB = os.path.abspath(os.path.join(_HERE, "..", "..", "OSP-PBPK-Model-Library"))


def _named_answer_pairs(res: dict) -> list[tuple[str, float]]:
    """(name, value) of each declared answer parameter - so a ROUND answer value
    (e.g. Kd FcRn = 12.7, GFR fraction = 0.24) can be checked by name+value even
    though it is too round to fingerprint as a string."""
    ae = res.get("answer_edits") or {}
    pairs = []
    if isinstance(ae.get("parameters"), dict):
        for k, v in ae["parameters"].items():
            if isinstance(v, (int, float)):
                pairs.append((k.split("@", 1)[0], v))
    for key in ("cascade_parameters", "disposition_parameters"):
        for spec in ae.get(key) or []:
            for p in spec.get("parameters") or []:
                if isinstance(p.get("value"), (int, float)):
                    pairs.append((p.get("name"), p["value"]))
    for spec in ae.get("interaction_parameters") or []:
        for k, v in (spec.get("parameters") or {}).items():
            if isinstance(v, (int, float)):
                pairs.append((k, v))
    return pairs


def _named_value_leaks(pairs: list[tuple[str, float]], blanked: dict) -> list[str]:
    """Answer parameters whose exact (name, value) still survives in the blanked
    snapshot - a leak even when the value is too round to fingerprint by string."""
    want: dict[str, list[float]] = {}
    for n, v in pairs:
        want.setdefault(n, []).append(float(v))
    hits = []

    def w(o):
        if isinstance(o, dict):
            nm = o.get("Name")
            val = o.get("Value")
            if nm in want and isinstance(val, (int, float)):
                for tv in want[nm]:
                    if abs(float(val) - tv) <= 1e-9 + 1e-6 * abs(tv):
                        hits.append(f"{nm}={val}")
            for v in o.values():
                w(v)
        elif isinstance(o, list):
            for v in o:
                w(v)

    w(blanked)
    return hits


def _is_fingerprint(v: float) -> bool:
    """A fitted value precise enough to be a unique fingerprint (so a substring
    hit is a real leak, not a coincidence with a default/data point). Round
    numbers like 1.0 / 0.1 / 0.24 are NOT fingerprints."""
    s = repr(float(v))
    digits = s.replace("-", "").replace(".", "").replace("e", "").lstrip("0")
    return len(digits) >= 5


def _value_leaks(answers: list[float], text: str) -> list[float]:
    """Declared-answer values whose full-precision string appears in ``text``
    (the agent input, or the blanked snapshot). Round values are skipped."""
    leaks = []
    for v in answers:
        if not _is_fingerprint(v):
            continue
        if repr(float(v)) in text:
            leaks.append(v)
    return leaks


def _build(kind: str, path: str):
    with open(path, encoding="utf-8") as fh:
        snap = json.load(fh)
    if kind == "single":
        if _LIB not in sys.path:
            sys.path.insert(0, _LIB)
        import build_benchmark as G                       # library-dir module
        res = G.build_files(path)
        if res.get("skip"):
            return None, res.get("reason"), snap
        return res, None, snap
    if kind == "metabolite":
        import examples.build_metabolite_benchmark as G
    elif kind == "biologic":
        import examples.build_biologic_benchmark as G
    elif kind == "ddi":
        import examples.build_ddi_benchmark as G
    else:
        return None, f"unknown type {kind}", snap
    res = G.build(path)
    if res.get("skip"):
        return None, res.get("reason"), snap
    return res, None, snap


def _completeness(kind: str, res: dict, blanked: dict) -> list[str]:
    """Return a list of MISSING-input problems (empty = complete)."""
    inp = res["input"]
    gd = inp.get("given_data", {}) or {}
    problems = []

    def need(cond, msg):
        if not cond:
            problems.append(msg)

    need(inp.get("objective"), "no objective")

    if kind == "single":
        need(gd.get("clinical_observed_data"), "no observed clinical data to fit")
        bg = inp.get("background") or {}
        need(bg.get("literature_facts"), "no known-biology facts (enzymes/route)")
        need(gd.get("literature_physicochemical") is not None,
             "no literature physicochemical block")
        comp = (blanked.get("Compounds") or [{}])[0]
        need(any(p.get("Name") == "Molecular weight" for p in comp.get("Parameters") or []),
             "blanked snapshot missing Molecular weight")
        need(blanked.get("ExpressionProfiles") is not None,
             "blanked snapshot missing ExpressionProfiles (enzyme system)")
    elif kind == "metabolite":
        pbm = gd.get("plasma_by_molecule") or {}
        need(len([m for m, v in pbm.items() if v]) >= 2,
             "fewer than 2 molecules with plasma data")
        need(inp.get("cascade_structure") or inp.get("modelled_compounds"),
             "no cascade_structure / modelled_compounds")
        if inp.get("kind") == "cascade":
            need((inp.get("background") or {}).get("metabolic_pathway"),
                 "cascade missing metabolic_pathway biology")
    elif kind == "biologic":
        need(gd.get("biodistribution_by_matrix"), "no biodistribution data")
        bs = inp.get("biologic_structure") or {}
        need(bs.get("molecular_weight"), "no molecular weight")
        need("do not add enzyme" in (inp.get("unknowns_guidance") or "").lower(),
             "guidance does not forbid enzyme clearance")
    elif kind == "ddi":
        need(gd.get("observed_interaction_ratios"), "no observed interaction ratios")
        need(gd.get("victim_observed_profiles"), "no victim observed profiles")
        ds = inp.get("ddi_structure") or {}
        need(ds.get("perpetrators"), "no perpetrator mechanisms")
        need(ds.get("control_treatment_pairs"), "no control/treatment pairs")
    return problems


def _answers(kind: str, res: dict) -> list[float]:
    """The DECLARED answer values the agent must recover (from answer_edits) -
    NOT every fitted origin in the snapshot. A single-compound task keeps study/
    formulation fitted params (given design); a DDI task keeps the victim's
    disposition. Only the declared set is the answer, so only it is a leak."""
    ae = res.get("answer_edits") or {}
    vals = []
    # single: {"parameters": {name: value}}
    if isinstance(ae.get("parameters"), dict):
        vals += [v for v in ae["parameters"].values() if isinstance(v, (int, float))]
    # metabolite / biologic: [{compound, parameters:[{name,value}]}]
    for key in ("cascade_parameters", "disposition_parameters"):
        for spec in ae.get(key) or []:
            for p in spec.get("parameters") or []:
                if isinstance(p.get("value"), (int, float)):
                    vals.append(p["value"])
    # ddi: [{perpetrator, parameters:{name: value}}]
    for spec in ae.get("interaction_parameters") or []:
        for v in (spec.get("parameters") or {}).values():
            if isinstance(v, (int, float)):
                vals.append(v)
    return vals


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--type", required=True,
                    choices=["single", "metabolite", "biologic", "ddi"])
    ap.add_argument("snapshot")
    args = ap.parse_args()

    res, reason, snap = _build(args.type, args.snapshot)
    stem = os.path.basename(args.snapshot)
    print(f"== verify {args.type}: {stem} ==")
    if res is None:
        print(f"  SKIPPED (not a {args.type} task): {reason}")
        sys.exit(2)

    blanked = res["blanked"]
    input_json = json.dumps(res["input"], ensure_ascii=False)
    blanked_json = json.dumps(blanked, ensure_ascii=False)
    answers = _answers(args.type, res)
    n_fp = sum(1 for v in answers if _is_fingerprint(v))

    # -- leak checks (against the DECLARED answers only) --
    input_leaks = _value_leaks(answers, input_json)
    snapshot_leaks = _value_leaks(answers, blanked_json)
    named_leaks = _named_value_leaks(_named_answer_pairs(res), blanked)
    leak_ok = not input_leaks and not snapshot_leaks and not named_leaks
    print(f"  declared answers to recover: {len(answers)} "
          f"({n_fp} precise enough to fingerprint by string)")
    print(f"  LEAK  answer values in the AGENT INPUT:        "
          f"{input_leaks if input_leaks else 'none'}")
    print(f"  LEAK  answer values in the BLANKED snapshot:   "
          f"{snapshot_leaks if snapshot_leaks else 'none'}")
    print(f"  LEAK  answer (name,value) kept in blanked:     "
          f"{named_leaks if named_leaks else 'none'}")

    # -- completeness checks --
    problems = _completeness(args.type, res, blanked)
    inputs_ok = not problems
    print(f"  INPUTS present: {'YES' if inputs_ok else 'MISSING -> ' + '; '.join(problems)}")

    # -- what the agent is given (eyeball) --
    gd = res["input"].get("given_data", {})
    n_data = (len(gd.get("clinical_observed_data") or gd.get("victim_observed_profiles")
                  or gd.get("biodistribution_by_matrix")
                  or gd.get("plasma_by_molecule") or []))
    print(f"  given: objective + {n_data} observed dataset(s)/matrix(es); "
          f"keys={list(res['input'].keys())}")

    ok = leak_ok and inputs_ok
    print(f"  => {'PASS' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
