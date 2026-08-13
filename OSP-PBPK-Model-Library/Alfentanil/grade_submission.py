"""Grade an agent's PBPK submission against a task input (and optional answer key).

Two layers, matching the rubric in the agent input:

  1. NUMERICAL (deterministic, stdlib only) --------------------------------
     * data_fit             : pairs predicted vs observed plasma concentrations
                              (interpolating predictions onto the observed time
                              grid), then GMFE and % within 2-fold, reported
                              overall / by route / by study.
     * parameter_plausibility: rule-based physiological bound checks on every
                              parameter the agent reported.
     * output_plausibility  : rule-based checks on the predicted profiles
                              (non-negativity, IV decline / oral rise-then-fall,
                              dose-ordering of Cmax & AUC, terminal-half-life
                              spread).

  2. AGENTIC (optional) ----------------------------------------------------
     A Claude judge reads the numerical scorecard plus the agent's structural
     choices and rationales and renders the PHYSICAL-REASONING verdict the
     numbers can't: is the model mechanistically sound, are the flags real
     problems or acceptable, and what should the agent change. Runs only if
     `anthropic` is installed and ANTHROPIC_API_KEY is set (or --reason given);
     otherwise the numerical scorecard is produced alone.

Usage:
    python grade_submission.py --input json_input/Alfentanil-Model.input.json \\
                               --submission demo_submission.json
    # options:
    #   --key answer_key/Alfentanil-Model.answer_key.json   (auxiliary closeness)
    #   --reason            force the agentic layer (errors if it can't run)
    #   --no-reason         skip the agentic layer even if a key is present
    #   --model claude-opus-5   --effort medium
    #   --out scorecard.json
    #   --selftest          run the GMFE math self-check and exit
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
from typing import Any


# --------------------------------------------------------------------------- #
# small helpers
# --------------------------------------------------------------------------- #

def _finite(x: Any) -> float | None:
    try:
        f = float(x)
    except (TypeError, ValueError):
        return None
    return None if (f != f or f in (float("inf"), float("-inf"))) else f


def _interp(xp: list[float], fp: list[float], x: float) -> float | None:
    """Linear interpolation of (xp, fp) at x; None if out of range or degenerate."""
    pts = sorted((a, b) for a, b in zip(xp, fp)
                 if _finite(a) is not None and _finite(b) is not None)
    if len(pts) < 2 or x < pts[0][0] or x > pts[-1][0]:
        # allow exact endpoint match with a single point
        for a, b in pts:
            if abs(a - x) < 1e-9:
                return b
        return None
    for i in range(1, len(pts)):
        x0, y0 = pts[i - 1]
        x1, y1 = pts[i]
        if x0 <= x <= x1:
            if x1 == x0:
                return y0
            t = (x - x0) / (x1 - x0)
            return y0 + t * (y1 - y0)
    return None


def _parse_dose(dose: Any) -> tuple[float | None, str]:
    if not isinstance(dose, str):
        return None, ""
    m = re.match(r"\s*([0-9.eE+-]+)\s*(\S+)?", dose)
    if not m:
        return None, ""
    return _finite(m.group(1)), (m.group(2) or "").strip()


def _auc(t: list[float], c: list[float]) -> float | None:
    pts = sorted((a, b) for a, b in zip(t, c)
                 if _finite(a) is not None and _finite(b) is not None)
    if len(pts) < 2:
        return None
    return sum((pts[i][0] - pts[i - 1][0]) * (pts[i][1] + pts[i - 1][1]) / 2
               for i in range(1, len(pts)))


def _t_half(t: list[float], c: list[float]) -> float | None:
    pts = [(a, b) for a, b in zip(t, c)
           if _finite(a) is not None and _finite(b) is not None and b > 0][-3:]
    if len(pts) < 2 or pts[0][1] <= pts[-1][1] or pts[-1][0] <= pts[0][0]:
        return None
    k = (math.log(pts[0][1]) - math.log(pts[-1][1])) / (pts[-1][0] - pts[0][0])
    return math.log(2) / k if k > 0 else None


# --------------------------------------------------------------------------- #
# load
# --------------------------------------------------------------------------- #

def load_task(input_path: str) -> dict[str, Any]:
    with open(input_path, encoding="utf-8") as fh:
        d = json.load(fh)
    obs = {}
    for o in d.get("given_data", {}).get("clinical_observed_data", []):
        obs[o["dataset"]] = {
            "time_h": o.get("time_h", []),
            "conc_mg_L": o.get("conc_mg_L", []),
            "route": o.get("route"),
            "dose": o.get("dose"),
            "study": o.get("study"),
        }
    return {
        "compound": d.get("compound"),
        "objective": d.get("objective"),
        "rubric": d.get("evaluation_rubric"),
        "observed": obs,
    }


def load_submission(path: str) -> dict[str, Any]:
    with open(path, encoding="utf-8") as fh:
        d = json.load(fh)
    return d.get("submission", d)          # accept either wrapped or bare


# --------------------------------------------------------------------------- #
# 1a. data fit
# --------------------------------------------------------------------------- #

def grade_data_fit(observed: dict, submission: dict) -> dict[str, Any]:
    preds = {p.get("dataset"): p for p in submission.get("predicted_profiles") or []}
    pairs_all: list[tuple[float, float, str, str]] = []   # (obs, pred, route, study)
    per_dataset = []
    missing = []

    for name, ob in observed.items():
        pr = preds.get(name)
        if not pr:
            missing.append(name)
            continue
        pt, pc = pr.get("time_h", []), pr.get("pred_conc_mg_L", [])
        n = 0
        ds_pairs = []
        for t, o in zip(ob["time_h"], ob["conc_mg_L"]):
            of = _finite(o)
            tf = _finite(t)
            if of is None or of <= 0 or tf is None:
                continue
            # if pred shares the exact grid, pair by value; else interpolate
            p = None
            if len(pt) == len(ob["time_h"]):
                idx = ob["time_h"].index(t) if t in ob["time_h"] else None
                if idx is not None and idx < len(pc):
                    p = _finite(pc[idx])
            if p is None:
                p = _interp(pt, pc, tf)
            if p is None or p <= 0:
                continue
            ds_pairs.append((of, p))
            pairs_all.append((of, p, ob["route"] or "NA", ob["study"] or "NA"))
            n += 1
        per_dataset.append({"dataset": name, "n_paired": n,
                            **_fit_metrics([(a, b) for a, b in ds_pairs])})

    def agg(sel):
        return _fit_metrics([(o, p) for o, p, r, s in pairs_all if sel(r, s)])

    by_route = {}
    for r in sorted({r for _, _, r, _ in pairs_all}):
        by_route[r] = agg(lambda rr, ss, want=r: rr == want)
    by_study = {}
    for s in sorted({s for _, _, _, s in pairs_all}):
        by_study[s] = agg(lambda rr, ss, want=s: ss == want)

    return {
        "overall": _fit_metrics([(o, p) for o, p, r, s in pairs_all]),
        "by_route": by_route,
        "by_study": by_study,
        "per_dataset": per_dataset,
        "missing_predictions": missing,
        "n_datasets_scored": len(observed) - len(missing),
    }


def _fit_metrics(pairs: list[tuple[float, float]]) -> dict[str, Any]:
    fe = [p / o for o, p in pairs if o > 0 and p > 0]
    if not fe:
        return {"n": 0, "gmfe": None, "pct_within_2fold": None}
    log_abs = [abs(math.log(f)) for f in fe]
    gmfe = math.exp(sum(log_abs) / len(log_abs))
    within = sum(1 for f in fe if 0.5 <= f <= 2.0) / len(fe)
    return {"n": len(fe), "gmfe": round(gmfe, 3),
            "pct_within_2fold": round(100 * within, 1)}


# --------------------------------------------------------------------------- #
# 1b. parameter plausibility (rule-based)
# --------------------------------------------------------------------------- #

def grade_parameters(submission: dict) -> dict[str, Any]:
    checks = []

    def add(param, value, unit, status, message):
        checks.append({"parameter": param, "value": value, "unit": unit,
                       "status": status, "message": message})

    for p in submission.get("parameters") or []:
        name = (p.get("parameter") or "").lower()
        v = _finite(p.get("value"))
        unit = (p.get("unit") or "")
        if v is None:
            add(p.get("parameter"), p.get("value"), unit, "warn",
                "non-numeric value; cannot check")
            continue
        u = unit.lower()

        if "unbound" in name or name in ("fu", "fup") or "fraction unbound" in name:
            add(p.get("parameter"), v, unit,
                "ok" if 0 < v <= 1 else "flag",
                "fraction unbound must be in (0, 1]")
        elif "lipophilic" in name or name in ("logp", "logd", "logma"):
            add(p.get("parameter"), v, unit,
                "ok" if -2 <= v <= 7 else "flag",
                "lipophilicity (logP/logD) expected within about [-2, 7]")
        elif "clearance" in name or re.search(r"\bcl\b", name):
            if v <= 0:
                add(p.get("parameter"), v, unit, "flag", "clearance must be > 0")
            elif "l/min" in u and v > 1.5:
                add(p.get("parameter"), v, unit, "warn",
                    "clearance > ~1.5 L/min exceeds adult hepatic blood flow — check")
            elif ("l/h" in u) and v > 90:
                add(p.get("parameter"), v, unit, "warn",
                    "clearance > ~90 L/h exceeds adult hepatic blood flow — check")
            else:
                add(p.get("parameter"), v, unit, "ok", "clearance positive & bounded")
        elif "permeab" in name:
            add(p.get("parameter"), v, unit,
                "ok" if v > 0 else "flag", "permeability must be > 0")
        elif "volume" in name or re.search(r"\bvd?\b", name) or "distribution" in name:
            if "l/kg" in u:
                add(p.get("parameter"), v, unit,
                    "ok" if 0.05 <= v <= 50 else "flag",
                    "volume of distribution expected ~0.05-50 L/kg")
            else:
                add(p.get("parameter"), v, unit,
                    "ok" if v > 0 else "flag", "volume must be > 0")
        else:
            add(p.get("parameter"), v, unit,
                "ok" if v > 0 else "warn",
                "no specific bound; checked finite/positive")

    n_flag = sum(1 for c in checks if c["status"] == "flag")
    n_warn = sum(1 for c in checks if c["status"] == "warn")
    return {"n_parameters": len(checks), "n_flag": n_flag, "n_warn": n_warn,
            "checks": checks}


# --------------------------------------------------------------------------- #
# 1c. output plausibility (rule-based)
# --------------------------------------------------------------------------- #

def grade_outputs(observed: dict, submission: dict) -> dict[str, Any]:
    findings = []
    profiles = submission.get("predicted_profiles") or []
    route_of = {n: (observed.get(n, {}) or {}).get("route") for n in observed}
    dose_of = {n: _parse_dose((observed.get(n, {}) or {}).get("dose")) for n in observed}

    summ = {}   # dataset -> {cmax, tmax, auc, thalf, route}
    for pr in profiles:
        name = pr.get("dataset")
        t = pr.get("time_h", [])
        c = pr.get("pred_conc_mg_L", [])
        cc = [_finite(x) for x in c]
        if any(x is not None and x < 0 for x in cc):
            findings.append({"dataset": name, "status": "flag",
                             "message": "negative predicted concentration"})
        pos = [(a, b) for a, b in zip(t, cc) if _finite(a) is not None and b is not None]
        if len(pos) < 2:
            continue
        vals = [b for _, b in pos]
        imax = max(range(len(vals)), key=lambda i: vals[i])
        cmax = vals[imax]
        route = (route_of.get(name) or "").upper()

        if "IV" in route:
            # IV: a bolus peaks at t0; an infusion peaks at end-of-infusion. Both
            # are valid, so we don't warn on a late peak. We only flag a genuine
            # rebound: a secondary rise well AFTER the peak (non-monotone tail).
            rebound = sum(1 for i in range(imax + 2, len(vals))
                          if vals[i] > vals[i - 1] * 1.10)
            if rebound >= 2:
                findings.append({"dataset": name, "status": "warn",
                                 "message": "IV profile rises again after the peak "
                                            "(non-monotone elimination)"})
        else:
            # oral: expect a rise to Cmax then a decline
            if imax == 0:
                findings.append({"dataset": name, "status": "warn",
                                 "message": "oral profile peaks at first point "
                                            "(no absorption phase captured)"})
            if imax == len(vals) - 1:
                findings.append({"dataset": name, "status": "flag",
                                 "message": "oral profile still rising at last point "
                                            "(no elimination captured)"})

        summ[name] = {"route": route, "cmax": cmax,
                      "auc": _auc([a for a, _ in pos], vals),
                      "thalf": _t_half([a for a, _ in pos], vals)}

    # dose ordering of Cmax & AUC within each (route, dose-unit) group
    dose_checks = []
    groups: dict[tuple, list] = {}
    for name, s in summ.items():
        dv, du = dose_of.get(name, (None, ""))
        if dv is None:
            continue
        groups.setdefault((s["route"], du), []).append((dv, name, s))
    for (route, du), items in groups.items():
        if len({dv for dv, _, _ in items}) < 2:
            continue
        items.sort(key=lambda x: x[0])
        auc_ok = all(items[i][2]["auc"] is None or items[i - 1][2]["auc"] is None
                     or items[i][2]["auc"] >= items[i - 1][2]["auc"] * 0.9
                     for i in range(1, len(items)))
        cmax_ok = all(items[i][2]["cmax"] >= items[i - 1][2]["cmax"] * 0.9
                      for i in range(1, len(items)))
        dose_checks.append({"route": route, "dose_unit": du,
                            "doses": [dv for dv, _, _ in items],
                            "cmax_increases_with_dose": cmax_ok,
                            "auc_increases_with_dose": auc_ok})

    # terminal half-life spread per route
    thalf_checks = []
    by_route: dict[str, list] = {}
    for name, s in summ.items():
        if s["thalf"]:
            by_route.setdefault(s["route"], []).append(s["thalf"])
    for route, hs in by_route.items():
        if len(hs) >= 2:
            spread = max(hs) / min(hs) if min(hs) > 0 else None
            thalf_checks.append({
                "route": route, "n": len(hs),
                "t_half_h_min": round(min(hs), 2), "t_half_h_max": round(max(hs), 2),
                "spread_ratio": round(spread, 2) if spread else None,
                "status": "ok" if spread and spread <= 3 else "warn",
            })

    n_flag = sum(1 for f in findings if f["status"] == "flag")
    n_warn = sum(1 for f in findings if f["status"] == "warn")
    return {"n_flag": n_flag, "n_warn": n_warn, "profile_findings": findings,
            "dose_ordering": dose_checks, "terminal_half_life": thalf_checks}


# --------------------------------------------------------------------------- #
# auxiliary: closeness to the reference answer key (not a primary grade)
# --------------------------------------------------------------------------- #

def compare_to_key(submission: dict, key_path: str) -> dict[str, Any]:
    with open(key_path, encoding="utf-8") as fh:
        key = json.load(fh)
    ref = {(_norm(p["parameter"])): p for p in key.get("estimated_parameters") or []}
    rows = []
    for p in submission.get("parameters") or []:
        rp = ref.get(_norm(p.get("parameter") or ""))
        v = _finite(p.get("value"))
        if rp and v is not None and _finite(rp.get("value")):
            rv = float(rp["value"])
            fold = (v / rv) if rv else None
            rows.append({"parameter": p.get("parameter"),
                         "submitted": v, "reference": rv,
                         "fold_vs_reference": round(fold, 2) if fold else None})
    return {"note": "auxiliary only — many valid models differ from the reference",
            "matched_parameters": rows,
            "reference_methods": key.get("structural_choices", {}).get("calculation_methods")}


def _norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]", "", (s or "").lower())


# --------------------------------------------------------------------------- #
# 2. agentic physical-reasoning layer
# --------------------------------------------------------------------------- #

JUDGE_TOOL = {
    "name": "submit_judgment",
    "description": "Record the physical-reasoning verdict on the submission.",
    "input_schema": {
        "type": "object",
        "properties": {
            "per_dimension": {
                "type": "object",
                "properties": {
                    "data_fit": {"type": "object", "properties": {
                        "assessment": {"type": "string"},
                        "score_0_5": {"type": "integer"}}},
                    "parameter_plausibility": {"type": "object", "properties": {
                        "assessment": {"type": "string"},
                        "score_0_5": {"type": "integer"}}},
                    "output_plausibility": {"type": "object", "properties": {
                        "assessment": {"type": "string"},
                        "score_0_5": {"type": "integer"}}},
                },
            },
            "mechanistic_soundness": {"type": "string",
                "description": "Do the structure and choices make biological sense?"},
            "overall_verdict": {"type": "string", "enum": ["pass", "revise", "fail"]},
            "actionable_feedback": {"type": "array", "items": {"type": "string"},
                "description": "Concrete changes the agent could make to improve."},
            "summary": {"type": "string"},
        },
        "required": ["per_dimension", "overall_verdict", "summary"],
    },
}


def agentic_review(task: dict, submission: dict, numeric: dict,
                   model: str, effort: str) -> dict[str, Any]:
    import anthropic
    client = anthropic.Anthropic()

    payload = {
        "objective": task.get("objective"),
        "rubric": task.get("rubric"),
        "submission_structural_model": submission.get("structural_model"),
        "submission_parameters": submission.get("parameters"),
        "submission_self_assessment": submission.get("self_assessment"),
        "numerical_scorecard": {
            "data_fit": numeric["data_fit"]["overall"],
            "data_fit_by_route": numeric["data_fit"]["by_route"],
            "parameter_flags": [c for c in numeric["parameters"]["checks"]
                                if c["status"] != "ok"],
            "output_flags": numeric["outputs"]["profile_findings"],
            "dose_ordering": numeric["outputs"]["dose_ordering"],
            "terminal_half_life": numeric["outputs"]["terminal_half_life"],
        },
    }
    system = (
        "You are a senior pharmacometrician grading a PBPK model submission. "
        "The numbers (GMFE, % within 2-fold, rule-based flags) are given to you. "
        "Your job is the PHYSICAL REASONING the numbers can't do: judge whether "
        "the chosen structure and parameters are mechanistically sound, whether "
        "each flag is a genuine problem or acceptable given the drug's biology, "
        "and what the agent should change. Do NOT reward matching any specific "
        "reference model — reward a model that is both well-fitting and "
        "physically coherent. Call submit_judgment with your verdict.")
    msg = client.messages.create(
        model=model,
        max_tokens=4000,
        thinking={"type": "adaptive"},
        output_config={"effort": effort},
        tools=[JUDGE_TOOL],
        system=system,
        messages=[{"role": "user",
                   "content": "Grade this submission.\n\n"
                              + json.dumps(payload, ensure_ascii=False, indent=2)}],
    )
    for block in msg.content:
        if getattr(block, "type", None) == "tool_use" and block.name == "submit_judgment":
            return block.input
    # fallback: parse a JSON object from any text block
    for block in msg.content:
        if getattr(block, "type", None) == "text":
            m = re.search(r"\{.*\}", block.text, re.DOTALL)
            if m:
                try:
                    return json.loads(m.group(0))
                except json.JSONDecodeError:
                    pass
    return {"error": "judge did not return a structured verdict"}


# --------------------------------------------------------------------------- #
# orchestration
# --------------------------------------------------------------------------- #

def grade(input_path: str, submission_path: str, key_path: str | None,
          reason: bool | None, model: str, effort: str) -> dict[str, Any]:
    task = load_task(input_path)
    submission = load_submission(submission_path)

    numeric = {
        "data_fit": grade_data_fit(task["observed"], submission),
        "parameters": grade_parameters(submission),
        "outputs": grade_outputs(task["observed"], submission),
    }
    scorecard: dict[str, Any] = {
        "schema": "osp-scorecard/v1",
        "compound": task["compound"],
        "input": os.path.basename(input_path),
        "submission": os.path.basename(submission_path),
        "numerical": numeric,
    }
    if key_path:
        scorecard["reference_closeness"] = compare_to_key(submission, key_path)

    # decide whether to run the agentic layer
    want = reason
    if want is None:
        want = bool(os.environ.get("ANTHROPIC_API_KEY"))
    if want:
        try:
            scorecard["agentic_review"] = agentic_review(
                task, submission, numeric, model, effort)
        except Exception as exc:                       # noqa: BLE001
            scorecard["agentic_review"] = {"skipped": f"{type(exc).__name__}: {exc}"}
    else:
        scorecard["agentic_review"] = {
            "skipped": "no ANTHROPIC_API_KEY (numerical layer only); "
                       "use --reason to force"}
    return scorecard


def _print_summary(sc: dict[str, Any]) -> None:
    df = sc["numerical"]["data_fit"]["overall"]
    pr = sc["numerical"]["parameters"]
    op = sc["numerical"]["outputs"]
    print("=" * 66)
    print(f"scorecard: {sc['compound']}  ({sc['submission']})")
    print("=" * 66)
    print(f"DATA FIT   overall GMFE={df.get('gmfe')}  "
          f"within-2fold={df.get('pct_within_2fold')}%  (n={df.get('n')})")
    for r, m in sc["numerical"]["data_fit"]["by_route"].items():
        print(f"           {r:4} GMFE={m.get('gmfe')}  within-2fold={m.get('pct_within_2fold')}%")
    miss = sc["numerical"]["data_fit"]["missing_predictions"]
    if miss:
        print(f"           MISSING predictions for {len(miss)} dataset(s)")
    print(f"PARAMS     {pr['n_flag']} flag / {pr['n_warn']} warn / {pr['n_parameters']} total")
    for c in pr["checks"]:
        if c["status"] != "ok":
            print(f"           [{c['status']}] {c['parameter']}={c['value']} — {c['message']}")
    print(f"OUTPUTS    {op['n_flag']} flag / {op['n_warn']} warn")
    for f in op["profile_findings"]:
        print(f"           [{f['status']}] {f['dataset']}: {f['message']}")
    ar = sc.get("agentic_review", {})
    if "skipped" in ar:
        print(f"REASONING  skipped: {ar['skipped']}")
    elif "overall_verdict" in ar:
        print(f"REASONING  verdict={ar['overall_verdict'].upper()} — {ar.get('summary','')}")
        for fb in ar.get("actionable_feedback", []):
            print(f"           -> {fb}")
    print("=" * 66)


def _selftest() -> None:
    obs = {"D": {"time_h": [1, 2, 3], "conc_mg_L": [10, 5, 2.5],
                 "route": "IV", "dose": "1 mg", "study": "S"}}
    perfect = {"predicted_profiles": [{"dataset": "D", "time_h": [1, 2, 3],
                                       "pred_conc_mg_L": [10, 5, 2.5]}]}
    off = {"predicted_profiles": [{"dataset": "D", "time_h": [1, 2, 3],
                                   "pred_conc_mg_L": [20, 10, 5]}]}
    g1 = grade_data_fit(obs, perfect)["overall"]
    g2 = grade_data_fit(obs, off)["overall"]
    assert abs(g1["gmfe"] - 1.0) < 1e-6, g1
    assert abs(g2["gmfe"] - 2.0) < 1e-6, g2
    assert g1["pct_within_2fold"] == 100.0 and g2["pct_within_2fold"] == 100.0
    print("selftest OK: GMFE(pred=obs)=1.0, GMFE(pred=2*obs)=2.0")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", help="task input JSON (json_input/*.input.json)")
    ap.add_argument("--submission", help="the agent's submission JSON")
    ap.add_argument("--key", default=None, help="answer key (auxiliary closeness)")
    ap.add_argument("--reason", action="store_true", help="force the agentic layer")
    ap.add_argument("--no-reason", action="store_true", help="skip the agentic layer")
    ap.add_argument("--model", default=os.environ.get("ANTHROPIC_MODEL", "claude-opus-5"))
    ap.add_argument("--effort", default="medium")
    ap.add_argument("--out", default=None, help="write full scorecard JSON here")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    if args.selftest:
        _selftest()
        return
    if not args.input or not args.submission:
        ap.error("--input and --submission are required (or use --selftest)")

    reason = True if args.reason else (False if args.no_reason else None)
    sc = grade(args.input, args.submission, args.key, reason, args.model, args.effort)
    _print_summary(sc)
    out = args.out or os.path.join(
        "scorecards",
        os.path.basename(args.submission).replace(".json", "") + ".scorecard.json")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(sc, fh, indent=2, ensure_ascii=False)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
