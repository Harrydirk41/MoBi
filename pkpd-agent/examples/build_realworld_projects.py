r"""Build the 'real-world' project tree: raw report + data, like a modeler gets.

The benchmark's input.json is a PRE-DIGESTED context lake - we already read the
PDF and extracted logP / fu / Km / enzyme into clean JSON. A real modeler is not
handed that: they get a report and data files and must extract the inputs
themselves. This script rebuilds every OSP compound into that raw form, in a NEW
tree (pbpk-realworld/), leaving OSP-PBPK-Model-Library/ untouched.

Each compound gets BOTH versions:

  pbpk-realworld/<Compound>/
    context/            # this compound AS a reference for others (full disclosure)
      report/<C>.md     #   the complete evaluation report
      data/*.csv        #   every clinical study as raw time,conc
    test/               # this compound AS the target (what the agent is handed)
      report/<C>_inputs.md   the report with the ANSWER redacted:
                             keeps intro, physchem, clinical description, and the
                             literature in-vitro kinetics; drops the modeling
                             strategy, structure/method choices, and all results.
      data/*.csv        #   the BUILDING studies only (held-out studies withheld)
      README.md         #   "what you are handed" note
    answer/             # sealed - judge only
      held_out.txt      #   the verification studies used for grading
      reference.json    #   the reference's fitted parameters + chosen methods

Leave-one-out is enforced by the harness at run time (when testing X, read every
OTHER compound's context/, never X's own).

    python -m examples.build_realworld_projects            # dry preview
    python -m examples.build_realworld_projects --write
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import re
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_PKG = os.path.dirname(_HERE)
_LIB = os.path.abspath(os.path.join(_PKG, "..", "OSP-PBPK-Model-Library"))
_OUT = os.path.abspath(os.path.join(_PKG, "..", "pbpk-realworld"))

_FIT_ORIGINS = ("parameteridentification", "optimized", "fitted")


def _slug(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s).strip("_")[:80]


# --------------------------------------------------------------------------- #
# 1. clinical data -> CSV
# --------------------------------------------------------------------------- #
def _observed_csvs(model: dict) -> list[dict]:
    """One dict per ObservedData block: {name, study, route, dose, csv}."""
    out = []
    for od in model.get("ObservedData") or []:
        name = od.get("Name") or "dataset"
        ep = {e.get("Name"): e.get("Value") for e in od.get("ExtendedProperties") or []}
        base = od.get("BaseGrid") or {}
        col = (od.get("Columns") or [{}])[0]
        t = base.get("Values") or []
        c = col.get("Values") or []
        n = min(len(t), len(c))
        if n < 2:
            continue
        tu, cu = base.get("Unit", "h"), col.get("Unit", "mg/l")
        rows = [f"time[{tu}],concentration[{cu}]"]
        for i in range(n):
            rows.append(f"{t[i]},{c[i]}")
        # route/dose from the dataset name (e.g. "... - PO - 0.5 mg - Plasma ...")
        route = "IV" if re.search(r"\bIV\b", name) else ("PO" if re.search(r"\bPO\b", name) else "")
        dose = (re.search(r"([\d.]+ *mg)", name) or [None, ""])[1]
        out.append({"name": name, "study": ep.get("Study Id") or "study",
                    "reference": ep.get("Reference") or "", "route": route, "dose": dose,
                    "csv": "\n".join(rows)})
    return out


# --------------------------------------------------------------------------- #
# 2. report redaction  (keep inputs, drop the answer)
# --------------------------------------------------------------------------- #
def _split_sections(md: str) -> list[tuple[int, str, str]]:
    """Return [(level, heading, body)] for each '# '/'## ' section, in order."""
    lines = md.splitlines()
    secs, cur = [], None
    for ln in lines:
        m = re.match(r"^(#{1,2}) +(.*)$", ln)
        if m:
            if cur:
                secs.append(cur)
            cur = [len(m.group(1)), m.group(2).strip(), []]
        elif cur:
            cur[2].append(ln)
    if cur:
        secs.append(cur)
    return [(lvl, h, "\n".join(b).strip()) for lvl, h, b in secs]


def _scrub_method_rows(body: str) -> str:
    """Drop any line that names a chosen distribution/permeability method - some
    physchem tables list 'Partition coefficients | Rodgers & Rowland' as a row,
    which is a structural CHOICE (the answer), not an input."""
    bad = [m.lower() for m in _METHODS] + ["partition coefficient", "rodgers"]
    out = []
    for ln in body.splitlines():
        low = ln.lower()
        if ln.lstrip().startswith("|") and any(b in low for b in bad):
            continue
        # a stray inline method mention outside a table
        if any(m.lower() in low for m in _METHODS):
            continue
        out.append(ln)
    return "\n".join(out)


def _invitro_table(model: dict) -> str:
    """Literature in-vitro kinetics a modeler would have: enzyme identity + the
    MEASUREMENT-origin Km/Vmax (never the fitted kcat/CLint)."""
    rows = []
    comp = (model.get("Compounds") or [{}])[0]
    for p in comp.get("Processes") or []:
        mol = p.get("Molecule")
        typ = p.get("InternalName") or ""
        for par in p.get("Parameters") or []:
            nm = par.get("Name", "")
            origin = ((par.get("ValueOrigin") or {}).get("Source") or "").lower()
            if origin in _FIT_ORIGINS:
                continue                                # fitted -> that's the answer
            if any(k in nm for k in ("Km", "In vitro Vmax", "Vmax")):
                rows.append((mol or "-", typ, nm, par.get("Value"), par.get("Unit", "")))
    if not rows:
        return ""
    out = ["## Literature in-vitro kinetics (given)\n",
           "These are the measured in-vitro values from the literature - the fitted",
           "catalytic rates are NOT given; you must determine them.\n",
           "| molecule | process | parameter | value | unit |",
           "| --- | --- | --- | --- | --- |"]
    for mol, typ, nm, v, u in rows:
        out.append(f"| {mol} | {typ} | {nm} | {v} | {u} |")
    return "\n".join(out)


def redact_report(md: str, model: dict) -> str:
    """Keep intro + the DATA section + references; add the literature in-vitro
    table; drop modeling strategy, assumptions, results, and conclusion (the
    answer). Never emit a fitted value or a chosen calculation method."""
    secs = _split_sections(md)
    title = secs[0][1] if secs else "PBPK modeling task"
    keep = []
    for lvl, h, body in secs:
        hl = h.lower()
        if re.search(r"\bintroduction\b", hl):
            keep.append(("intro", "# Background", body))
        elif re.search(r"^2\.2\b|clinical data|in vitro / physico|in vitro / physicochem|physicochemical data", hl):
            # the DATA section and its subsections (physchem + clinical studies)
            keep.append(("data", f"## {h}", body))
    parts = [f"# {title}\n",
             "> You are handed this report and the clinical data files. Extract the",
             "> physicochemical and in-vitro inputs yourself, decide the model",
             "> structure, and fit the parameters. The reference model's chosen",
             "> methods and fitted values are withheld.\n"]
    for _kind, head, body in keep:
        parts.append(head + "\n\n" + _scrub_method_rows(body))
    iv = _invitro_table(model)
    if iv:
        parts.append(iv)
    return "\n\n".join(parts).strip() + "\n"


# --------------------------------------------------------------------------- #
# 3. answer (sealed)
# --------------------------------------------------------------------------- #
def _reference_answer(comp_dir: str, stem: str) -> dict:
    ak = os.path.join(_LIB, comp_dir, "answer_key", stem + ".answer_key.json")
    if os.path.isfile(ak):
        try:
            return json.load(open(ak, encoding="utf-8"))
        except Exception:
            pass
    return {}


# --------------------------------------------------------------------------- #
# leak guard
# --------------------------------------------------------------------------- #
_METHODS = ["Rodgers and Rowland", "Poulin and Theil", "Berezhkovskiy",
            "Charge dependent Schmitt"]


def _leak_check(test_report: str, answer: dict) -> list[str]:
    hits = []
    low = test_report.lower()
    for m in _METHODS:                                  # a chosen method must not appear
        if m.lower() in low:
            hits.append(f"method:{m}")
    if "gmfe" in low:
        hits.append("gmfe")
    for p in answer.get("estimated_parameters", []):    # no fitted value verbatim
        name = (p.get("parameter") or "")
        # fraction unbound is a MEASURED physchem input a modeler always has, even
        # when the reference also re-fits it - giving the measured fu is not a leak.
        if "unbound" in name.lower():
            continue
        v = p.get("value")
        if isinstance(v, (int, float)) and v not in (0.0, 1.0):
            if f"{v:.6g}" in test_report or repr(v) in test_report:
                hits.append(f"fitted:{name}={v}")
    return hits


# --------------------------------------------------------------------------- #
def build_one(comp_dir: str, write: bool) -> dict | None:
    reports = glob.glob(os.path.join(_LIB, comp_dir, "*_evaluation_report.md"))
    models = glob.glob(os.path.join(_LIB, comp_dir, "json", "*-Model.json")) \
        or glob.glob(os.path.join(_LIB, comp_dir, "json", "*.json"))
    models = [m for m in models if not m.endswith((".input.json", ".blanked.json"))]
    if not reports or not models:
        return None
    report_md = open(reports[0], encoding="utf-8").read()
    model = json.load(open(models[0], encoding="utf-8"))
    stem = os.path.basename(models[0])[:-5]

    csvs = _observed_csvs(model)
    if not csvs:
        return None
    from pkpd_agent.engines import osp_split
    split = osp_split.split_studies(model, [{"dataset": c["name"]} for c in csvs])
    build_names = set(split.get("building") or [c["name"] for c in csvs])
    held = split.get("held_out_studies") or []

    test_report = redact_report(report_md, model)
    answer = _reference_answer(comp_dir, stem)
    leaks = _leak_check(test_report, answer)

    res = {"compound": comp_dir, "n_studies": len(csvs), "n_build": len(build_names),
           "n_held": len(held), "leaks": leaks, "invitro": "Literature in-vitro" in test_report}

    if write and not leaks:
        root = os.path.join(_OUT, comp_dir)
        for sub in ("context/report", "context/data", "test/report", "test/data", "answer"):
            os.makedirs(os.path.join(root, sub), exist_ok=True)
        # context = full disclosure
        open(os.path.join(root, "context/report", f"{comp_dir}.md"), "w", encoding="utf-8").write(report_md)
        # test = redacted report
        open(os.path.join(root, "test/report", f"{comp_dir}_inputs.md"), "w", encoding="utf-8").write(test_report)
        # data
        man_c, man_t = [], []
        for c in csvs:
            fn = _slug(c["name"]) + ".csv"
            open(os.path.join(root, "context/data", fn), "w", encoding="utf-8").write(c["csv"])
            man_c.append({"file": fn, "study": c["study"], "route": c["route"], "dose": c["dose"]})
            if c["name"] in build_names:
                open(os.path.join(root, "test/data", fn), "w", encoding="utf-8").write(c["csv"])
                man_t.append({"file": fn, "study": c["study"], "route": c["route"], "dose": c["dose"]})
        json.dump(man_c, open(os.path.join(root, "context/data", "_manifest.json"), "w"), indent=2)
        json.dump(man_t, open(os.path.join(root, "test/data", "_manifest.json"), "w"), indent=2)
        # answer (sealed)
        open(os.path.join(root, "answer", "held_out.txt"), "w", encoding="utf-8").write("\n".join(held))
        json.dump(answer, open(os.path.join(root, "answer", "reference.json"), "w"), indent=2)
        open(os.path.join(root, "test", "README.md"), "w", encoding="utf-8").write(
            f"# {comp_dir} - modeling task\n\nYou are given `report/` (physchem, in-vitro, clinical "
            f"description) and `data/` ({len(man_t)} clinical studies as time,concentration CSVs). "
            "Extract the inputs, build the whole-body PBPK model, and fit the identifiable parameters. "
            "The reference model's methods and fitted values are withheld; your model is graded on how "
            "well it predicts held-out studies.\n")
    return res


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--write", action="store_true", help="write pbpk-realworld/ (default: preview)")
    args = ap.parse_args()

    comps = sorted(d for d in os.listdir(_LIB)
                   if os.path.isdir(os.path.join(_LIB, d))
                   and glob.glob(os.path.join(_LIB, d, "*_evaluation_report.md")))
    ok = leaked = skipped = 0
    print(f"{'compound':16s}{'studies':8s}{'build':7s}{'held':6s}{'in-vitro':9s}status")
    for c in comps:
        r = build_one(c, args.write)
        if r is None:
            skipped += 1
            continue
        if r["leaks"]:
            leaked += 1
            print(f"{c:16s}{r['n_studies']:<8d}{r['n_build']:<7d}{r['n_held']:<6d}"
                  f"{'yes' if r['invitro'] else 'no':9s}LEAK -> {r['leaks']}")
        else:
            ok += 1
            print(f"{c:16s}{r['n_studies']:<8d}{r['n_build']:<7d}{r['n_held']:<6d}"
                  f"{'yes' if r['invitro'] else 'no':9s}{'WROTE' if args.write else 'ok'}")
    print(f"\n{ok} ok · {leaked} leaked (skipped write) · {skipped} no report/data"
          + (f"  ->  {_OUT}" if args.write else "   (dry run; pass --write)"))


if __name__ == "__main__":
    if _PKG not in sys.path:
        sys.path.insert(0, _PKG)
    main()
