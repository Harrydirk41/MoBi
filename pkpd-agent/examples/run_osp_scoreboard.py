r"""Run the PBPK build agent across many library models and collect ONE honest scoreboard.

For each model it auto-resolves the four files from the benchmark base name (no hand-typing):
  benchmark/<base>.blanked.json   -> --snapshot   (agent starts here)
  json_input/<base>.input.json    -> --input      (the agent-facing task)
  json/<base>.json                -> --reference   (ground truth, for the comparison)
  answer_key/<base>.answer_edits.json -> --answer-edits (parameter comparison)
and writes the report (with its machine-readable .json) into <model>/report/.

Two phases (default: both):
  --run        execute run_llm_build on each model (needs PK-Sim CLI + ANTHROPIC_API_KEY)
  --aggregate  read every model's report/*.json and print+write one markdown scoreboard,
               ordered by the difficulty rank below, showing agent vs reference GMFE,
               structure match, and parameter recovery (good / soft / bad).

    python -m examples.run_osp_scoreboard --run --aggregate ^
        --models Vancomycin,Tizanidine,Sufentanil,Montelukast,Raltegravir,Mexiletine,^
                 Alfentanil,Sildenafil,Midazolam,Digoxin --target 1.6 --max-steps 6

Difficulty ranking is an INFORMED ESTIMATE (mechanism complexity + #unknowns + data breadth),
not a measured truth - the scoreboard is exactly what tests whether that ordering holds.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import subprocess
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_LIB = os.path.abspath(os.path.join(_HERE, "..", "..", "OSP-PBPK-Model-Library"))

# easy -> hard (informed estimate; the run is what checks it)
DEFAULT_ORDER = ["Vancomycin", "Tizanidine", "Sufentanil", "Montelukast", "Raltegravir",
                 "Mexiletine", "Alfentanil", "Sildenafil", "Midazolam", "Digoxin"]


def _resolve(model: str, hard: bool = False) -> "dict | None":
    """Find the primary (non-pediatric) benchmark base for a model dir and resolve its files.

    ``hard=True`` resolves the mechanism-DISCOVERY variant (``*.hard_blanked.json`` +
    ``*.hard.input.json``): structure stripped and the metabolizing enzyme withheld, so
    the agent must build the model and discover the mechanism. The answer key and reference
    are shared; the report is written to a separate ``*.hard.*`` path so it never clobbers
    the easy-mode report."""
    d = os.path.join(_LIB, model)
    suffix = ".hard_blanked.json" if hard else ".blanked.json"
    # exclude the OTHER variant's snapshots from the glob
    cand = [p for p in glob.glob(os.path.join(d, "benchmark", "*" + suffix))
            if hard or not p.endswith(".hard_blanked.json")]
    blanked = sorted(cand, key=lambda p: ("pediatric" in p.lower(), len(p)))
    for b in blanked:
        base = os.path.basename(b)[: -len(suffix)]
        inp = os.path.join(d, "json_input", base + (".hard.input.json" if hard else ".input.json"))
        ref = os.path.join(d, "json", base + ".json")
        ans = os.path.join(d, "answer_key", base + ".answer_edits.json")
        if all(os.path.isfile(p) for p in (inp, ref, ans)):
            tag = ".hard" if hard else ""
            return {"model": model, "base": base, "snapshot": b, "input": inp,
                    "reference": ref, "answer": ans, "hard": hard,
                    "report": os.path.join(d, "report", base + tag + ".html"),
                    "report_json": os.path.join(d, "report", base + tag + ".json")}
    return None


def _run_one(f: dict, target: float, max_steps: int, pksim: "str | None") -> bool:
    cmd = [sys.executable, "-m", "examples.run_llm_build",
           "--snapshot", f["snapshot"], "--input", f["input"],
           "--reference", f["reference"], "--answer-edits", f["answer"],
           "--report", f["report"], "--target", str(target), "--max-steps", str(max_steps)]
    if pksim:
        cmd += ["--pksim", pksim]
    print(f"\n{'='*70}\n== RUN {f['model']} ({f['base']}) ==\n{'='*70}", flush=True)
    try:
        return subprocess.run(cmd, cwd=os.path.join(_HERE, "..")).returncode == 0
    except Exception as e:                              # noqa: BLE001
        print(f"  {f['model']} FAILED to launch: {e}")
        return False


def _cell(x, nd=3):
    return "-" if x is None else (f"{x:.{nd}g}" if isinstance(x, (int, float)) else str(x))


def _aggregate(files: list, order: list, out: str) -> None:
    rows = []
    for f in files:
        j = f.get("report_json")
        if not (j and os.path.isfile(j)):
            rows.append((f["model"], None)); continue
        try:
            rows.append((f["model"], json.load(open(j, encoding="utf-8"))))
        except Exception:                              # noqa: BLE001
            rows.append((f["model"], None))
    rank = {m: i + 1 for i, m in enumerate(order)}
    rows.sort(key=lambda r: rank.get(r[0], 99))

    hard = any(f.get("hard") for f in files)
    struct_note = ("Structure = did the agent DISCOVER the reference's mechanism "
                   "(metabolizing enzyme + methods) - it was withheld in this mode."
                   if hard else
                   "Structure = did the agent pick the reference's methods/processes.")
    L = [f"# PBPK agent scoreboard{' - HARD (mechanism discovery)' if hard else ''} "
         "(agent vs ground-truth model)", "",
         "PRIMARY grade is HELD-OUT: the agent fits the building studies and is graded on "
         "how well it PREDICTS the held-out verification studies, RELATIVE to the reference "
         f"model on the same held-out data (PASS = agent held-out GMFE <= 1.25x reference's). "
         f"{struct_note} 'all-data GMFE' is the in-sample fit, shown for context.", "",
         "| # | model | held-out agent | held-out ref | ratio | structure | all-data agent | verdict |",
         "|---|---|---|---|---|---|---|---|"]
    for i, (model, d) in enumerate(rows, 1):
        if not d:
            L.append(f"| {i} | {model} | - | - | - | - | - | (no report - run failed / not run) |")
            continue
        fit = d.get("fit") or {}
        ho = d.get("heldout") or {}
        sm = d.get("structure_match")
        smtxt = "match" if sm is True else ("differs" if sm is False else "-")
        ag, rg = fit.get("agent_gmfe"), fit.get("reference_gmfe")
        a_ho, r_ho, ratio = ho.get("agent_gmfe"), ho.get("reference_gmfe"), ho.get("ratio")
        if ho and ho.get("pass") is not None:
            verdict = "PASS (predicts held-out ~= reference)" if ho.get("pass") \
                else "FAIL (worse than reference on held-out)"
        else:
            verdict = "no held-out split (graded on all data)"
        L.append(f"| {i} | {model} | {_cell(a_ho)} | {_cell(r_ho)} | "
                 f"{_cell(ratio,2) if ratio else '-'} | {smtxt} | {_cell(ag)} | {verdict} |")
    # honest aggregate - the primary bar is the HELD-OUT relative grade
    scored = [d for _, d in rows if d and isinstance((d.get('fit') or {}).get('agent_gmfe'), (int, float))]
    graded = [d for d in scored if (d.get("heldout") or {}).get("pass") is not None]
    n_pass = sum(1 for d in graded if d["heldout"]["pass"])
    n_none = len(scored) - len(graded)
    L += ["", f"**{len(scored)}/{len(rows)} models produced a scored run; of those, {len(graded)} had "
          f"a held-out split and {n_pass} PASSED (predicted held-out data within 1.25x of the "
          f"reference); {n_none} had too few studies to hold out (graded on all data).** Missing "
          f"rows are runs that did not finish - report them as gaps, not passes."]
    text = "\n".join(L)
    open(out, "w", encoding="utf-8").write(text)
    print("\n" + text + f"\n\n(scoreboard -> {os.path.abspath(out)})")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--models", default=",".join(DEFAULT_ORDER),
                    help="comma-separated model dir names, easy->hard")
    ap.add_argument("--run", action="store_true", help="execute the agent on each model")
    ap.add_argument("--aggregate", action="store_true", help="collect reports into a scoreboard")
    ap.add_argument("--target", type=float, default=1.6)
    ap.add_argument("--max-steps", type=int, default=6)
    ap.add_argument("--pksim", default=None)
    ap.add_argument("--hard", action="store_true",
                    help="run the mechanism-DISCOVERY variant (structure stripped, "
                         "metabolizing enzyme withheld) instead of the fill-in-the-blank one")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    if not (args.run or args.aggregate):
        args.run = args.aggregate = True
    if args.out is None:
        args.out = "osp_scoreboard_hard.md" if args.hard else "osp_scoreboard.md"

    order = [m.strip() for m in args.models.split(",") if m.strip()]
    files, missing = [], []
    for m in order:
        f = _resolve(m, hard=args.hard)
        (files.append(f) if f else missing.append(m))
    if missing:
        print(f"[skipped - no complete benchmark found]: {missing}")

    if args.run:
        ok = 0
        for f in files:
            ok += _run_one(f, args.target, args.max_steps, args.pksim)
            if args.aggregate:
                _aggregate(files, order, args.out)   # refresh the scoreboard after EACH model
        print(f"\n== ran {ok}/{len(files)} models ==")
    elif args.aggregate:
        _aggregate(files, order, args.out)


if __name__ == "__main__":
    main()
