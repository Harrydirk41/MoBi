r"""One-shot coverage matrix: which benchmark task each library snapshot supports.

Runs all four benchmark builders over the OSP library and reports, per snapshot,
which task type(s) it can generate - single-compound, metabolite cascade, DDI
(victim), or biologic - plus the parameters/quantities each would recover. This
is the capability map: what the agent benchmark now covers end-to-end.

    python -m examples.benchmark_coverage ..\OSP-PBPK-Model-Library
"""

from __future__ import annotations

import argparse
import glob
import json
import os

from pkpd_agent.engines import osp_ddi, osp_metabolite, osp_biologic

import examples.build_biologic_benchmark as BIO
import examples.build_metabolite_benchmark as MET
import examples.build_ddi_benchmark as DDI


def _single_ok(snap: dict) -> bool:
    """A single-compound recovery task: one small molecule, fitted params, not a
    DDI/metabolite/biologic."""
    comps = snap.get("Compounds") or []
    if len(comps) != 1 or comps[0].get("IsSmallMolecule") is False:
        return False
    if any(s.get("Interactions") for s in snap.get("Simulations") or []):
        return False

    def fitted(o, acc):
        if isinstance(o, dict):
            if (o.get("ValueOrigin") or {}).get("Source") == "ParameterIdentification" \
                    and isinstance(o.get("Value"), (int, float)):
                acc.append(1)
            for v in o.values():
                fitted(v, acc)
        elif isinstance(o, list):
            for v in o:
                fitted(v, acc)
    acc = []
    fitted(comps[0], acc)
    return bool(acc)


def classify_coverage(path: str) -> dict:
    with open(path, encoding="utf-8") as fh:
        snap = json.load(fh)
    tasks = []
    if _single_ok(snap):
        tasks.append("single")
    m = osp_metabolite.analyze_metabolites(snap)
    if m and (m.get("scorable_molecules")):
        tasks.append(f"metabolite({len(m['scorable_molecules'])})")
    d = osp_ddi.analyze_ddi(snap)
    if d and d.get("victims") and d.get("pairs"):
        r = DDI.build(path)
        if not r.get("skip"):
            tasks.append(f"ddi({r['n_pairs']}p)")
    b = osp_biologic.analyze_biologic(snap)
    if b and b.get("disposition_parameters"):
        tasks.append(f"biologic({len(b['observed_matrices'])}m)")
    return {"tasks": tasks}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("target", nargs="?", default="../OSP-PBPK-Model-Library")
    args = ap.parse_args()

    files = ([args.target] if os.path.isfile(args.target)
             else sorted(glob.glob(os.path.join(args.target, "*/json/*.json"))))
    covered, uncovered = [], []
    counts: dict[str, int] = {}
    print(f"{'snapshot':40} tasks")
    print("-" * 70)
    for f in files:
        stem = os.path.basename(f)
        try:
            c = classify_coverage(f)
        except Exception as exc:                       # noqa: BLE001
            print(f"{stem:40} ERROR {type(exc).__name__}: {exc}")
            continue
        tasks = c["tasks"]
        print(f"{stem:40} {', '.join(tasks) if tasks else '-'}")
        (covered if tasks else uncovered).append(stem)
        for t in tasks:
            key = t.split("(")[0]
            counts[key] = counts.get(key, 0) + 1
    print("-" * 70)
    print(f"covered: {len(covered)}/{len(covered) + len(uncovered)} snapshots")
    print(f"by task type: {counts}")
    if uncovered:
        print(f"not covered ({len(uncovered)}): {', '.join(uncovered)}")


if __name__ == "__main__":
    main()
