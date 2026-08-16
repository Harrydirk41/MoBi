r"""End-to-end validation of the metabolite / biologic / DDI task engines on PK-Sim.

The structure analyzers and scorers are unit-tested without PK-Sim, but three
things can only be confirmed on a machine that HAS PK-Sim, because they depend on
the exact Results.csv column names PK-Sim emits:

  * metabolite - are the PARENT and each DAUGHTER metabolite pulled from the one
    run (per-molecule column match)?
  * biologic   - do the (organ, compartment) column TOKENS match the real
    biodistribution column names (whole blood + each tissue)?  [main unknown]
  * ddi        - does the victim's plasma column resolve in a multi-compound run?

This runs one real snapshot per task type, reports how many molecules/matrices/
arms were matched, and prints the FIRST result CSV header so any column-name
mismatch is visible. Run it once on a PK-Sim box to confirm (or to read off the
exact tissue column names if a biologic matrix comes back unmatched).

    set PKPD_PKSIM_CLI=C:\Program Files\Open Systems Pharmacology\PK-Sim 12.3\PKSim.CLI.exe
    python -m examples.validate_new_tasks --lib ..\OSP-PBPK-Model-Library
    python -m examples.validate_new_tasks --only biologic --keep   # inspect columns
"""

from __future__ import annotations

import argparse
import glob
import json
import os

from pkpd_agent.config import AgentConfig
from pkpd_agent.engines.osp_cli import OSPCli
from pkpd_agent.engines import osp_metabolite, osp_biologic, osp_ddi


def _first(lib: str, compound: str) -> str | None:
    fs = sorted(glob.glob(os.path.join(lib, compound, "json", "*.json")))
    return fs[0] if fs else None


def _dump_header(cli: OSPCli, snapshot_path: str, sims: list | None) -> None:
    """Run once keeping the workdir, print the first Results.csv header so the
    real column names are visible (to confirm/fix the token matchers)."""
    cli2 = OSPCli(pksim_cli_path=cli.pksim_cli_path, timeout_s=cli.timeout_s,
                  keep_workdir=True)
    res = cli2.build_and_run(snapshot_path, simulations=sims,
                             prune_simulations=bool(sims))
    proj = res.get("project")
    if not proj:
        print("   (no project built - cannot show header)")
        return
    out_dir = os.path.join(os.path.dirname(os.path.dirname(proj)), "out")
    csvs = sorted(glob.glob(os.path.join(out_dir, "*", "*-Results.csv")))
    if csvs:
        with open(csvs[0], encoding="utf-8", errors="replace") as fh:
            header = fh.readline().strip()
        print(f"   first CSV columns:\n     " + header.replace(",", "\n     "))


def val_metabolite(cli, lib, keep):
    path = _first(lib, "Itraconazole")
    if not path:
        print("metabolite: Itraconazole snapshot not found"); return
    with open(path, encoding="utf-8") as fh:
        snap = json.load(fh)
    m = osp_metabolite.analyze_multicompound(snap)
    print(f"metabolite: {path}\n   cascade {' -> '.join(m['chain'])}, "
          f"scorable {m['scorable_molecules']}")
    out = osp_metabolite.run_metabolite_prediction(cli, path, m, snapshot=snap)
    if not out["ok"]:
        print(f"   FAILED: {out['message']}")
        if keep:
            _dump_header(cli, path, None)
        return
    for mol, v in out["score"]["per_molecule"].items():
        print(f"   {mol:28} GMFE {v['overall'].get('gmfe')} "
              f"(matched {v['n_matched']}/{v['n_datasets']} datasets)")
    print(f"   cascade GMFE: {(out['score'].get('cascade') or {}).get('gmfe')}")


def val_biologic(cli, lib, keep):
    path = _first(lib, "BAY794620")
    if not path:
        print("biologic: BAY794620 snapshot not found"); return
    with open(path, encoding="utf-8") as fh:
        snap = json.load(fh)
    b = osp_biologic.analyze_biologic(snap)
    n_specs = len(osp_biologic.matrix_specs(
        osp_biologic.biologic_observed(snap, b["molecule"]), b["molecule"]))
    print(f"biologic: {path}\n   {b['molecule']}, "
          f"{len(b['observed_matrices'])} organ/compartment matrices "
          f"({n_specs} study-split), recover "
          f"{[p['name'] for p in b['disposition_parameters']]}")
    out = osp_biologic.run_biologic_prediction(cli, path, b, snapshot=snap)
    if not out["ok"]:
        print(f"   FAILED: {out['message']}")
        if keep:
            _dump_header(cli, path, None)
        return
    sc = out["score"]
    print(f"   matched {out['n_matched']}/{out['n_matrices']} matrices, "
          f"overall GMFE {sc['overall'].get('gmfe')}")
    # per-matrix breakdown (bias > 1 = over-predicts): a SYSTEMATIC bias across
    # all matrices points to a quantity/basis mismatch; error concentrated in the
    # low-penetration tissues (fat/brain/muscle) is normal mAb-PBPK behaviour.
    rows = [(k, v) for k, v in sc["per_matrix"].items() if v["matched"]]
    rows.sort(key=lambda kv: (kv[1].get("gmfe") or 0), reverse=True)
    print(f"   {'matrix':24} {'GMFE':>7} {'bias':>7}")
    for k, v in rows:
        print(f"   {k:24} {str(v.get('gmfe')):>7} {str(v.get('bias')):>7}")
    unmatched = [k for k, v in sc["per_matrix"].items() if not v["matched"]]
    if unmatched:
        print(f"   UNMATCHED (no model organ - e.g. Tumor/Bone/Intestine have no "
              f"standard PK-Sim compartment): {unmatched}")
        if keep:
            _dump_header(cli, path, None)


def val_ddi(cli, lib, keep):
    path = _first(lib, "Erythromycin")
    if not path:
        print("ddi: Erythromycin snapshot not found"); return
    with open(path, encoding="utf-8") as fh:
        snap = json.load(fh)
    d = osp_ddi.analyze_ddi(snap)
    victim = d["victims"][0]
    print(f"ddi: {path}\n   perpetrator -> {victim}, {len(d['pairs'])} pairs")
    out = osp_ddi.run_ddi_prediction(cli, path, d, victim)
    if not out["ok"]:
        print(f"   FAILED: {out['message']}"); return
    for r in out["predicted_ratios"]:
        print(f"   {r['treatment']:34} AUCR {r.get('auc_ratio')}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lib", default="../OSP-PBPK-Model-Library")
    ap.add_argument("--pksim", default=None)
    ap.add_argument("--only", default=None,
                    help="metabolite | biologic | ddi (default: all three)")
    ap.add_argument("--keep", action="store_true",
                    help="on a miss, print the real CSV column header")
    args = ap.parse_args()

    cfg = AgentConfig(mock=False)
    cli = OSPCli(pksim_cli_path=args.pksim or cfg.pksim_cli_path,
                 timeout_s=cfg.pksim_timeout_s)
    if not cli.pksim_cli_path or not os.path.exists(cli.pksim_cli_path):
        print(f"PKSim.CLI not found at {cli.pksim_cli_path!r}; set PKPD_PKSIM_CLI "
              "or pass --pksim. (This validation needs a machine with PK-Sim.)")
        return

    runners = {"metabolite": val_metabolite, "biologic": val_biologic, "ddi": val_ddi}
    for name, fn in runners.items():
        if args.only and args.only != name:
            continue
        print("=" * 70)
        try:
            fn(cli, args.lib, args.keep)
        except Exception as exc:                       # noqa: BLE001
            print(f"{name}: ERROR {type(exc).__name__}: {exc}")
    print("=" * 70)


if __name__ == "__main__":
    main()
