r"""Harness-fidelity check: does re-running each REFERENCE model reproduce its
PUBLISHED GMFE?

Every OSP evaluation report states the reference model's own goodness-of-fit GMFE
(see ``examples.published_gmfe``). A faithful harness, given the SAME reference
snapshot and the SAME observed data, must re-run it to roughly that GMFE. If the
harness score is far from the published one, the harness is not modelling the
reference correctly - and then every agent-vs-reference comparison in the scoreboard
is meaningless. This is the prerequisite check for the whole benchmark.

For each model it: loads the reference snapshot + the observed data from the task
input, runs the reference through the SAME path the report uses
(``cli.build_and_run`` -> ``map_predictions`` -> ``score_fit``), and compares the
harness GMFE to the published one. Requires PK-Sim (``--pksim``); without it the
script still prints the published column so the yardstick is visible.

    python -m examples.check_harness_fidelity --pksim "<PKSim.CLI.exe>"
    python -m examples.check_harness_fidelity            # published column only
"""

from __future__ import annotations

import argparse
import glob
import json
import os

from pkpd_agent.engines import osp_score
from pkpd_agent.engines.osp_cli import OSPCli

_HERE = os.path.dirname(os.path.abspath(__file__))
_LIB = os.path.abspath(os.path.join(_HERE, "..", "..", "OSP-PBPK-Model-Library"))


def _observed_for(model_dir: str, base: str) -> list | None:
    """The observed clinical data the reference is scored against (from the task
    input - the same data the agent sees)."""
    for suffix in (".input.json", ".hard.input.json"):
        p = os.path.join(model_dir, "json_input", base + suffix)
        if os.path.exists(p):
            gd = json.load(open(p, encoding="utf-8")).get("given_data", {}) or {}
            obs = gd.get("clinical_observed_data")
            if obs:
                return obs
    return None


def _reference_gmfe(cli: OSPCli, ref_path: str, observed: list) -> "float | None":
    res = cli.build_and_run(ref_path)
    if not res.get("ok"):
        return None
    with open(ref_path, encoding="utf-8") as fh:
        linkage = osp_score.linkage_from_snapshot(json.load(fh))
    pred, _ = osp_score.map_predictions(res.get("profiles", []), observed, linkage=linkage)
    return osp_score.score_fit(observed, pred)["overall"]["gmfe"]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pksim", default=os.environ.get("PKPD_PKSIM_CLI"),
                    help="path to PKSim.CLI.exe (or set PKPD_PKSIM_CLI)")
    ap.add_argument("--models", default=None, help="comma-separated subset")
    ap.add_argument("--tol", type=float, default=1.25,
                    help="max harness/published ratio (either direction) still OK")
    args = ap.parse_args()

    pub_path = os.path.join(_LIB, "published_gmfe.json")
    if not os.path.exists(pub_path):
        raise SystemExit("run `python -m examples.published_gmfe` first to create "
                         "published_gmfe.json")
    published = json.load(open(pub_path, encoding="utf-8"))
    cli = OSPCli(pksim_cli_path=args.pksim) if args.pksim else None

    want = {m.strip() for m in args.models.split(",")} if args.models else None
    # ONE primary reference per model dir - the PRIMARY (non-pediatric) model, the one
    # the benchmark actually grades. The published GMFE parsed from the report is the
    # ADULT model's, so a pediatric variant compared against it is the wrong yardstick
    # (a false OFF), and the scoreboard never uses pediatric models anyway.
    by_dir: dict[str, str] = {}
    for ref in sorted(glob.glob(os.path.join(_LIB, "*", "json", "*.json"))):
        d = os.path.dirname(os.path.dirname(ref))
        cur = by_dir.get(d)
        if cur is None or ("pediatric" in os.path.basename(cur).lower()
                           and "pediatric" not in os.path.basename(ref).lower()):
            by_dir[d] = ref
    refs = [by_dir[d] for d in sorted(by_dir)]
    print(f"{'model':22} {'published':>9} {'harness':>8} {'ratio':>6}  verdict")
    n_ok = n_off = n_run = 0
    for ref in refs:
        model_dir = os.path.dirname(os.path.dirname(ref))
        model = os.path.basename(model_dir)
        base = os.path.splitext(os.path.basename(ref))[0]
        if want and model not in want:
            continue
        pub = (published.get(model) or {}).get("all")
        if pub is None:
            continue
        hg = None
        if cli is not None:
            observed = _observed_for(model_dir, base)
            if observed:
                try:
                    hg = _reference_gmfe(cli, ref, observed)
                except Exception as e:                      # noqa: BLE001
                    hg = None
                    print(f"{model:22} {pub:>9} {'ERR':>8}   {str(e)[:40]}")
                    continue
        if hg is None:
            print(f"{model:22} {pub:>9} {'-':>8} {'-':>6}  (no run)")
            continue
        n_run += 1
        ratio = hg / pub if pub else float("inf")
        ok = (1 / args.tol) <= ratio <= args.tol
        n_ok += ok; n_off += (not ok)
        print(f"{model:22} {pub:>9} {hg:>8.3g} {ratio:>6.2f}  {'OK' if ok else 'OFF - harness misfits the reference'}")
    if cli is not None:
        print(f"\n{n_ok}/{n_run} references reproduced within {args.tol}x of published; "
              f"{n_off} OFF (harness fidelity problems to fix before trusting agent-vs-reference).")
    else:
        print("\n(no --pksim: published yardstick shown only; pass --pksim to run the "
              "references and fill the harness column)")


if __name__ == "__main__":
    main()
