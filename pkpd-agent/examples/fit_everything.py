r"""'Fit everything' baseline - the deterministic anti-pattern, for contrast.

Instead of the agent selecting a minimal identifiable set, this frees EVERY non-given
(placeholder) parameter at once with wide bounds and lets the optimizer fit them all to
the BUILDING data. No LLM, no judgement. It then scores the result on the building data
AND on the held-out verification data, next to the reference model, so you can see the
over-fitting directly: freeing too many collinear/unidentifiable parameters usually fits
the building data as well or better, but predicts the HELD-OUT data worse - the model
memorised the fit instead of capturing the physiology.

    python -m examples.fit_everything ^
        --snapshot ..\OSP-PBPK-Model-Library\Triazolam\benchmark\Triazolam-Model.blanked.json ^
        --input    ..\OSP-PBPK-Model-Library\Triazolam\json_input\Triazolam-Model.input.json ^
        --reference ..\OSP-PBPK-Model-Library\Triazolam\json\Triazolam-Model.json ^
        --pksim "C:\...\PKSim.CLI.exe"
"""

from __future__ import annotations

import argparse
import json
import os

from pkpd_agent.engines import osp_optimize as OO
from pkpd_agent.engines import osp_score, osp_split
from pkpd_agent.engines.osp_catalog import describe_parameter
from pkpd_agent.engines.osp_cli import OSPCli
from pkpd_agent.tools.osp_loop_tools import _current_model


def _wide_bounds(name: str) -> list:
    """A GENEROUS physical bound per parameter - the point is to over-free."""
    base = name.split("@", 1)[0]
    pb = osp_score.physical_bounds(base)
    if pb:
        return [pb[0] if pb[0] > 0 else 1e-4, pb[1]]
    rng = (describe_parameter(base) or {}).get("range")
    if rng and rng[1] and rng[1] > (rng[0] or 0):
        return [rng[0] if rng[0] and rng[0] > 0 else rng[1] * 1e-6, rng[1]]
    n = base.lower()
    if "dissolution time" in n:
        return [0.1, 300.0]
    if "dissolution shape" in n:
        return [0.1, 5.0]
    if any(k in n for k in ("clear", "vmax", "kcat", "clspec", "clint")):
        return [1e-3, 1e4]
    if "km" in n:
        return [1e-2, 1e4]
    if "permeab" in n:
        return [1e-8, 1e-2]
    if "solub" in n:
        return [1e-3, 1e5]
    return [1e-4, 1e4]


def _gmfe_on(cli, snapshot, edits, observed, linkage):
    res = cli.build_and_run(snapshot, edits=edits)
    if not res.get("ok"):
        return None
    pred, _ = osp_score.map_predictions(res.get("profiles", []), observed, linkage=linkage)
    return osp_score.score_fit(observed, pred)["overall"]["gmfe"] if pred else None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--snapshot", required=True)
    ap.add_argument("--input", required=True)
    ap.add_argument("--reference")
    ap.add_argument("--pksim", default=os.environ.get("PKPD_PKSIM_CLI"))
    ap.add_argument("--max-evals", type=int, default=80)
    args = ap.parse_args()

    cli = OSPCli(pksim_cli_path=args.pksim)
    with open(args.snapshot, encoding="utf-8") as fh:
        snap = json.load(fh)
    observed = json.load(open(args.input, encoding="utf-8"))["given_data"]["clinical_observed_data"]
    link = osp_score.linkage_from_snapshot(snap)
    split = osp_split.split_studies(snap, observed)
    bset, vset = set(split["building"]), set(split["verification"])
    build_obs = [o for o in observed if o["dataset"] in bset]
    verif_obs = [o for o in observed if o["dataset"] in vset]

    model = _current_model(args.snapshot)
    estimate, fix = {}, {}
    for p in model["parameters"]:
        if p["value_status"] == "placeholder":
            estimate[p["name"]] = _wide_bounds(p["name"])
        elif p["value_status"] == "given":
            fix[p["name"]] = p["value"]

    print(f"== FIT EVERYTHING on {os.path.basename(args.snapshot)} ==")
    print(f"   freeing {len(estimate)} parameters at once (vs a modeler's ~2-5):")
    for k, b in estimate.items():
        print(f"     {k:44} in [{b[0]:.3g}, {b[1]:.3g}]")
    print(f"   holding {len(fix)} given/measured fixed.")
    print(f"   building={len(build_obs)} datasets, held-out verification={len(verif_obs)}\n")

    r = OO.run_optimization(cli, args.snapshot, build_obs, estimate=estimate, fix=fix,
                            max_evals=args.max_evals)
    if not r.get("ok"):
        print("optimization failed:", r.get("message")); return
    best = {**fix, **r["optimized"]}
    edits = {"parameters": best}
    b_gmfe = _gmfe_on(cli, args.snapshot, edits, build_obs, link)
    v_gmfe = _gmfe_on(cli, args.snapshot, edits, verif_obs, link) if verif_obs else None
    rv_gmfe = None
    if args.reference and verif_obs:
        rlink = osp_score.linkage_from_snapshot(json.load(open(args.reference, encoding="utf-8")))
        rv_gmfe = _gmfe_on(cli, args.reference, None, verif_obs, rlink)

    print("\n=== result ===")
    print(f"   BUILDING (in-sample) GMFE : {b_gmfe}")
    print(f"   HELD-OUT  GMFE (this)     : {v_gmfe}")
    print(f"   HELD-OUT  GMFE (reference): {rv_gmfe}")
    if v_gmfe and rv_gmfe:
        print(f"   held-out ratio vs reference: {v_gmfe / rv_gmfe:.2f}  "
              f"({'worse' if v_gmfe > rv_gmfe * 1.25 else 'comparable'})")
    if b_gmfe and v_gmfe and v_gmfe > b_gmfe * 1.2:
        print("   -> OVER-FIT signature: fits building better than it predicts held-out.")
    # flag railed params (a symptom of unidentifiability)
    railed = r.get("params_at_bound") or []
    if railed:
        print(f"   params that RAILED to a bound (unidentifiable): "
              f"{[p.get('parameter') for p in railed]}")


if __name__ == "__main__":
    main()
