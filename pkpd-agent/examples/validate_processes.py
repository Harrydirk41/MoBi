r"""Validate the addable process-type catalog against a real PK-Sim.

For each single-compound mechanism in ``osp_catalog.PROCESS_TYPES`` this adds it
to a suitable expressed molecule (or as a systemic process), builds the snapshot
with PKSim.CLI, and runs one simulation - reporting which mechanisms BUILD
cleanly. Use it to flip the ``validated`` flags in osp_catalog in one pass on a
machine that has PK-Sim.

    set PKPD_PKSIM_CLI=C:\Program Files\Open Systems Pharmacology\PK-Sim 12.3\PKSim.CLI.exe
    python -m examples.validate_processes ^
        --snapshot ..\OSP-PBPK-Model-Library\Alfentanil\benchmark\Alfentanil-Model.blanked.json

Options:
    --only <type,type>   validate just these process types
    --unvalidated        validate only the types currently marked validated=false
    --keep               keep the PK-Sim work folders for inspection

DDI / interaction mechanisms (inhibition, induction) are NOT tested here: they
need a multi-compound Interactions setup, not a single-compound add. They are
listed at the end as "needs DDI setup".
"""

from __future__ import annotations

import argparse
import json
import os

from pkpd_agent.config import AgentConfig
from pkpd_agent.engines.osp_cli import OSPCli
from pkpd_agent.engines import osp_catalog


def _expressed(snapshot_path: str) -> dict[str, str]:
    with open(snapshot_path, encoding="utf-8") as fh:
        data = json.load(fh)
    return {ep.get("Molecule"): (ep.get("Type") or "").lower()
            for ep in data.get("ExpressionProfiles") or [] if ep.get("Molecule")}


def _has_process_for(snapshot_path: str) -> set[str]:
    """molecules that already carry a process (avoid 'already_present')."""
    with open(snapshot_path, encoding="utf-8") as fh:
        comp = (json.load(fh).get("Compounds") or [{}])[0]
    return {p.get("Molecule") for p in comp.get("Processes") or [] if p.get("Molecule")}


def _pick_molecule(spec: dict, expressed: dict[str, str], taken: set[str]):
    """A free expressed molecule of the tier the process needs (None for system)."""
    at = spec["applies_to"]
    if at == "system":
        return None, True
    want = {"enzyme": "enzyme", "transporter": "transporter",
            "target": None}[at]
    for mol, mtype in expressed.items():
        if mol in taken:
            continue
        if want is None or mtype == want:
            return mol, True
    return None, False


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--snapshot", required=True)
    ap.add_argument("--pksim", default=None)
    ap.add_argument("--only", default=None, help="comma-separated process types")
    ap.add_argument("--unvalidated", action="store_true")
    ap.add_argument("--keep", action="store_true")
    args = ap.parse_args()

    cfg = AgentConfig(mock=False)
    cli = OSPCli(pksim_cli_path=args.pksim or cfg.pksim_cli_path,
                 timeout_s=cfg.pksim_timeout_s, keep_workdir=args.keep)
    if not cli.pksim_cli_path or not os.path.exists(cli.pksim_cli_path):
        print(f"PKSim.CLI not found at {cli.pksim_cli_path!r}; "
              "set PKPD_PKSIM_CLI or pass --pksim.")
        return

    expressed = _expressed(args.snapshot)
    already = _has_process_for(args.snapshot)
    sims = cli.simulation_names(args.snapshot)
    one = sims[:1]                      # a single simulation is enough to confirm a build
    print(f"expressed molecules: {expressed}")
    print(f"validating against 1 simulation: {one}\n")

    types = list(osp_catalog.PROCESS_TYPES)
    if args.only:
        types = [t.strip() for t in args.only.split(",")]
    elif args.unvalidated:
        types = [t for t, s in osp_catalog.PROCESS_TYPES.items()
                 if not s.get("validated")]

    results = []
    for t in types:
        spec = osp_catalog.PROCESS_TYPES.get(t)
        if not spec:
            results.append((t, "SKIP", "unknown type")); continue
        # each type is validated in its OWN independent build, so they may reuse
        # the same free molecule; only avoid molecules that already carry a
        # process in the base snapshot (which would report 'already_present').
        mol, ok = _pick_molecule(spec, expressed, set(already))
        if not ok:
            results.append((t, "SKIP", f"no free expressed {spec['applies_to']}"))
            continue
        add = {"type": t}
        if mol:
            add["molecule"] = mol
        was = "valid" if spec.get("validated") else "UNVAL"
        print(f"[{was}] {t}  ->  {mol or '(systemic)'} ... ", end="", flush=True)
        res = cli.build_and_run(args.snapshot, edits={"add_processes": [add]},
                                simulations=one, prune_simulations=True)
        if res["ok"] and res.get("profiles"):
            print("BUILDS \u2713")
            results.append((t, "BUILDS", f"on {mol or 'systemic'}"))
        else:
            msg = (res.get("message") or "no result").splitlines()[0][:120]
            print("FAILED \u2717")
            results.append((t, "FAILED", msg))

    # summary + suggested flag flips
    print("\n================ SUMMARY ================")
    flip_on, keep_off = [], []
    for t, status, detail in results:
        cur = osp_catalog.PROCESS_TYPES.get(t, {}).get("validated")
        print(f"  {status:7} {t:34} {detail}")
        if status == "BUILDS" and not cur:
            flip_on.append(t)
        if status == "FAILED":
            keep_off.append(t)
    if flip_on:
        print("\nSet validated=True for (confirmed to build):")
        for t in flip_on:
            print(f"    osp_catalog.PROCESS_TYPES['{t}']['validated'] = True")
    if keep_off:
        print("\nLeave validated=False / fix these (failed to build):")
        for t in keep_off:
            print(f"    {t}")

    ddi = [r["type"] for r in osp_catalog.interaction_process_types()]
    print(f"\nDDI / interaction mechanisms (need a multi-compound setup, not "
          f"tested here): {', '.join(ddi)}")


if __name__ == "__main__":
    main()
