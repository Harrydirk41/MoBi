r"""Predict a victim DDI ratio from a multi-compound OSP snapshot (headless).

For a perpetrator+victim model (e.g. Erythromycin + Midazolam) this:
  1. reads the DDI structure (perpetrator, victim, mechanisms, control/treatment
     pairs) with osp_ddi.analyze_ddi;
  2. runs each control (victim alone) and treatment (victim + perpetrator) arm
     with PKSim.CLI, extracting the VICTIM's plasma;
  3. computes the predicted interaction ratio AUCR = AUC_treatment / AUC_control
     - e.g. does the model reproduce erythromycin ~quadrupling midazolam AUC.

Optionally tune the perpetrator's interaction parameters and re-check, and score
against observed ratios.

    set PKPD_PKSIM_CLI=C:\Program Files\Open Systems Pharmacology\PK-Sim 12.3\PKSim.CLI.exe
    python -m examples.run_ddi ^
        --snapshot ..\OSP-PBPK-Model-Library\Erythromycin\json\Erythromycin-Model.json ^
        --victim Midazolam
"""

from __future__ import annotations

import argparse
import json
import os

from pkpd_agent.config import AgentConfig
from pkpd_agent.engines.osp_cli import OSPCli
from pkpd_agent.engines import osp_ddi


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--snapshot", required=True)
    ap.add_argument("--victim", default=None, help="victim compound (default: inferred)")
    ap.add_argument("--pksim", default=None)
    args = ap.parse_args()

    with open(args.snapshot, encoding="utf-8") as fh:
        snap = json.load(fh)
    ddi = osp_ddi.analyze_ddi(snap)
    if not ddi:
        print("Not a DDI snapshot (need >1 compound + a perpetrator with "
              "inhibition/induction and treatment Interactions).")
        return

    victim = args.victim or (ddi["victims"][0] if ddi["victims"] else None)
    print(f"perpetrator(s): {[p['name'] for p in ddi['perpetrators']]}")
    for p in ddi["perpetrators"]:
        for m in p["mechanisms"]:
            print(f"    {p['name']}: {m['internal_name']} on {m['target']} "
                  f"({m['data_source']}) {m['parameters']}")
    print(f"victim: {victim}")
    print(f"control/treatment pairs: {len(ddi['pairs'])}")
    for pr in ddi["pairs"]:
        print(f"    {pr['control']}  ->  {pr['treatment']}  [{pr['route']}]")
    print(f"(unpaired arms - e.g. the perpetrator's own PK: {len(ddi['unpaired_treatments'])})\n")

    cfg = AgentConfig(mock=False)
    cli = OSPCli(pksim_cli_path=args.pksim or cfg.pksim_cli_path,
                 timeout_s=cfg.pksim_timeout_s)
    if not cli.pksim_cli_path or not os.path.exists(cli.pksim_cli_path):
        print(f"PKSim.CLI not found at {cli.pksim_cli_path!r}; set PKPD_PKSIM_CLI. "
              "(Structure above was read without it.)")
        return

    print("running control + treatment arms (extracting the victim's plasma)...",
          flush=True)
    out = osp_ddi.run_ddi_prediction(cli, args.snapshot, ddi, victim)
    if not out["ok"]:
        print(f"failed: {out['message']}")
        return
    print(f"\nran {out['n_ran']} simulations for {out['n_pairs']} pairs\n")
    print("Predicted interaction ratios (victim exposure with vs without perpetrator):")
    for r in out["predicted_ratios"]:
        aucr = r.get("auc_ratio")
        cmaxr = r.get("cmax_ratio")
        aucr_s = f"AUCR={aucr:.2f}" if isinstance(aucr, (int, float)) else "AUCR=n/a"
        cmaxr_s = f"CmaxR={cmaxr:.2f}" if isinstance(cmaxr, (int, float)) else ""
        print(f"    {r['treatment']:34} [{r['route']}]  {aucr_s}  {cmaxr_s}")


if __name__ == "__main__":
    main()
