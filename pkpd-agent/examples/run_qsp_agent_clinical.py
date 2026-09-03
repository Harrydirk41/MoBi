r"""Drive the agent-network -> clinical-sbproj bridge from ONE PowerShell/terminal command.

Just like run_qsp_paper_pipeline's PART II, this starts MATLAB headlessly through the SimBiology
engine (matlab.engine) and calls the sb_*.m helpers for you - so you never type MATLAB commands
yourself. It transplants the agent-built immune network (mynet.xml from run_qsp_build_network) into
your real .sbproj's immune block, keeping the paper's given DAS28/PK/dose shell.

    # 0. emit the agent network (needs your ANTHROPIC_API_KEY):
    python -m examples.run_qsp_build_network --model ra --live --prune --emit mynet.xml

    # 1. DRY RUN (default) - prints the transplant report, changes nothing:
    python -m examples.run_qsp_agent_clinical --sbproj "Vantage RA QSP Model v1.0.sbproj" ^
        --net mynet.xml

    # 2. APPLY - transplant + baseline sanity sim + save the agent-based sbproj:
    python -m examples.run_qsp_agent_clinical --sbproj "Vantage RA QSP Model v1.0.sbproj" ^
        --net mynet.xml --apply --out agent_clinical.sbproj

Then run your usual train/test/simulate on agent_clinical.sbproj (run_qsp_paper_pipeline --sbproj
agent_clinical.sbproj --vpop ... --matlab, or sb_fit / sb_run_vpop). Needs MATLAB + SimBiology +
the matlab.engine Python package installed (the same setup run_qsp_paper_pipeline --matlab uses).
"""

from __future__ import annotations

import argparse
import json
import os


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sbproj", required=True, help="path to the paper's real .sbproj")
    ap.add_argument("--net", default="mynet.xml", help="the agent-built SBML (run_qsp_build_network)")
    ap.add_argument("--apply", action="store_true",
                    help="actually transplant + save (default is a dry run that changes nothing)")
    ap.add_argument("--out", default="agent_clinical.sbproj",
                    help="where to save the agent-based clinical sbproj (with --apply)")
    ap.add_argument("--readout", default="DAS28_CRP", help="readout to sanity-check after apply")
    args = ap.parse_args()

    for p in (args.sbproj, args.net):
        if not os.path.isfile(p):
            raise SystemExit(f"file not found: {p}")

    from pkpd_agent.engines.simbiology import SimBiologyEngine
    sb = SimBiologyEngine()
    print("== starting MATLAB (headless) ==", flush=True)
    sb.start()
    try:
        if not sb.has_simbiology():
            raise SystemExit("MATLAB started but SimBiology is not licensed here.")
        if not args.apply:
            print("== DRY RUN: loading sbproj + transplant report (nothing changed) ==", flush=True)
            sb.load_project(os.path.abspath(args.sbproj))
            rep = sb.eng.sb_transplant_immune(os.path.abspath(args.net), True, nargout=1)
            # the .m already prints a human report via the captured MATLAB stdout; also dump the
            # struct fields so you (and I) can eyeball the seams programmatically.
            print("\n== transplant report (struct) ==")
            for k in ("immuneShared", "uncovered", "shellSpeciesInRemoved", "clinicalCouplings",
                      "removeMixedDetails"):
                v = rep.get(k) if isinstance(rep, dict) else None
                print(f"  {k}: {json.dumps(v, default=str)[:600] if v is not None else '(n/a)'}")
            print("\nDRY RUN only. Re-run with --apply once the report looks right.")
        else:
            print("== APPLY: transplant + baseline sim + save agent-based sbproj ==", flush=True)
            sb.eng.sb_agent_clinical(os.path.abspath(args.sbproj), os.path.abspath(args.net),
                                     os.path.abspath(args.out), args.readout, "", nargout=0)
            print(f"\n  saved -> {os.path.abspath(args.out)}")
            print("  next: run_qsp_paper_pipeline --sbproj <that> --vpop Vpop1.xlsx --matlab, "
                  "or sb_fit / sb_run_vpop for train / test / simulate.")
    finally:
        sb.stop()


if __name__ == "__main__":
    main()
