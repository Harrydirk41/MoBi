r"""Run the paper-vs-my-model comparison end to end, from the Python side.

sb_paper_compare is a MATLAB function - it runs INSIDE the MATLAB Engine, not at a
PowerShell/conda prompt. This driver starts the engine (as the other run_* scripts
do), calls sb_paper_compare to produce the three full-Vpop arm CSVs, then prints the
three-column model-vs-model table with paper_compare.report().

    set ANTHROPIC_API_KEY not needed - this is pure simulation, no LLM.

    python -m examples.run_paper_compare ^
        --sbproj "..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj" ^
        --vpop   "..\RA-QSP-Model\Vpop1.xlsx" ^
        --out    examples\matlab\out

This is ~900 simulations (3 arms x 300 patients, full population - no subsampling),
so it takes a while. That is deliberate: matching the paper's n is the whole point.
"""

from __future__ import annotations

import argparse
import os

from pkpd_agent.engines.simbiology import SimBiologyEngine
from examples import paper_compare


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sbproj", required=True)
    ap.add_argument("--vpop", required=True)
    ap.add_argument("--out", default="paper_compare_out",
                    help="folder for the three arm CSVs (created if missing)")
    args = ap.parse_args()

    out = os.path.abspath(args.out)
    os.makedirs(out, exist_ok=True)

    sb = SimBiologyEngine()
    try:
        print("== starting MATLAB/SimBiology engine ==", flush=True)
        sb.start()
        print("== running the three full-Vpop arms (this is the slow part) ==",
              flush=True)
        # sb_paper_compare loads the project itself and calls sb_run_vpop x3.
        sb.eng.sb_paper_compare(os.path.abspath(args.sbproj),
                                os.path.abspath(args.vpop), out, nargout=0)
    finally:
        sb.stop()

    print("\n== MODEL-vs-MODEL COMPARISON ==\n", flush=True)
    paper_compare.report(out)


if __name__ == "__main__":
    main()
