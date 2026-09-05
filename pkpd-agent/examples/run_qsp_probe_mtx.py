r"""Probe HOW MTX's pharmacodynamic effect is wired in the paper clinical model, so fit_clinical can
target the real free POTENCY constant instead of guessing by name suffix. Prints:
  - sb_drug_mechanism('MTX'): the reactions/rate-laws/RULES that reference MTX,
  - every MTX-named parameter with its value and whether it is CONSTANT (rule-driven outputs are not
    fittable; PK disposition constants must not be touched by an efficacy fit).

    python -m examples.run_qsp_probe_mtx --sbproj "..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj" \
        --modeldir "..\RA-QSP-Model"
"""

from __future__ import annotations

import argparse
import os


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sbproj", required=True)
    ap.add_argument("--modeldir", required=True)
    ap.add_argument("--drug", default="MTX")
    args = ap.parse_args()

    from pkpd_agent.engines.simbiology import SimBiologyEngine
    sb = SimBiologyEngine()
    sb.start()
    try:
        sb.eng.addpath(os.path.abspath(args.modeldir), nargout=0)
        sb.load_project(os.path.abspath(args.sbproj))
        print(f"\n== sb_drug_mechanism('{args.drug}') ==", flush=True)
        sb.eng.sb_drug_mechanism(args.drug, nargout=0)

        print(f"\n== all parameters matching /{args.drug}/ (name | value | constant) ==", flush=True)
        params = sb.list_parameters()["parameters"]
        for p in params:
            if args.drug.lower() in p["name"].lower():
                print(f"  {p['name']:28} {p.get('value')!s:>14}   constant={p.get('constant')}")
    finally:
        sb.stop()


if __name__ == "__main__":
    main()
