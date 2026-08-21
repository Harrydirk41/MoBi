r"""Convert a standard SBML export to network.json - no MATLAB needed.

SimBiology (File > Export > SBML, or `sbmlexport(model,'model.xml')`) and every other QSP
tool can emit SBML; this parses it into the same network.json shape the benchmarks consume,
so the whole Stage-1 suite can run without the bespoke MATLAB dumper.

    python -m examples.sbml_to_network_json --sbml model.xml --out network.json
    python -m examples.run_llm_qsp_all --network network.json --infer
"""

from __future__ import annotations

import argparse
import json

from pkpd_agent.engines.sbml_import import sbml_to_network


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sbml", required=True)
    ap.add_argument("--out", default="network.json")
    args = ap.parse_args()

    data = sbml_to_network(args.sbml)
    with open(args.out, "w", encoding="utf-8") as fh:
        json.dump(data, fh)
    print(f"wrote {args.out}: {len(data['species'])} species, {len(data['reactions'])} "
          f"reactions, {len(data['rules'])} rules, {len(data['parameters'])} parameters")
    print("next:  python -m examples.run_llm_qsp_all --network "
          f"{args.out} --infer")


if __name__ == "__main__":
    main()
