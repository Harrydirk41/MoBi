r"""Generate the Stage-1 answer key: dump the RA model's full wiring to network.json.

Runs the MATLAB dumper (sb_network_json.m) via the engine so you never touch a raw
MATLAB prompt. Do this ONCE; then run the reconstruction benchmark against the result:

    python -m examples.dump_network ^
        --sbproj "..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj" ^
        --out network.json
    python -m examples.run_llm_qsp_all --network network.json --model ra

No API key needed here - it is pure model extraction, no LLM.
"""

from __future__ import annotations

import argparse

from pkpd_agent.engines.qsp_model import QSPModel
from pkpd_agent.engines.simbiology import SimBiologyEngine


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sbproj", required=True)
    ap.add_argument("--out", default="network.json")
    args = ap.parse_args()

    sb = SimBiologyEngine()
    try:
        print("== starting MATLAB/SimBiology engine ==", flush=True)
        sb.start()
        sb.load_project(args.sbproj)
        s = sb.network_json(args.out)
    finally:
        sb.stop()

    c = s.get("counts", {})
    print(f"wrote {args.out}: {c.get('species')} species, {c.get('reactions')} "
          f"reactions, {c.get('rules')} rules, {c.get('parameters')} parameters")
    model = QSPModel.inferred(args.out, "auto")
    print(f"derived {len(model.nodes)} nodes, {len(model.edges)} regulatory edges, "
          f"{len(model.readout_drivers)} readout drivers (the reconstruction answer key)")


if __name__ == "__main__":
    main()
