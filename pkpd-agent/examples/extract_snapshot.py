"""Extract agent-ready data from an OSP snapshot JSON.

    python -m examples.extract_snapshot <snapshot.json> [outdir]

Writes <Compound>_observed.csv and <Compound>_compound_params.csv into
`outdir` (default: a "<snapshot>_extracted" folder next to the JSON) and
prints the summary, the modeling choices (the PBPK "answer key"), and a
per-study NCA table.

Example (Windows):
    python -m examples.extract_snapshot ^
        "C:\\Users\\harry\\PyCharmMiscProject\\MoBi\\temp_analysis\\Clarithromycin-Model.json"
"""

import os
import sys

from pkpd_agent.config import AgentConfig
from pkpd_agent.state import ModelingSession
from pkpd_agent.tools import build_default_registry


def main() -> None:
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    path = sys.argv[1]
    if not os.path.exists(path):
        print(f"file not found: {path}")
        sys.exit(1)

    outdir = sys.argv[2] if len(sys.argv) > 2 else (
        os.path.splitext(path)[0] + "_extracted")

    reg = build_default_registry(AgentConfig())
    session = ModelingSession(goal="extract snapshot")
    r = reg.dispatch("snapshot_extract", {"path": path, "outdir": outdir}, session)

    if not r.ok:
        print("extraction failed:", r.message)
        sys.exit(1)

    d = r.data
    print("=" * 70)
    print(r.message)
    print("=" * 70)
    print(f"compound              : {d['compound']}")
    print(f"observed datasets     : {d['n_observed_datasets']}")
    print(f"compound parameters   : {d['n_parameters']}")
    print(f"simulations/protocols : {d['n_simulations']} / {d['n_protocols']}")

    mc = d["modeling_choices"]
    print("\n--- modeling choices (the PBPK 'answer key') ---")
    print("  distribution/permeability:")
    for m in mc["calculation_methods"]:
        print(f"     {m}")
    procs = [p["molecule"] for p in mc["processes"] if p["molecule"]]
    print(f"  processes (enzymes/transporters): {procs}")
    print(f"  measured vs fitted: {mc['fit_vs_measured']}")
    print(f"  model type: {mc['model_type']}")

    print("\n--- NCA per study (from observed data) ---")
    print(f"  {'study':34}{'Cmax':>8}{'Tmax':>7}{'AUC':>9}{'t½(h)':>7}")
    for row in d["nca_summary"]:
        th = row.get("t_half_h")
        print(f"  {row['study'][:34]:34}{row.get('c_max_mg_L',0):>8.3f}"
              f"{row.get('t_max_h',0):>7.2f}{row.get('auc_mg_h_L',0):>9.2f}"
              f"{(th if th else 0):>7.2f}")

    print("\nsaved files:")
    for k, p in d["files"].items():
        print(f"  {k}: {p}")


if __name__ == "__main__":
    main()
