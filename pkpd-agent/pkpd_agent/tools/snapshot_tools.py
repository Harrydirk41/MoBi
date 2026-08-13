"""Tool: extract agent-ready data from an OSP snapshot JSON.

Turns any PK-Sim/MoBi snapshot into observed data, compound parameters, the
modeling choices, and an NCA summary - the data/analysis half of the PBPK
pipeline, with no OSP install required. (Simulation still needs a .pkml +
ospsuite; this tool does not run the model.)
"""

from __future__ import annotations

from typing import Any

from ..engines.osp_snapshot import OSPSnapshot
from .registry import Tool, ToolRegistry, ToolResult


def register_snapshot_tools(registry: ToolRegistry, config) -> None:

    def extract(args: dict[str, Any], session) -> ToolResult:
        snap = OSPSnapshot.from_file(args["path"])
        session.put("osp_snapshot", snap)

        summary = snap.summary()
        nca = snap.nca_table()
        choices = snap.modeling_choices()
        session.put("snapshot_observed", snap.observed_profiles())
        session.put("snapshot_nca", nca)

        written = {}
        if args.get("outdir"):
            written = snap.write_csvs(args["outdir"])

        return ToolResult.success(
            f"extracted {summary.get('compound')}: "
            f"{summary.get('n_observed_datasets')} observed datasets, "
            f"{len(snap.compound_parameters())} compound parameters",
            **summary,
            nca_summary=nca,
            modeling_choices=choices,
            n_parameters=len(snap.compound_parameters()),
            files=written,
        )

    registry.register(Tool(
        name="snapshot_extract",
        description=(
            "Extract agent-ready data from an OSP (PK-Sim/MoBi) snapshot JSON "
            "file: the observed clinical concentration-time datasets, the "
            "compound parameter file (physchem/ADME/metabolism), the modeling "
            "choices (distribution method, enzymes/transporters, which inputs "
            "were measured vs fitted), and a per-study NCA summary "
            "(Cmax/Tmax/AUC/half-life). OBSERVE step. Pass outdir to also write "
            "observed + parameter CSVs. Does NOT run a simulation - use "
            "osp_simulate on a .pkml for that."
        ),
        input_schema={
            "type": "object",
            "properties": {
                "path": {"type": "string", "description": "Path to the snapshot .json"},
                "outdir": {"type": "string",
                           "description": "Optional folder to write CSVs into"},
            },
            "required": ["path"],
        },
        handler=extract,
        phase="observe",
    ))
