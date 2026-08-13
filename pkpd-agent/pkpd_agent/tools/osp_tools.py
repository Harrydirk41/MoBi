"""OSP engine adapter + tools (mechanistic PBPK / QSP via MoBi / PK-Sim).

OSP is a .NET/Windows stack driven headlessly through one of:
  * the R package ``ospsuite`` (Rscript subprocess), or
  * ``MoBi.CLI`` operating on snapshot JSON.

The model is exchanged as a *snapshot* (JSON) - the machine-readable
representation an LLM can read and write. In ``mock`` mode the handlers return
synthetic simulation results (including mass-balance and dimensional info the
verification gates check). The real subprocess construction is written out and
guarded.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import tempfile
from typing import Any

from .registry import Tool, ToolRegistry, ToolResult

_OSP_WORKER = os.path.join(os.path.dirname(os.path.dirname(__file__)),
                           "engines", "r_workers", "osp_sim.R")


class OSPEngine:
    def __init__(self, config) -> None:
        self.config = config

    @staticmethod
    def _seed(*parts: str) -> int:
        h = hashlib.sha256("|".join(parts).encode()).hexdigest()
        return int(h[:8], 16)

    def load_snapshot(self, path: str) -> dict[str, Any]:
        if self.config.mock:
            return {
                "snapshot_id": f"snap::{path}",
                "molecules": ["Drug", "Target", "Complex"],
                "compartments": ["Plasma", "Liver", "Kidney", "Muscle"],
                "reactions": ["binding", "internalization", "target_turnover"],
                "source": "mock",
            }
        with open(path, "r", encoding="utf-8") as fh:
            snap = json.load(fh)
        return {
            "snapshot_id": f"snap::{path}",
            "molecules": [m.get("Name") for m in snap.get("Molecules", [])],
            "compartments": [c.get("Name") for c in snap.get("SpatialStructures", [])],
            "source": "osp",
        }

    def simulate(self, snapshot_id: str, output: str) -> dict[str, Any]:
        if self.config.mock:
            s = self._seed("sim", snapshot_id, output)
            # A plausible PK-ish curve summary plus the invariants the gates check.
            cmax = round(4.0 + (s % 60) / 10, 3)
            return {
                "snapshot_id": snapshot_id,
                "output": output,
                "t_max_h": round(1.0 + (s % 30) / 10, 2),
                "c_max": cmax,
                "auc": round(cmax * (8 + s % 12), 2),
                "all_values_finite": True,
                "min_concentration": 0.0,       # non-negative -> physical
                "mass_balance_residual": round((s % 5) / 1e6, 9),
                "source": "mock",
            }
        return self._simulate_real(snapshot_id, output)

    def _simulate_real(self, snapshot_id: str, output: str) -> dict[str, Any]:
        """Drive OSP headlessly via the ospsuite R worker (osp_sim.R).

        snapshot_id encodes the path to a .pkml simulation. Requires an Rscript
        with ospsuite installed (see SETUP_WINDOWS.md)."""
        path = snapshot_id.replace("snap::", "", 1)
        rs = self.config.rscript_path
        if not (rs and (os.path.exists(rs) or shutil.which(rs))):
            raise RuntimeError("Rscript not found for OSP (set config.rscript_path)")
        with tempfile.TemporaryDirectory() as d:
            out_path = os.path.join(d, "out.json")
            proc = subprocess.run([rs, _OSP_WORKER, path, out_path, output or ""],
                                  capture_output=True, text=True, encoding="utf-8",
                                  errors="replace", timeout=1200)
            if not os.path.exists(out_path):
                raise RuntimeError(f"OSP worker produced no output: {proc.stderr[:400]}")
            with open(out_path, encoding="utf-8") as fh:
                result = json.load(fh)
        if not result.get("ok", False):
            raise RuntimeError(result.get("message", "OSP simulation failed"))
        result["snapshot_id"] = snapshot_id
        return result

    def set_parameter(self, snapshot_id: str, path: str, value: float) -> dict[str, Any]:
        if self.config.mock:
            return {"snapshot_id": snapshot_id, "path": path, "value": value,
                    "applied": True, "source": "mock"}
        raise NotImplementedError("real snapshot editing not wired in this skeleton")


# --------------------------------------------------------------------------- #
# Tool registration
# --------------------------------------------------------------------------- #

def register_osp_tools(registry: ToolRegistry, config) -> None:
    engine = OSPEngine(config)

    def load_snapshot(args: dict[str, Any], session) -> ToolResult:
        info = engine.load_snapshot(args["path"])
        session.put(info["snapshot_id"], {"kind": "osp_snapshot", **info})
        return ToolResult.success(f"loaded {info['snapshot_id']}", **info)

    def set_parameter(args: dict[str, Any], session) -> ToolResult:
        res = engine.set_parameter(args["snapshot_id"], args["path"], float(args["value"]))
        return ToolResult.success("parameter set", **res)

    def simulate(args: dict[str, Any], session) -> ToolResult:
        res = engine.simulate(args["snapshot_id"], args.get("output", "Plasma|Drug"))
        session.put(f"sim::{args['snapshot_id']}", res)
        return ToolResult.success("simulation complete", **res)

    registry.register(Tool(
        name="osp_load_snapshot",
        description=(
            "Load an OSP (MoBi/PK-Sim) model snapshot (JSON) into the session. "
            "OBSERVE step. Returns molecules, compartments and reactions - the "
            "structure an LLM can reason over before editing."
        ),
        input_schema={
            "type": "object",
            "properties": {"path": {"type": "string"}},
            "required": ["path"],
        },
        handler=load_snapshot,
        phase="observe",
    ))

    registry.register(Tool(
        name="osp_set_parameter",
        description=(
            "Set a parameter (by container path) in a loaded snapshot, e.g. a "
            "clearance, partition coefficient, or dose. An ACT step that "
            "mutates the mechanistic model."
        ),
        input_schema={
            "type": "object",
            "properties": {
                "snapshot_id": {"type": "string"},
                "path": {"type": "string", "description": "Container path of the parameter."},
                "value": {"type": "number"},
            },
            "required": ["snapshot_id", "path", "value"],
        },
        handler=set_parameter,
        phase="act",
        destructive=True,
    ))

    registry.register(Tool(
        name="osp_simulate",
        description=(
            "Integrate the ODE system for a loaded snapshot and summarize an "
            "output curve (Tmax, Cmax, AUC) plus physical-sanity invariants "
            "(all values finite, minimum concentration, mass-balance residual). "
            "An ACT/EVALUATE step - the verification gates read the invariants."
        ),
        input_schema={
            "type": "object",
            "properties": {
                "snapshot_id": {"type": "string"},
                "output": {"type": "string", "description": "Observer path, e.g. 'Plasma|Drug'."},
            },
            "required": ["snapshot_id"],
        },
        handler=simulate,
        phase="act",
    ))
