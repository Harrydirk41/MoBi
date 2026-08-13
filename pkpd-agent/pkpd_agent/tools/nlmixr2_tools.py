"""Real nlmixr2 backend (NLME population PK) via an Rscript subprocess.

nlmixr2 does the thing pkfit cannot: true nonlinear mixed-effects estimation
with random effects (SAEM / FOCEi), no NONMEM required. It is cross-platform
(needs R + a C/C++ compiler). This adapter:

  * builds a NONMEM-style CSV from the session dataset,
  * runs pkpd_agent/engines/r_workers/nlmixr2_fit.R,
  * parses the JSON result.

If R / nlmixr2 is not available it returns a clear error (or a synthetic result
in mock mode) rather than crashing the loop.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile
from typing import Any

from .registry import Tool, ToolRegistry, ToolResult

_WORKER = os.path.join(os.path.dirname(os.path.dirname(__file__)),
                       "engines", "r_workers", "nlmixr2_fit.R")


class Nlmixr2Engine:
    def __init__(self, config) -> None:
        self.config = config

    def _rscript(self) -> str | None:
        cand = self.config.rscript_path
        if cand and (os.path.exists(cand) or shutil.which(cand)):
            return cand
        return shutil.which("Rscript")

    def available(self) -> tuple[bool, str]:
        rs = self._rscript()
        if not rs:
            return False, "Rscript not found (set config.rscript_path)"
        try:
            p = subprocess.run(
                [rs, "-e", "cat(requireNamespace('nlmixr2', quietly=TRUE))"],
                capture_output=True, text=True, encoding="utf-8",
                errors="replace", timeout=60)
            if "TRUE" in p.stdout:
                return True, rs
            return False, "R found but nlmixr2 package not installed"
        except Exception as exc:  # noqa: BLE001
            return False, f"{type(exc).__name__}: {exc}"

    # -- data -> NONMEM-style CSV --------------------------------------- #
    @staticmethod
    def _write_csv(data, path: str) -> None:
        import numpy as np  # local import; only needed on the real path
        subjects = np.unique(data.subject)
        rows = ["ID,TIME,DV,AMT,EVID,CMT"]
        for sid in subjects:
            mask = data.subject == sid
            dose = float(data.dose[mask][0])
            rows.append(f"{int(sid)},0,0,{dose},1,1")  # dosing event into depot
            for t, dv in zip(data.time[mask], data.dv[mask]):
                rows.append(f"{int(sid)},{float(t)},{float(dv)},0,0,2")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("\n".join(rows))

    def fit(self, data, model: str) -> dict[str, Any]:
        if self.config.mock:
            return {"model": model, "ofv": 340.0, "aic": 350.0,
                    "parameter_estimates": {"tcl": 1.2, "tv": 3.8, "tka": 0.1},
                    "relative_standard_errors": {"tcl": 0.06, "tv": 0.07},
                    "minimization_successful": True, "source": "nlmixr2-mock"}
        ok, info = self.available()
        if not ok:
            raise RuntimeError(f"nlmixr2 backend unavailable: {info} "
                               "(see SETUP_WINDOWS.md)")
        rs = info
        with tempfile.TemporaryDirectory() as d:
            csv_path = os.path.join(d, "data.csv")
            out_path = os.path.join(d, "out.json")
            self._write_csv(data, csv_path)
            proc = subprocess.run(
                [rs, _WORKER, csv_path, out_path, model, self.config.nlmixr2_est],
                capture_output=True, text=True, encoding="utf-8",
                errors="replace", timeout=1200)
            if not os.path.exists(out_path):
                raise RuntimeError(f"nlmixr2 worker produced no output: "
                                   f"{proc.stderr[:400]}")
            with open(out_path, encoding="utf-8") as fh:
                result = json.load(fh)
        if not result.get("ok", False):
            raise RuntimeError(result.get("message", "nlmixr2 fit failed"))
        return result


def register_nlmixr2_tools(registry: ToolRegistry, config) -> None:
    engine = Nlmixr2Engine(config)

    def fit(args: dict[str, Any], session) -> ToolResult:
        data = session.get("dataset")
        if data is None:
            return ToolResult.error("no dataset loaded - call pkfit_load_data first")
        res = engine.fit(data, args.get("model", "1cpt_oral"))
        fid = f"nlmixr2::{res.get('model')}"
        session.put(fid, res)
        content = {k: v for k, v in res.items()}
        content["fit_id"] = fid
        return ToolResult.success(f"nlmixr2 fit complete: {fid}", **content)

    registry.register(Tool(
        name="nlmixr2_fit",
        description=(
            "Fit a TRUE population NLME model (random effects on CL and V) to "
            "the loaded dataset with nlmixr2 (FOCEi/SAEM). This is the real "
            "mixed-effects estimation that pkfit's naive-pooling cannot do. "
            "model is '1cpt_oral' or '2cpt_oral'. Returns OFV, AIC, fixed-effect "
            "estimates and their RSEs. ACT. Requires R + nlmixr2 (see "
            "SETUP_WINDOWS.md); errors clearly if unavailable."
        ),
        input_schema={
            "type": "object",
            "properties": {
                "model": {"type": "string", "enum": ["1cpt_oral", "2cpt_oral"]},
            },
            "required": ["model"],
        },
        handler=fit,
        phase="act",
    ))
