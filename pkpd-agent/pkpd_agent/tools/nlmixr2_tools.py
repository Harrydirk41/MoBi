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
    def _write_csv(data, path: str, covariate: dict | None = None) -> None:
        import numpy as np  # local import; only needed on the real path
        subjects = np.unique(data.subject)
        header = "ID,TIME,DV,AMT,EVID,CMT"
        cov_vals = None
        ref = 1.0
        if covariate:
            cov_name = covariate["cov"]
            ref = float(covariate.get("ref", 1.0))
            cov_vals = data.covariates[cov_name]
            header += ",WTR"
        rows = [header]
        for sid in subjects:
            mask = data.subject == sid
            dose = float(data.dose[mask][0])
            wtr = f",{float(cov_vals[mask][0]) / ref}" if cov_vals is not None else ""
            rows.append(f"{int(sid)},0,0,{dose},1,1{wtr}")  # dosing event
            for t, dv in zip(data.time[mask], data.dv[mask]):
                rows.append(f"{int(sid)},{float(t)},{float(dv)},0,0,2{wtr}")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("\n".join(rows))

    def fit(self, data, model: str, covariate: dict | None = None) -> dict[str, Any]:
        covparam = covariate["param"] if covariate else "none"
        if self.config.mock:
            out = {"model": model, "covariate": None if not covariate else
                   f"WT_on_{covparam}", "ofv": 340.0, "aic": 350.0,
                   "parameter_estimates": {"tcl": 1.27, "tv": 3.75, "tka": 0.10,
                                           "prop.err": 0.15},
                   "relative_standard_errors": {"tcl": 0.06, "tv": 0.02},
                   "iiv_cv_percent": {"eta.cl": 21.0, "eta.v": 19.0},
                   "shrinkage_percent": [8.0, 11.0],
                   "pct_observations_within_90_pi": 90.5,
                   "minimization_successful": True, "source": "nlmixr2-mock"}
            if covariate:
                out["parameter_estimates"][f"cov_{covparam.lower()}"] = 0.78
            return out
        ok, info = self.available()
        if not ok:
            raise RuntimeError(f"nlmixr2 backend unavailable: {info} "
                               "(see SETUP_WINDOWS.md)")
        rs = info
        with tempfile.TemporaryDirectory() as d:
            csv_path = os.path.join(d, "data.csv")
            out_path = os.path.join(d, "out.json")
            self._write_csv(data, csv_path, covariate)
            proc = subprocess.run(
                [rs, _WORKER, csv_path, out_path, model,
                 self.config.nlmixr2_est, covparam],
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
        covariate = None
        if args.get("covariate_param"):
            covariate = {"param": args["covariate_param"],
                         "cov": args.get("covariate", "WT"),
                         "ref": float(args.get("covariate_ref", 70.0))}
        res = engine.fit(data, args.get("model", "1cpt_oral"), covariate)
        cov = res.get("covariate")
        fid = f"nlmixr2::{res.get('model')}" + (f"+{cov}" if cov else "")
        session.put(fid, res)
        content = dict(res)
        content["fit_id"] = fid
        return ToolResult.success(f"nlmixr2 fit complete: {fid}", **content)

    def vpc(args: dict[str, Any], session) -> ToolResult:
        fit_obj = session.get(args["fit_id"])
        if fit_obj is None:
            return ToolResult.error(f"unknown fit_id {args['fit_id']!r}")
        cov = fit_obj.get("pct_observations_within_90_pi")
        if cov is None:
            return ToolResult.error("no VPC coverage was computed for this fit "
                                    "(nlmixr2 VPC is available for 1cpt_oral).")
        return ToolResult.success(
            "VPC (from the nlmixr2 fit)",
            pct_observations_within_90_pi=cov, model=fit_obj.get("model"),
            source="nlmixr2-vpc")

    registry.register(Tool(
        name="nlmixr2_fit",
        description=(
            "Fit a TRUE population NLME model (random effects on CL and V) to "
            "the loaded dataset with nlmixr2 (FOCEi/SAEM) - the real "
            "mixed-effects estimation pkfit's naive-pooling cannot do. model is "
            "'1cpt_oral' or '2cpt_oral'. Optionally add an allometric covariate "
            "estimated INSIDE the NLME fit via covariate_param ('CL' or 'V'), "
            "covariate ('WT'), covariate_ref. Returns OFV/AIC, fixed-effect "
            "estimates + RSEs, the covariate coefficient, IIV as CV%, "
            "shrinkage, and a model-based VPC coverage. Prefer this over pkfit "
            "for the final model and for estimating covariate exponents "
            "(naive-pooling biases them). ACT. Requires R + nlmixr2."
        ),
        input_schema={
            "type": "object",
            "properties": {
                "model": {"type": "string", "enum": ["1cpt_oral", "2cpt_oral"]},
                "covariate_param": {"type": "string", "enum": ["CL", "V"]},
                "covariate": {"type": "string", "description": "e.g. 'WT'"},
                "covariate_ref": {"type": "number"},
            },
            "required": ["model"],
        },
        handler=fit,
        phase="act",
    ))

    registry.register(Tool(
        name="nlmixr2_vpc",
        description=(
            "Report the model-based VPC coverage (fraction of observations "
            "within the 90% prediction interval) for a fitted nlmixr2 model, "
            "by fit_id. Unlike the pkfit VPC, this one includes the estimated "
            "IIV, so it attributes variability correctly. EVALUATE."
        ),
        input_schema={
            "type": "object",
            "properties": {"fit_id": {"type": "string"}},
            "required": ["fit_id"],
        },
        handler=vpc,
        phase="evaluate",
    ))
