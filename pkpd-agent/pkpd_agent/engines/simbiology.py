"""Drive a SimBiology (MATLAB) QSP/PKPD model headlessly from Python.

The SimBiology analogue of OSPCli. A persistent MATLAB Engine session loads a
.sbproj once and simulates it many times. Because this MATLAB Engine build crashes
when marshaling ARRAYS back to Python, all bulk data crosses the boundary via
FILES (the .m helpers write JSON/CSV; Python reads them) - exactly the PKSim.CLI
pattern. Only scalars/strings go through the engine directly.

    from pkpd_agent.engines.simbiology import SimBiologyEngine
    with SimBiologyEngine() as sb:
        sb.load_project(r"C:\\...\\Vantage RA QSP Model v1.0.sbproj")
        info = sb.model_info()                 # names of species/params/doses/variants
        base = sb.simulate()                   # baseline
        tcz  = sb.simulate(dose="TCZ8mgkg_Q4W_IV_t200")
"""

from __future__ import annotations

import csv
import json
import os
import tempfile
from dataclasses import dataclass
from typing import Any

# where the sb_*.m helpers live (examples/matlab, alongside this package)
_DEFAULT_MATLAB_DIR = os.path.abspath(os.path.join(
    os.path.dirname(__file__), "..", "..", "examples", "matlab"))


def _tmp(suffix: str) -> str:
    fd, path = tempfile.mkstemp(suffix=suffix, prefix="sb_")
    os.close(fd)
    return path


def _read_csv(path: str) -> dict[str, Any]:
    """A results CSV (header row of names, then numeric rows) -> {names, time,
    columns:{name:[floats]}}."""
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.reader(fh))
    if not rows:
        return {"names": [], "time": [], "columns": {}}
    header = rows[0]
    cols: dict[str, list] = {h: [] for h in header}
    for r in rows[1:]:
        for h, v in zip(header, r):
            try:
                cols[h].append(float(v))
            except (ValueError, TypeError):
                cols[h].append(None)
    time_key = header[0] if header else "Time"
    return {"names": header, "time": cols.get(time_key, []), "columns": cols}


@dataclass
class SimBiologyEngine:
    matlab_dir: str = _DEFAULT_MATLAB_DIR
    eng: Any = None

    # -- lifecycle ------------------------------------------------------- #
    def start(self) -> "SimBiologyEngine":
        import matlab.engine
        self.eng = matlab.engine.start_matlab()
        self.eng.addpath(self.matlab_dir, nargout=0)
        return self

    def stop(self) -> None:
        if self.eng is not None:
            try:
                self.eng.quit()
            except Exception:                          # noqa: BLE001
                pass
            self.eng = None

    def __enter__(self) -> "SimBiologyEngine":
        return self.start()

    def __exit__(self, *exc) -> None:
        self.stop()

    # -- model ----------------------------------------------------------- #
    def load_project(self, path: str) -> None:
        """Load a .sbproj; the model is kept in MATLAB as 'sbmodel'."""
        self.eng.sb_load(path, nargout=0)

    def has_simbiology(self) -> bool:
        try:
            return bool(self.eng.license("test", "SimBiology", nargout=1))
        except Exception:                              # noqa: BLE001
            return False

    def model_info(self) -> dict[str, Any]:
        """Names/counts of species, parameters, reactions, doses, variants."""
        out = _tmp(".json")
        try:
            self.eng.sb_model_json(out, nargout=0)
            with open(out, encoding="utf-8") as fh:
                return json.load(fh)
        finally:
            _quiet_rm(out)

    def simulate(self, dose: str = "", variant: str = "",
                 stop_time: float = 0.0) -> dict[str, Any]:
        """Simulate the loaded model, optionally with a named dose and/or variant.
        Returns {names, time, columns:{state:[...]}} read from a results CSV."""
        out = _tmp(".csv")
        try:
            n = self.eng.sb_simulate_csv(dose or "", variant or "",
                                         float(stop_time), out, nargout=1)
            res = _read_csv(out)
            res["n_columns"] = int(n)
            return res
        finally:
            _quiet_rm(out)

    def add_drug(self, target: str, efficacy: float,
                 start_day: float = 200.0) -> None:
        """Add a designed anti-cytokine drug to the loaded model: from ``start_day``
        it suppresses the ``target`` disease-driver parameter to (1-efficacy) of its
        baseline. Reload the project (load_project) to reset before a new design."""
        self.eng.sb_add_drug(target, float(efficacy), float(start_day), nargout=0)

    def sample_vpop(self, param_spec: str, n_samples: int = 60,
                    baseline_day: float = 200.0, seed: int = 1) -> dict[str, Any]:
        """Generate a virtual population by sampling disease-driver parameters and
        simulating each candidate to its untreated disease baseline. ``param_spec``
        is ';'-joined 'name,lo,hi,scale' entries (scale 'lin'/'log'). Returns the
        per-candidate parameter draws and baseline DAS28-CRP (column 'DAS28_base')."""
        import io
        out = _tmp(".csv")
        so = io.StringIO()
        try:
            self.eng.sb_sample_vpop(param_spec, float(n_samples), float(baseline_day),
                                    out, float(seed), nargout=0, stdout=so, stderr=so)
            res = _read_csv(out)
            res["matlab_log"] = so.getvalue()
            return res
        finally:
            _quiet_rm(out)

    def run_vpop(self, vpop_xlsx: str, dose: str = "", stop_time: float = 700.0,
                 baseline_day: float = 200.0, readout_day: float = 284.0,
                 limit: int = 0, param_overrides: str = "") -> dict[str, Any]:
        """Run the whole virtual population (an .xlsx of patient parameter sets)
        under one or more named doses, and return each patient's MODEL-COMPUTED
        clinical-response flags. The Vantage RA model encodes the trial as events:
        it captures DAS28_BL (day 199), sets ACR20/50/70/Remission (day 284,
        first-line week-12 readout) and - for MTX_NonResp patients - sets
        MTX_NonResp_TCZ_ACR20/50/70/Rem (day 600, the second-line TCZ readout). So
        the response is read from the model, not recomputed.

        ``dose`` may join several dose names with ';' (a SimBiology dose array is
        applied), e.g. 'MTX_15mg_Q1W_SC_t200;TCZ8mgkg_Q4W_IV_t200'. ``stop_time``
        forces the sim end (needs >=601 days to capture the second-line flags).
        ``limit`` runs a representative EVENLY-SPACED subsample of N patients (the
        Vpop rows are ordered by severity, so the first N would be the sickest
        slice). MATLAB reads the Excel and writes
        the results CSV; only paths cross the boundary. Returned columns include
        DAS28_BL/base/read/end and the flag columns ACR20, ACR50, ACR70, Rem,
        MTX_NonResp, TCZ_ACR20, TCZ_ACR50, TCZ_ACR70, TCZ_Rem."""
        import io
        out = _tmp(".csv")
        so = io.StringIO()
        try:
            self.eng.sb_run_vpop(vpop_xlsx, dose or "", float(stop_time),
                                 float(baseline_day), float(readout_day), out,
                                 float(limit), param_overrides or "",
                                 nargout=0, stdout=so, stderr=so)
            res = _read_csv(out)
            res["matlab_log"] = so.getvalue()
            return res
        finally:
            _quiet_rm(out)


def _quiet_rm(path: str) -> None:
    try:
        os.remove(path)
    except OSError:
        pass
