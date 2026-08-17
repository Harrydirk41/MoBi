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

    def run_vpop(self, vpop_xlsx: str, dose: str = "", baseline_day: float = 200.0,
                 readout_day: float = 284.0, limit: int = 0) -> dict[str, Any]:
        """Run the whole virtual population (an .xlsx of patient parameter sets)
        under an optional named dose, and return each patient's DAS28-CRP at the
        treatment-start BASELINE day and at the READOUT day. The clinical response
        (ACR = % DAS28 improvement from the patient's own baseline) is computed by
        the caller. ``limit`` runs only the first N patients. MATLAB reads the
        Excel and writes the results CSV; only paths cross the boundary."""
        import io
        out = _tmp(".csv")
        so = io.StringIO()
        try:
            self.eng.sb_run_vpop(vpop_xlsx, dose or "", float(baseline_day),
                                 float(readout_day), out, float(limit),
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
