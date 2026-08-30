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
import io
import json
import os
import sys
import tempfile
from dataclasses import dataclass
from typing import Any


class _TeeStringIO(io.StringIO):
    """An io.StringIO that ALSO echoes each write to the real console. matlab.engine
    requires stdout to be an io.StringIO (a subclass passes the isinstance check), so this
    lets a long MATLAB call show its fprintf progress live instead of only on return."""

    def write(self, s):
        try:
            sys.__stdout__.write(s)
            sys.__stdout__.flush()
        except Exception:
            pass
        return super().write(s)

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

    def network_json(self, out_path: str) -> dict[str, Any]:
        """Dump the full network structure (species, reactions with rate laws, rules,
        parameters) to ``out_path`` and return it. The Stage-1 reconstruction answer
        key. Writes to a real path (kept), not a temp file."""
        self.eng.sb_network_json(os.path.abspath(out_path), nargout=0)
        with open(out_path, encoding="utf-8") as fh:
            return json.load(fh)

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
                 limit: int = 0, param_overrides: str = "",
                 states: "list | None" = None) -> dict[str, Any]:
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
            state_spec = ";".join(states) if states else ""
            self.eng.sb_run_vpop(vpop_xlsx, dose or "", float(stop_time),
                                 float(baseline_day), float(readout_day), out,
                                 float(limit), param_overrides or "", state_spec,
                                 nargout=0, stdout=so, stderr=so)
            res = _read_csv(out)
            res["matlab_log"] = so.getvalue()
            return res
        finally:
            _quiet_rm(out)

    def knockout_readout(self, param_names: "list | None", dose: str,
                         readout_day: float, readout_state: str = "DAS28_CRP",
                         stop_time: float = 700.0) -> float:
        """Ablate a regulatory edge and read the clinical readout off the REAL model. Freezes
        the named rule-parameters (deactivates their assigning rules, holds them constant),
        severing their source-dependence, then simulates under ``dose`` and returns
        ``readout_state`` (default DAS28_CRP) at ``readout_day``. Empty ``param_names`` = the
        baseline run (nothing knocked out). The model is restored afterwards, so this can be
        called repeatedly against one loaded project. This is the knob behind
        ``llm_topology.functional_weights``: |readout(knockout) - readout(baseline)| weights an
        edge by its clinical impact."""
        import io
        names = list(param_names or [])
        pcsv = ""
        if names:
            pcsv = _tmp(".txt")
            with open(pcsv, "w", encoding="utf-8") as fh:
                fh.write("\n".join(names))
        so = io.StringIO()
        try:
            val = self.eng.sb_knockout_readout(pcsv, dose or "", float(stop_time),
                                               float(readout_day), readout_state,
                                               nargout=1, stdout=so, stderr=so)
            return float(val)
        finally:
            if pcsv:
                _quiet_rm(pcsv)

    def cohort_multi_arm(self, param_spec: str, arms_spec: str, baseline_day: float,
                         readout_day: float, n_samples: int, seed: int = 1,
                         states: "list | None" = None, n_extra: int = 2,
                         stream: bool = False, seed_csv: str = "") -> dict[str, Any]:
        """Sample a virtual cohort and record each candidate's untreated baseline
        severity AND its primary response flag under several therapy arms - the rich
        cohort the multi-anchor Vpop selection needs. ``arms_spec`` is arms joined by
        ';;', each 'label:dose1,dose2'. Returns {columns:{sample, <params>, sev_base,
        <arm labels>}, matlab_log}. Each candidate is simulated once per arm in MATLAB
        (no bridge in the loop); the weight optimization afterwards is cheap."""
        import time
        out = _tmp(".csv")
        prog = out + ".prog"
        so = io.StringIO()
        args = (param_spec, arms_spec or "", float(baseline_day), float(readout_day),
                float(n_samples), float(seed), ";".join(states) if states else "", out,
                float(n_extra), seed_csv or "")
        try:
            if stream:
                # matlab.engine buffers stdout until the call returns, so live progress
                # can't come through it; instead run the call in the BACKGROUND and poll
                # the .prog file sb_cohort writes as it goes.
                fut = self.eng.sb_cohort(*args, nargout=0, stdout=so, stderr=so,
                                         background=True)
                last = ""
                while not fut.done():
                    try:
                        with open(prog, encoding="utf-8") as fh:
                            cur = fh.read().strip()
                        if cur and cur != last:
                            print(f"   cohort progress: {cur}", flush=True)
                            last = cur
                    except OSError:
                        pass
                    time.sleep(2)
                fut.result()
            else:
                self.eng.sb_cohort(*args, nargout=0, stdout=so, stderr=so)
            res = _read_csv(out)
            res["matlab_log"] = so.getvalue()
            return res
        finally:
            _quiet_rm(out)
            _quiet_rm(prog)

    def select_ga(self, cohort_columns: dict, anchor_spec: str,
                  pop_target: float = 0) -> dict[str, Any]:
        """Select a virtual population from an already-simulated cohort with SimBiology's
        native genetic algorithm (sb_select_ga.m) - the paper's Vpop method. It picks a
        SUBSET of real candidates whose aggregate matches the anchors, returning an actual
        population (not fractional weights). ``cohort_columns`` is the {name:[...]} dict
        from ``cohort_multi_arm``; ``anchor_spec`` is ';'-joined 'moment:COL:MEAN:SD' /
        'rate:COL:TARGET'. Returns {columns (the selected rows), n_selected, matlab_log}."""
        import io
        cohort = _tmp(".csv")
        out = _tmp(".csv")
        so = io.StringIO()
        try:
            names = list(cohort_columns.keys())
            nrows = max((len(v) for v in cohort_columns.values()), default=0)
            with open(cohort, "w", newline="", encoding="utf-8") as fh:
                w = csv.writer(fh)
                w.writerow(names)
                for i in range(nrows):
                    w.writerow(["" if i >= len(cohort_columns[nm]) or cohort_columns[nm][i]
                                is None else cohort_columns[nm][i] for nm in names])
            self.eng.sb_select_ga(cohort, anchor_spec, float(pop_target), out,
                                  nargout=0, stdout=so, stderr=so)
            res = _read_csv(out)
            res["n_selected"] = len(next(iter(res["columns"].values()), [])) \
                if res.get("columns") else 0
            res["matlab_log"] = so.getvalue()
            return res
        finally:
            _quiet_rm(cohort)
            _quiet_rm(out)

    def set_parameter(self, name: str, value: float) -> float:
        """Set a model parameter's value, returning its previous value. A general
        perturb/restore primitive (used by the calibration recover demo)."""
        return float(self.eng.sb_set_param(name, float(value), nargout=1))

    def perturb_response(self, species: str, high_value: float, readout_state: str,
                         readout_day: float, stop_time: float = 200.0,
                         decouple: "list | None" = None, decouple_value: float = 0.0) -> float:
        """Isolating perturbation: clamp ``species`` (a cytokine) at ``high_value`` and read
        ``readout_state`` (a cell) at ``readout_day`` - the single-cytokine experiment that
        pins one coupled regulator. ``decouple`` names other species to hold at
        ``decouple_value`` (~0), reproducing an in-vitro condition (the cell + only this
        cytokine, no network feedback) so the response isolates one regulator cleanly.
        Everything else stays at the model's current values; the model is restored afterwards."""
        return float(self.eng.sb_perturb_response(
            species, float(high_value), readout_state, float(readout_day), float(stop_time),
            ";".join(decouple or []), float(decouple_value), nargout=1))

    def list_parameters(self) -> dict[str, Any]:
        """Enumerate every parameter in the loaded model (name, value, constant) via
        sb_params.m - the full candidate set for driver selection, straight from the
        model rather than a pre-curated config. Returns {parameters:[{name,value,
        constant}], matlab_log}."""
        import io
        out = _tmp(".csv")
        so = io.StringIO()
        try:
            self.eng.sb_params(out, nargout=0, stdout=so, stderr=so)
            params = []
            try:
                with open(out, newline="", encoding="utf-8") as fh:
                    rdr = csv.reader(fh)
                    next(rdr, None)
                    for row in rdr:
                        if not row or not row[0].strip():
                            continue
                        val = None
                        const = None
                        if len(row) > 1:
                            try:
                                val = float(row[1])
                            except (ValueError, TypeError):
                                pass
                        if len(row) > 2:
                            try:
                                const = bool(int(row[2])) if int(row[2]) >= 0 else None
                            except (ValueError, TypeError):
                                pass
                        params.append({"name": row[0].strip(), "value": val,
                                       "constant": const})
            except OSError:
                pass
            return {"parameters": params, "matlab_log": so.getvalue()}
        finally:
            _quiet_rm(out)

    def gsa(self, param_spec: str, observable: str, readout_day: float,
            n_samples: int = 1000) -> dict[str, Any]:
        """Rank candidate parameters by global (Sobol) sensitivity on a readout, computed
        natively (sb_gsa.m) - the numerical half of choosing which parameters to vary.
        ``param_spec`` is ';'-joined 'name,lo,hi'; ``observable`` is the model readout.
        Returns {columns: {parameter, first_order, total_order}, matlab_log}, ranked by
        total order (most influential first)."""
        import io
        out = _tmp(".csv")
        so = io.StringIO()
        try:
            self.eng.sb_gsa(param_spec, observable, float(readout_day), float(n_samples),
                            out, nargout=0, stdout=so, stderr=so)
            names, fo, to = [], [], []
            try:
                with open(out, newline="", encoding="utf-8") as fh:
                    rdr = csv.reader(fh)
                    next(rdr, None)
                    for row in rdr:
                        if len(row) < 3:
                            continue
                        names.append(row[0].strip())
                        for dst, cell in ((fo, row[1]), (to, row[2])):
                            try:
                                dst.append(float(cell))
                            except (ValueError, TypeError):
                                dst.append(None)
            except OSError:
                pass
            return {"columns": {"parameter": names, "first_order": fo, "total_order": to},
                    "matlab_log": so.getvalue()}
        finally:
            _quiet_rm(out)

    def fit_native(self, param_spec: str, data_csv: str, response_map: str,
                   method: str = "lsqnonlin", doses: str = "") -> dict[str, Any]:
        """Calibrate parameters using SimBiology's NATIVE estimator (sbiofit) - the
        numerical optimization runs inside MATLAB, not over the bridge. ``param_spec``
        is ';'-joined 'name,lo,hi,scale' (scale 'lin'/'log'); ``data_csv`` is the
        observed data table (a Time column + one column per observed component);
        ``response_map`` is ';'-joined 'ModelComponent = DataColumn'; ``method`` is the
        sbiofit optimizer ('lsqnonlin'/'fmincon'/'particleswarm'/'ga'/...). Returns
        {columns: {parameter: [...], estimate: [...]}, matlab_log}."""
        import io
        out = _tmp(".csv")
        so = io.StringIO()
        try:
            self.eng.sb_fit(param_spec, os.path.abspath(data_csv), response_map,
                            method, doses or "", out, nargout=1, stdout=so, stderr=so)
            # The result table is "parameter,estimate": the first column is a NAME
            # (a string), so it must not be float-coerced the way _read_csv does.
            names, ests = [], []
            try:
                with open(out, newline="", encoding="utf-8") as fh:
                    rdr = csv.reader(fh)
                    next(rdr, None)  # header
                    for row in rdr:
                        if len(row) >= 2:
                            names.append(row[0].strip())
                            try:
                                ests.append(float(row[1]))
                            except (ValueError, TypeError):
                                ests.append(None)
            except OSError:
                pass
            return {"columns": {"parameter": names, "estimate": ests},
                    "matlab_log": so.getvalue()}
        finally:
            _quiet_rm(out)


def _quiet_rm(path: str) -> None:
    try:
        os.remove(path)
    except OSError:
        pass
