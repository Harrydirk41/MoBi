"""Headless PBPK engine via PK-Sim's command-line interface (PKSim.CLI).

This is the real OSP model-building engine, no GUI and no ospsuite required. It
reproduces exactly what the PK-Sim "Load project from snapshot -> Run" dialog
does, but scriptably, so the agent can EDIT the snapshot and get simulated
concentration-time curves back:

    snapshot JSON  --(PKSim.CLI snap -p)-->  .pksim5 project
    .pksim5        --(PKSim.CLI export -r -c)-->  per-simulation Results.csv
    Results.csv    --(this module)-->  predicted profiles (time h, conc mg/L)

Because the input is the snapshot (which we fully control as JSON), the agent's
action space is the whole PBPK model-choice space: partition method, processes
(enzymes/transporters), and every compound parameter -- not just a fixed .pkml.

Verified against PK-Sim 12.3:
  * snap   -i <in-folder> -o <out-folder> -p       (snapshot JSON -> .pksim5)
  * export -p <project.pksim5> -o <out> -r -c       (run + results CSV)
  * Results.csv columns: "IndividualId","Time [min]",
        "Organism|PeripheralVenousBlood|<drug>|Plasma (Peripheral Venous Blood) [µmol/l]", ...
"""

from __future__ import annotations

import csv
import glob
import json
import os
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from typing import Any


# --------------------------------------------------------------------------- #
# unit parsing
# --------------------------------------------------------------------------- #

_MASS_PER_L = {"mg/l": 1.0, "µg/l": 1e-3, "ug/l": 1e-3, "ng/ml": 1e-3,
               "ng/l": 1e-6, "g/l": 1e3}


def _unit_in_brackets(header: str) -> str:
    m = re.findall(r"\[([^\]]*)\]", header)
    return (m[-1] if m else "").strip().lower()


def _conc_to_mg_L(value: float, unit: str, mol_weight: float | None) -> float | None:
    u = unit.strip().lower()
    if u in _MASS_PER_L:
        return value * _MASS_PER_L[u]
    if "mol" in u and mol_weight:                        # µmol/l, nmol/l ...
        f = mol_weight / 1000.0                           # µmol/L * g/mol / 1000 = mg/L
        if u.startswith("nmol"):
            f /= 1000.0
        return value * f
    return None                                          # unknown -> caller skips


def _time_to_h(value: float, unit: str) -> float:
    u = unit.strip().lower()
    if u.startswith("min"):
        return value / 60.0
    if u.startswith("day"):
        return value * 24.0
    if u.startswith("s"):
        return value / 3600.0
    return value                                          # assume hours


# --------------------------------------------------------------------------- #
# result
# --------------------------------------------------------------------------- #

@dataclass
class PredictedProfile:
    simulation: str
    study: str | None
    route: str | None
    dose: str | None
    time_h: list[float]
    conc_mg_L: list[float]


@dataclass
class OSPCli:
    pksim_cli_path: str
    timeout_s: int = 900
    keep_workdir: bool = False

    # -- public API ------------------------------------------------------ #
    def simulation_names(self, snapshot_path: str) -> list[str]:
        with open(snapshot_path, encoding="utf-8") as fh:
            data = json.load(fh)
        return [s.get("Name") for s in (data.get("Simulations") or [])
                if s.get("Name")]

    def build_and_run(self, snapshot_path: str,
                      edits: dict | None = None,
                      param_overrides: dict[str, float] | None = None,
                      simulations: list[str] | None = None,
                      workdir: str | None = None) -> dict[str, Any]:
        """snapshot (optionally edited) -> .pksim5 -> run -> predicted profiles.

        ``edits`` is the structured spec from ``snapshot_edit`` (parameters,
        calculation_methods, processes). ``param_overrides`` is a shorthand for
        ``edits={"parameters": ...}`` kept for backwards compatibility.

        Returns {ok, message, profiles:[PredictedProfile...], project, logs}.
        """
        if not self.pksim_cli_path or not os.path.exists(self.pksim_cli_path):
            return {"ok": False, "message":
                    f"PKSim.CLI not found at {self.pksim_cli_path!r}; "
                    "set config.pksim_cli_path or PKPD_PKSIM_CLI.", "profiles": []}

        from .snapshot_edit import apply_edits
        with open(snapshot_path, encoding="utf-8") as fh:
            snap = json.load(fh)
        mol_weight = self._mol_weight(snap)

        edit_spec = dict(edits or {})
        if param_overrides:
            edit_spec.setdefault("parameters", {})
            edit_spec["parameters"] = {**param_overrides, **edit_spec["parameters"]}
        snap, applied = apply_edits(snap, edit_spec)

        wd = workdir or tempfile.mkdtemp(prefix="ospcli_")
        in_dir = os.path.join(wd, "in")
        proj_dir = os.path.join(wd, "proj")
        out_dir = os.path.join(wd, "out")
        for d in (in_dir, proj_dir, out_dir):
            os.makedirs(d, exist_ok=True)

        stem = os.path.splitext(os.path.basename(snapshot_path))[0]
        snap_in = os.path.join(in_dir, f"{stem}.json")
        with open(snap_in, "w", encoding="utf-8") as fh:
            json.dump(snap, fh, ensure_ascii=False)

        logs = []
        # 1. snapshot -> project
        r1 = self._run(["snap", "-i", in_dir, "-o", proj_dir, "-p"])
        logs.append(r1)
        if r1["returncode"] != 0:
            return self._fail("snap (snapshot->project) failed", logs, wd)
        projects = glob.glob(os.path.join(proj_dir, "*.pksim5"))
        project = next((p for p in projects
                        if os.path.splitext(os.path.basename(p))[0] == stem),
                       projects[0] if projects else None)
        if not project:
            return self._fail("no .pksim5 produced by snap", logs, wd)

        # 2. project -> run + results CSV (optionally a subset of simulations)
        export_args = ["export", "-p", project, "-o", out_dir, "-r", "-c"]
        if simulations:
            export_args += ["-s"] + list(simulations)
        r2 = self._run(export_args)
        logs.append(r2)
        if r2["returncode"] != 0:
            return self._fail("export (run) failed", logs, wd)

        profiles = self._parse_results(out_dir, mol_weight)
        if not self.keep_workdir and not workdir:
            shutil.rmtree(wd, ignore_errors=True)

        n_edits = (len(applied.get("parameters", {}))
                   + len(applied.get("calculation_methods", {}))
                   + len([v for v in applied.get("processes", {}).values()
                          if v == "disabled"]))
        return {"ok": bool(profiles),
                "message": (f"ran {len(profiles)} simulation(s); "
                            f"applied {n_edits} edit(s)")
                           if profiles else "no result CSVs parsed",
                "profiles": profiles, "project": project,
                "edits_applied": applied, "logs": logs, "mol_weight": mol_weight}

    # -- snapshot helpers ------------------------------------------------ #
    @staticmethod
    def _mol_weight(snap: dict) -> float | None:
        comp = (snap.get("Compounds") or [{}])[0]
        for p in comp.get("Parameters") or []:
            if p.get("Name") == "Molecular weight":
                return p.get("Value")
        return None

    # -- CSV parsing ----------------------------------------------------- #
    def _parse_results(self, out_dir: str, mol_weight: float | None) \
            -> list[PredictedProfile]:
        profiles = []
        for csv_path in sorted(glob.glob(os.path.join(out_dir, "*", "*-Results.csv"))):
            sim_name = os.path.basename(os.path.dirname(csv_path))
            prof = self._parse_one(csv_path, sim_name, mol_weight)
            if prof:
                profiles.append(prof)
        return profiles

    def _parse_one(self, csv_path: str, sim_name: str,
                   mol_weight: float | None) -> PredictedProfile | None:
        with open(csv_path, encoding="utf-8", errors="replace", newline="") as fh:
            rows = list(csv.reader(fh))
        if len(rows) < 2:
            return None
        header = rows[0]
        # time column
        t_idx = next((i for i, h in enumerate(header)
                      if h.lower().startswith("time")), None)
        if t_idx is None:
            return None
        t_unit = _unit_in_brackets(header[t_idx])
        # concentration column: peripheral venous plasma, a concentration unit
        c_idx = self._pick_conc_column(header)
        if c_idx is None:
            return None
        c_unit = _unit_in_brackets(header[c_idx])

        # aggregate across IndividualId by mean-per-time (single-individual sims
        # are a no-op; population sims collapse to the mean profile)
        acc: dict[float, list[float]] = {}
        for row in rows[1:]:
            if len(row) <= max(t_idx, c_idx):
                continue
            try:
                t = float(row[t_idx]); c = float(row[c_idx])
            except ValueError:
                continue
            acc.setdefault(t, []).append(c)

        times = sorted(acc)
        time_h, conc_mg_L = [], []
        for t in times:
            mean_c = sum(acc[t]) / len(acc[t])
            mg = _conc_to_mg_L(mean_c, c_unit, mol_weight)
            if mg is None:
                continue
            time_h.append(round(_time_to_h(t, t_unit), 6))
            conc_mg_L.append(round(mg, 10))
        if not time_h:
            return None

        study, route, dose = self._parse_sim_name(sim_name)
        return PredictedProfile(sim_name, study, route, dose, time_h, conc_mg_L)

    @staticmethod
    def _pick_conc_column(header: list[str]) -> int | None:
        conc_cols = []
        for i, h in enumerate(header):
            u = _unit_in_brackets(h)
            if u in _MASS_PER_L or "mol/l" in u:
                conc_cols.append((i, h))
        if not conc_cols:
            return None
        # prefer peripheral venous plasma (matches observed clinical matrix)
        for i, h in conc_cols:
            hl = h.lower()
            if "peripheralvenousblood" in hl and "plasma" in hl:
                return i
        for i, h in conc_cols:
            if "plasma" in h.lower():
                return i
        return conc_cols[0][0]

    @staticmethod
    def _parse_sim_name(name: str) -> tuple[str | None, str | None, str | None]:
        study = name.split(",")[0].strip() or None
        route = None
        m = re.search(r"(?:^|[\s_,])(iv|po|oral)(?:$|[\s_,])", name, re.IGNORECASE)
        if m:
            route = "PO" if m.group(1).lower() in ("po", "oral") else "IV"
        dose = None
        # weight-based e.g. "0.05 mg_kg" / "20µg_kg"  (underscore = per)
        m = re.search(r"([\d.]+)\s*(µg|ug|mg|g)\s*[_/]\s*kg", name, re.IGNORECASE)
        if m:
            dose = f"{m.group(1)} {m.group(2).replace('ug','µg')}/kg"
        else:                                             # absolute e.g. "1 mg"
            m = re.search(r"([\d.]+)\s*(µg|ug|mg|g)\b", name, re.IGNORECASE)
            if m:
                dose = f"{m.group(1)} {m.group(2).replace('ug','µg')}"
        return study, route, dose

    # -- subprocess ------------------------------------------------------ #
    def _run(self, args: list[str]) -> dict[str, Any]:
        cmd = [self.pksim_cli_path] + args
        try:
            p = subprocess.run(cmd, capture_output=True, text=True,
                               encoding="utf-8", errors="replace",
                               timeout=self.timeout_s)
            return {"cmd": " ".join(args), "returncode": p.returncode,
                    "stdout": (p.stdout or "")[-2000:], "stderr": (p.stderr or "")[-2000:]}
        except subprocess.TimeoutExpired:
            return {"cmd": " ".join(args), "returncode": -1,
                    "stdout": "", "stderr": f"timeout after {self.timeout_s}s"}
        except Exception as exc:                          # noqa: BLE001
            return {"cmd": " ".join(args), "returncode": -1,
                    "stdout": "", "stderr": f"{type(exc).__name__}: {exc}"}

    def _fail(self, msg: str, logs: list, wd: str) -> dict[str, Any]:
        last = logs[-1] if logs else {}
        detail = (last.get("stderr") or last.get("stdout") or "").strip()
        if not self.keep_workdir:
            shutil.rmtree(wd, ignore_errors=True)
        return {"ok": False, "message": f"{msg}: {detail[:400]}",
                "profiles": [], "logs": logs}
