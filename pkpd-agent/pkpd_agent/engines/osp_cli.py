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
                      workdir: str | None = None,
                      prune_simulations: bool = False,
                      target_molecule: str | None = None,
                      target_molecules: list[str] | None = None) -> dict[str, Any]:
        """snapshot (optionally edited) -> .pksim5 -> run -> predicted profiles.

        ``edits`` is the structured spec from ``snapshot_edit`` (parameters,
        calculation_methods, processes). ``param_overrides`` is a shorthand for
        ``edits={"parameters": ...}`` kept for backwards compatibility.

        ``simulations`` restricts which simulations are *run* (export -s).
        ``prune_simulations`` additionally drops the other simulations from the
        snapshot BEFORE snap, so PK-Sim only *builds* the ones requested - a big
        speedup during optimization, where snap otherwise rebuilds every
        simulation on every objective evaluation.

        ``target_molecules`` (metabolite cascades): pull SEVERAL molecules'
        plasma from the SAME run - parent + each daughter metabolite - each with
        its own molecular weight. When given, the result carries an extra
        ``profiles_by_molecule`` = {molecule: [PredictedProfile...]} so a single
        snap+export scores the whole cascade, not one molecule at a time.

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
        if prune_simulations and simulations:
            snap = self._prune_simulations(snap, simulations)

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

        profiles = self._parse_results(out_dir, mol_weight, target_molecule)
        by_molecule = None
        if target_molecules:
            mw_by_name = self._mol_weight_by_name(snap)
            by_molecule = self._parse_results_multi(out_dir, mw_by_name,
                                                    target_molecules, mol_weight)
        if not self.keep_workdir and not workdir:
            shutil.rmtree(wd, ignore_errors=True)

        n_edits = (len(applied.get("parameters", {}))
                   + len(applied.get("calculation_methods", {}))
                   + len([v for v in applied.get("processes", {}).values()
                          if v == "disabled"]))
        result = {"ok": bool(profiles),
                  "message": (f"ran {len(profiles)} simulation(s); "
                              f"applied {n_edits} edit(s)")
                             if profiles else "no result CSVs parsed",
                  "profiles": profiles, "project": project,
                  "edits_applied": applied, "logs": logs, "mol_weight": mol_weight}
        if by_molecule is not None:
            result["profiles_by_molecule"] = by_molecule
            # a cascade run "ok" if at least one requested molecule was parsed
            result["ok"] = result["ok"] or any(by_molecule.values())
        return result

    # -- snapshot helpers ------------------------------------------------ #
    @staticmethod
    def _prune_simulations(snap: dict, keep: list[str]) -> dict:
        """Return a shallow copy of the snapshot whose Simulations are limited to
        ``keep`` (so snap builds only those). Other top-level blocks (Individuals,
        Protocols, ObservedData, ...) are referenced by name and harmless to keep.
        Falls back to the full set if the filter would empty it."""
        keepset = set(keep)
        sims = [s for s in (snap.get("Simulations") or [])
                if s.get("Name") in keepset]
        if not sims:
            return snap
        pruned = dict(snap)
        pruned["Simulations"] = sims
        # ParameterIdentifications reference simulations by name; drop them so a
        # dangling reference cannot trip the mapper (we identify parameters here).
        if pruned.get("ParameterIdentifications"):
            pruned["ParameterIdentifications"] = []
        return pruned

    @staticmethod
    def _mol_weight(snap: dict) -> float | None:
        comp = (snap.get("Compounds") or [{}])[0]
        for p in comp.get("Parameters") or []:
            if p.get("Name") == "Molecular weight":
                return p.get("Value")
        return None

    @staticmethod
    def _mol_weight_by_name(snap: dict) -> dict[str, float | None]:
        """Molecular weight per compound name - a metabolite cascade has several
        compounds, each with its own MW (needed to convert its mol/L result to
        mg/L). Falls back to None (no conversion) when absent."""
        out: dict[str, float | None] = {}
        for comp in snap.get("Compounds") or []:
            nm = comp.get("Name")
            if not nm:
                continue
            mw = None
            for p in comp.get("Parameters") or []:
                if p.get("Name") == "Molecular weight":
                    mw = p.get("Value")
                    break
            out[nm] = mw
        return out

    # -- CSV parsing ----------------------------------------------------- #
    def _parse_results(self, out_dir: str, mol_weight: float | None,
                       target_molecule: str | None = None) \
            -> list[PredictedProfile]:
        profiles = []
        for csv_path in sorted(glob.glob(os.path.join(out_dir, "*", "*-Results.csv"))):
            sim_name = os.path.basename(os.path.dirname(csv_path))
            prof = self._parse_one(csv_path, sim_name, mol_weight, target_molecule)
            if prof:
                profiles.append(prof)
        return profiles

    def _parse_results_multi(self, out_dir: str,
                             mw_by_name: dict[str, float | None],
                             molecules: list[str],
                             default_mw: float | None = None) \
            -> dict[str, list[PredictedProfile]]:
        """{molecule: [PredictedProfile...]} - one concentration column per
        requested molecule per simulation CSV, each converted with that
        molecule's own molecular weight. This is how a metabolite cascade is
        scored from ONE run: parent + every daughter, side by side."""
        by_mol: dict[str, list[PredictedProfile]] = {m: [] for m in molecules}
        for csv_path in sorted(glob.glob(os.path.join(out_dir, "*", "*-Results.csv"))):
            sim_name = os.path.basename(os.path.dirname(csv_path))
            for m in molecules:
                mw = mw_by_name.get(m, default_mw)
                if mw is None:
                    mw = default_mw
                prof = self._parse_one(csv_path, sim_name, mw, target_molecule=m,
                                       require_molecule_match=True)
                if prof:
                    by_mol[m].append(prof)
        return by_mol

    def _parse_one(self, csv_path: str, sim_name: str,
                   mol_weight: float | None,
                   target_molecule: str | None = None,
                   require_molecule_match: bool = False) -> PredictedProfile | None:
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
        # (the target/victim molecule's, when a DDI run has several compounds)
        c_idx = self._pick_conc_column(header, target_molecule,
                                       require_molecule_match)
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
    def _pick_conc_column(header: list[str], target_molecule: str | None = None,
                          require_molecule_match: bool = False) -> int | None:
        conc_cols = []
        for i, h in enumerate(header):
            u = _unit_in_brackets(h)
            if u in _MASS_PER_L or "mol/l" in u:
                conc_cols.append((i, h))
        if not conc_cols:
            return None
        # in a DDI run the CSV has a column per compound - pick the target
        # (victim) molecule's plasma, else the perpetrator's would be scored.
        if target_molecule:
            tl = target_molecule.lower()
            for i, h in conc_cols:
                hl = h.lower()
                if tl in hl and "peripheralvenousblood" in hl and "plasma" in hl:
                    return i
            for i, h in conc_cols:
                if tl in h.lower():
                    return i
            # metabolite cascade: each molecule must match its OWN column; do NOT
            # fall back to the parent's plasma (that would score every daughter
            # against the parent). No name match -> this molecule isn't here.
            if require_molecule_match:
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
