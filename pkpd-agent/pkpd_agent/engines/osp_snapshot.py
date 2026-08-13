"""Parse an OSP (PK-Sim/MoBi) snapshot JSON into agent-ready data.

A snapshot is a *recipe* (see the discussion in the README): it does not carry
the runnable ODE model, but it DOES carry, in plain readable form:

  * observed clinical data      -> tidy time/concentration profiles
  * the compound parameter file -> physchem / ADME / metabolism parameters
  * the modeling choices        -> distribution method, processes (enzymes/
                                   transporters), and which inputs were
                                   measured vs fitted (the "answer key")

This module extracts all three with the standard library only, so it runs
anywhere and needs no OSP install. Simulation still needs a .pkml + ospsuite;
this is the data/analysis half of the pipeline.
"""

from __future__ import annotations

import csv
import json
import math
import os
from dataclasses import dataclass, field
from typing import Any


# --------------------------------------------------------------------------- #
# Observed profiles
# --------------------------------------------------------------------------- #

@dataclass
class ObservedProfile:
    name: str
    time_h: list[float]
    conc: list[float]           # in the snapshot's native concentration unit
    conc_unit: str
    mol_weight: float | None

    def conc_mg_L(self) -> list[float]:
        """Convert to mg/L. Molar units (µmol/l) need the molecular weight;
        mass units are returned as-is."""
        u = (self.conc_unit or "").lower()
        if self.mol_weight and "mol" in u:
            factor = self.mol_weight / 1000.0        # µmol/L * g/mol / 1000 = mg/L
            if u.startswith("nmol"):
                factor /= 1000.0
            return [c * factor for c in self.conc]
        return list(self.conc)

    def nca(self) -> dict[str, Any]:
        pts = sorted(zip(self.time_h, self.conc_mg_L()))
        t = [p[0] for p in pts]
        c = [p[1] for p in pts]
        if not t:
            return {"study": self.name, "n": 0}
        cmax = max(c)
        tmax = t[c.index(cmax)]
        auc = sum((t[i] - t[i - 1]) * (c[i] + c[i - 1]) / 2 for i in range(1, len(t)))
        thalf = None
        tail = [(t[i], c[i]) for i in range(len(t)) if c[i] > 0][-3:]
        if len(tail) >= 2 and tail[0][1] > tail[-1][1] and tail[-1][0] > tail[0][0]:
            k = (math.log(tail[0][1]) - math.log(tail[-1][1])) / (tail[-1][0] - tail[0][0])
            if k > 0:
                thalf = round(math.log(2) / k, 3)
        return {
            "study": self.name, "n": len(t),
            "c_max_mg_L": round(cmax, 4), "t_max_h": round(tmax, 3),
            "auc_mg_h_L": round(auc, 4), "t_half_h": thalf,
        }


# --------------------------------------------------------------------------- #
# Snapshot
# --------------------------------------------------------------------------- #

@dataclass
class OSPSnapshot:
    data: dict[str, Any]

    @classmethod
    def from_file(cls, path: str) -> "OSPSnapshot":
        with open(path, encoding="utf-8") as fh:
            return cls(json.load(fh))

    # -- inventory ------------------------------------------------------- #
    def summary(self) -> dict[str, Any]:
        comp = self.data.get("Compounds") or []
        return {
            "compound": comp[0].get("Name") if comp else None,
            "n_observed_datasets": len(self.data.get("ObservedData") or []),
            "n_compounds": len(comp),
            "n_simulations": len(self.data.get("Simulations") or []),
            "n_protocols": len(self.data.get("Protocols") or []),
            "n_expression_profiles": len(self.data.get("ExpressionProfiles") or []),
        }

    # -- observed data --------------------------------------------------- #
    def observed_profiles(self) -> list[ObservedProfile]:
        out: list[ObservedProfile] = []
        for od in self.data.get("ObservedData") or []:
            base = od.get("BaseGrid") or {}
            times = list(base.get("Values") or [])
            tunit = base.get("Unit", "")
            if tunit and tunit.lower().startswith("min"):
                times = [t / 60.0 for t in times]        # min -> h
            cols = od.get("Columns") or []
            if not cols:
                continue
            col = cols[0]
            vals = list(col.get("Values") or [])
            n = min(len(times), len(vals))
            if n == 0:
                continue
            mw = (col.get("DataInfo") or {}).get("MolWeight")
            out.append(ObservedProfile(
                name=od.get("Name", "unnamed"),
                time_h=times[:n], conc=vals[:n],
                conc_unit=col.get("Unit", ""), mol_weight=mw,
            ))
        return out

    # -- compound parameters -------------------------------------------- #
    def compound_parameters(self) -> list[dict[str, Any]]:
        comp = self.data.get("Compounds") or []
        if not comp:
            return []
        found: list[dict[str, Any]] = []

        def walk(obj: Any) -> None:
            if isinstance(obj, dict):
                if "Name" in obj and isinstance(obj.get("Value"), (int, float)):
                    found.append({"parameter": obj["Name"], "value": obj["Value"],
                                  "unit": obj.get("Unit", "")})
                for v in obj.values():
                    walk(v)
            elif isinstance(obj, list):
                for v in obj:
                    walk(v)

        walk(comp[0])
        return found

    # -- modeling choices (the "answer key") ---------------------------- #
    def modeling_choices(self) -> dict[str, Any]:
        comp = (self.data.get("Compounds") or [{}])[0]
        # measured-vs-fitted tags on the key input properties
        fit_vs_measured: dict[str, str] = {}
        for prop in ("Lipophilicity", "FractionUnbound", "Solubility",
                     "IntestinalPermeability"):
            entries = comp.get(prop) or []
            if entries and isinstance(entries[0], dict):
                fit_vs_measured[prop] = entries[0].get("Name", "")
        # processes (enzymes / transporters) and their data source
        processes = []
        for p in comp.get("Processes") or []:
            processes.append({
                "molecule": p.get("Molecule") or p.get("MoleculeName"),
                "data_source": p.get("DataSource") or p.get("InternalName"),
            })
        # distribution / permeability calculation methods
        methods = []
        for m in comp.get("CalculationMethods") or []:
            methods.append(m if isinstance(m, str) else m.get("Name"))
        sims = self.data.get("Simulations") or []
        return {
            "calculation_methods": methods,
            "processes": processes,
            "fit_vs_measured": fit_vs_measured,
            "model_type": sims[0].get("Model") if sims else None,
        }

    # -- NCA ------------------------------------------------------------- #
    def nca_table(self) -> list[dict[str, Any]]:
        return [p.nca() for p in self.observed_profiles()]

    # -- write CSVs ------------------------------------------------------ #
    def write_csvs(self, outdir: str) -> dict[str, str]:
        os.makedirs(outdir, exist_ok=True)
        comp_name = (self.summary().get("compound") or "compound").replace(" ", "_")

        obs_path = os.path.join(outdir, f"{comp_name}_observed.csv")
        with open(obs_path, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["study", "time_h", "conc_native", "conc_unit", "conc_mg_L"])
            for p in self.observed_profiles():
                mg = p.conc_mg_L()
                for t, c, m in zip(p.time_h, p.conc, mg):
                    w.writerow([p.name, round(t, 4), c, p.conc_unit, round(m, 6)])

        par_path = os.path.join(outdir, f"{comp_name}_compound_params.csv")
        with open(par_path, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["parameter", "value", "unit"])
            for p in self.compound_parameters():
                w.writerow([p["parameter"], p["value"], p["unit"]])

        return {"observed_csv": obs_path, "params_csv": par_path}
