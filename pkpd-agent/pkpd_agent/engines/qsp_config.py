"""The per-model adapter for the downstream (2-7) tasks, loaded from DATA.

The engine (SimBiologyEngine), the loop, and the numerical routines (qsp_tasks:
numeric_fit_1d, select_to_moments, ir_mask, ...) are model-agnostic. What is
specific to a given QSP model is: which internal states are the clinical readouts,
which CSV column plays which role, which doses and parameters exist, the trial
timeline, and the real-world reference data. All of that lives in
projects/<name>/tasks.json; this module loads it into a ``QSPTaskConfig`` -- no
model's specifics are hardcoded here.

To port the agent to a NEW QSP model: add a projects/<name>/ folder with a
tasks.json (and a matching MATLAB readout script that emits the run_columns). No
engine, loop, or numerical-routine code changes.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from typing import Any, Optional

from . import qsp_tasks

# the fixed CSV column roles the readout script writes; a config maps each to a state
READOUT_ROLES = (
    "first_line_1", "first_line_2", "first_line_3", "first_line_4",
    "latched_1", "latched_2", "latched_3", "latched_4", "latched_5",
    "trajectory", "baseline",
)

_PROJECTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                             "..", "..", "projects"))
_ALIASES = {"ra": "vantage_ra", "vantage_ra": "vantage_ra"}


def _int_keys(d: dict) -> dict:
    """JSON object keys are strings; the clinical-trial 'weeks' maps use int weeks."""
    out = {}
    for k, v in (d or {}).items():
        if isinstance(v, dict) and "weeks" in v:
            v = dict(v)
            v["weeks"] = {int(w): arm for w, arm in v["weeks"].items()}
        out[k] = v
    return out


@dataclass
class QSPTaskConfig:
    """The model-specific adapter for the downstream tasks. Everything here is data
    from tasks.json; nothing in the engine, loop, or numerical routines is."""
    name: str
    disease: str
    readout_desc: str
    severity_readout: str
    readout_states: list[str]
    timeline: dict[str, float]
    run_columns: dict[str, Any]
    drugs: dict[str, Any]
    vpop_drivers: dict[str, Any]
    vpop_target: dict[str, Any]
    fit_params: dict[str, Any]
    design_targets: dict[str, Any]
    clinical_trials: dict[str, Any]
    refractory_target: dict[str, Any]
    validate_arms: dict[str, Any]
    flagship_protocol: dict[str, Any]
    calibrated_arms: list = field(default_factory=list)
    trial_objective: str = ""
    design_background: str = ""
    fit_default_arm: str = ""
    fit_target: dict[str, Any] = field(default_factory=dict)

    # -- convenience: bind the general engine functions to this config's columns -- #
    def summarize_run(self, res: dict) -> dict:
        return qsp_tasks.summarize_run(res, self.run_columns)

    def ir_mask(self, run: dict, acr_key: str = None, threshold: float = 3.2) -> dict:
        return qsp_tasks.ir_mask(run, self.run_columns, acr_key, threshold)

    def response_in_subgroup(self, run: dict, ids: set, roles=None) -> dict:
        return qsp_tasks.response_in_subgroup(run, ids, self.run_columns, roles)

    def trial_target(self, drug: str, week: int, correction: str = "raw"):
        return qsp_tasks.trial_target(self.clinical_trials, drug, week, correction)

    def validate(self) -> list[str]:
        """Cheap structural checks; returns a list of problems (empty = ok)."""
        problems = []
        if len(self.readout_states) != len(READOUT_ROLES):
            problems.append(f"readout_states needs {len(READOUT_ROLES)} names, "
                            f"got {len(self.readout_states)}")
        for key in ("baseline_day", "first_line_readout_day",
                    "second_line_readout_day"):
            if key not in self.timeline:
                problems.append(f"timeline missing '{key}'")
        return problems


def config_from_dict(d: dict) -> QSPTaskConfig:
    return QSPTaskConfig(
        name=d.get("name", "QSP model"),
        disease=d.get("disease", ""),
        readout_desc=d.get("readout_desc", "the clinical response"),
        severity_readout=d.get("severity_readout", "the severity score"),
        readout_states=list(d.get("readout_states", [])),
        timeline={k: float(v) for k, v in (d.get("timeline") or {}).items()},
        run_columns=dict(d.get("run_columns", {})),
        drugs=dict(d.get("drugs", {})),
        vpop_drivers=dict(d.get("vpop_drivers", {})),
        vpop_target=dict(d.get("vpop_target", {})),
        fit_params=dict(d.get("fit_params", {})),
        design_targets=dict(d.get("design_targets", {})),
        clinical_trials=_int_keys(d.get("clinical_trials", {})),
        refractory_target=dict(d.get("refractory_target", {})),
        validate_arms=dict(d.get("validate_arms", {})),
        flagship_protocol=dict(d.get("flagship_protocol", {})),
        calibrated_arms=list(d.get("calibrated_arms", [])),
        trial_objective=d.get("trial_objective", ""),
        design_background=d.get("design_background", ""),
        fit_default_arm=d.get("fit_default_arm", ""),
        fit_target=dict(d.get("fit_target", {})))


def load_tasks(name: str, projects_dir: str = None) -> QSPTaskConfig:
    """Load a project's tasks.json into a QSPTaskConfig. `name` is the folder."""
    base = projects_dir or _PROJECTS_DIR
    with open(os.path.join(base, name, "tasks.json"), encoding="utf-8") as fh:
        return config_from_dict(json.load(fh))


def get(name: str = "vantage_ra", projects_dir: str = None) -> QSPTaskConfig:
    key = (name or "vantage_ra").lower()
    folder = _ALIASES.get(key, key)
    base = projects_dir or _PROJECTS_DIR
    if not os.path.isfile(os.path.join(base, folder, "tasks.json")):
        known = sorted(d for d in os.listdir(base)
                       if os.path.isfile(os.path.join(base, d, "tasks.json"))) \
            if os.path.isdir(base) else []
        raise KeyError(f"unknown project '{name}'. Known: {known}. "
                       "Add a projects/<name>/tasks.json.")
    return load_tasks(folder, base)
