"""One place for everything model-specific — the adapter that points the agent at a
particular QSP model.

The engine (SimBiologyEngine), the loop, and the numerical routines (numeric_fit_1d,
select_to_moments, ir_mask, ...) are model-agnostic. What is specific to a given
QSP model is: which internal states are the clinical readouts, which doses and
parameters exist, the trial timeline, and the real-world reference data. This module
captures that surface as a ``QSPModelConfig`` and provides the Vantage RA instance.

To port the agent to a NEW QSP model, build a new ``QSPModelConfig`` (no engine,
loop, or numerical-routine code changes):

  * ``readout_states`` - the 11 model state names, in role order, that
    ``sb_run_vpop`` reads (see ``READOUT_ROLES``); a new model has different events.
  * ``timeline``       - the readout days the model's events fire at.
  * ``drugs`` / ``vpop_drivers`` / ``fit_params`` / ``design_targets`` - the
    catalogs the task tools present to the agent (dose names, disease-driver
    parameters, calibratable parameters, druggable pathways).
  * ``clinical_reference`` - the module of real trial data to score against.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from . import osp_ra_trial
from . import ra_clinical_reference

# the fixed CSV column roles sb_run_vpop writes; a config maps each to a model state
READOUT_ROLES = (
    "first_line_1", "first_line_2", "first_line_3", "first_line_4",  # day-284 ACR flags
    "latched_1", "latched_2", "latched_3", "latched_4", "latched_5",  # day-600 flags
    "trajectory", "baseline",                                          # DAS28 states
)


@dataclass
class QSPModelConfig:
    """The model-specific adapter. Everything here changes per QSP model; nothing in
    the engine, loop, or numerical routines does."""
    name: str
    readout_states: list[str]          # 11 state names, role order (READOUT_ROLES)
    timeline: dict[str, float]         # baseline / first-line / second-line readout days
    drugs: dict[str, Any]              # drug formulary (code -> {drug, mechanism, doses})
    vpop_drivers: dict[str, Any]       # disease-driver params to sample for a Vpop
    vpop_target: dict[str, Any]        # target baseline-score distribution
    fit_params: dict[str, Any]         # calibratable PD parameters (name -> spec)
    design_targets: dict[str, Any]     # druggable pathways (param -> {pathway, ...})
    clinical_reference: Any = None     # module/object of real trial data to score on
    notes: str = ""

    def validate(self) -> list[str]:
        """Cheap structural checks; returns a list of problems (empty = ok)."""
        problems = []
        if len(self.readout_states) != len(READOUT_ROLES):
            problems.append(
                f"readout_states needs {len(READOUT_ROLES)} names, "
                f"got {len(self.readout_states)}")
        for key in ("baseline_day", "first_line_readout_day", "second_line_readout_day"):
            if key not in self.timeline:
                problems.append(f"timeline missing '{key}'")
        return problems


# --------------------------------------------------------------------------- #
# The Vantage RA QSP model — the one adapter that exists today. Its catalogs live
# in osp_ra_trial / ra_clinical_reference; this gathers them into the config.
# --------------------------------------------------------------------------- #

VANTAGE_RA = QSPModelConfig(
    name="Vantage RA QSP Model v1.0",
    readout_states=[
        "ACR20", "ACR50", "ACR70", "Remission",                       # first-line, day 284
        "MTX_NonResp", "MTX_NonResp_TCZ_ACR20", "MTX_NonResp_TCZ_ACR50",
        "MTX_NonResp_TCZ_ACR70", "MTX_NonResp_TCZ_Rem",               # latched, day 600
        "DAS28_CRP", "DAS28_BL",                                       # trajectory, baseline
    ],
    timeline={
        "baseline_captured_day": 199.0,
        "treatment_start_day": 200.0,
        "baseline_day": 200.0,
        "first_line_readout_day": 284.0,
        "second_line_readout_day": 600.0,
    },
    drugs=osp_ra_trial.DRUG_CATALOG,
    vpop_drivers=osp_ra_trial.VPOP_PARAMS,
    vpop_target=osp_ra_trial.VPOP_TARGET,
    fit_params=osp_ra_trial.FIT_PARAMS,
    design_targets=osp_ra_trial.DESIGN_TARGETS,
    clinical_reference=ra_clinical_reference,
    notes="Readout state names, drug catalog and references are all RA-specific; the "
          "engine/loop/numerical routines are not.",
)


REGISTRY: dict[str, QSPModelConfig] = {"vantage_ra": VANTAGE_RA}


def get(name: str = "vantage_ra") -> QSPModelConfig:
    return REGISTRY[name]
