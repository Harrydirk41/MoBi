"""Human-readable validation of a project's tasks.json / spec.json against the model.

The config loader silently defaults a missing/mistyped field to empty, so a typo in a
role name (a `vpop_driver` that is not actually a model parameter, a `readout_state`
that is not a species) makes a whole task quietly do nothing. For a traditional
modeler editing JSON by hand that is a trap. This module checks the config's names
against the real model (network.json) and returns plain-English problems - errors
(will not work) and warnings (probably not what you meant / still a stub).

Pure Python, no LLM, no MATLAB: it compares name sets. Run it after building a config
and before running the agent.
"""

from __future__ import annotations

from typing import Any

_STUB = "TODO"


def _names(network: dict, key: str) -> set:
    return {x.get("name") for x in network.get(key, []) if x.get("name")}


def validate_project(tasks: dict, spec: dict, network: dict) -> dict[str, list[str]]:
    """Check a tasks.json + spec.json against a network.json. Returns
    {'errors': [...], 'warnings': [...]} of plain-English messages (empty = clean)."""
    species = _names(network, "species")
    params = _names(network, "parameters")
    errors: list[str] = []
    warnings: list[str] = []

    def check_params(field: str, kind: str, warn_only: bool = False):
        for name in (tasks.get(field) or {}):
            if name not in params:
                msg = (f"{kind} '{name}' is not a parameter in the model - it will be "
                       f"silently ignored. Check the spelling against the model's "
                       f"parameter names, or remove it from '{field}'.")
                (warnings if warn_only else errors).append(msg)

    # 1. task parameter roles must be real model parameters (the main foot-gun)
    check_params("vpop_drivers", "vpop_driver")
    check_params("fit_params", "fit_param")
    check_params("design_targets", "design_target")

    # 2. readout states should be model species (flags/scores are species/derived states)
    for st in (tasks.get("readout_states") or []):
        if st not in species:
            warnings.append(f"readout_state '{st}' is not a model species - make sure "
                            "the readout script emits it (it may be a derived output).")

    # 3. run_columns structure
    rc = tasks.get("run_columns") or {}
    if not rc:
        errors.append("run_columns is missing - the agent cannot read any response "
                      "from a Vpop run without it.")
    else:
        if not (rc.get("severity") or {}).get("baseline") or \
                not (rc.get("severity") or {}).get("readout"):
            warnings.append("run_columns.severity needs both 'baseline' and 'readout' "
                            "column names for the Vpop / calibration tasks.")
        if not rc.get("first_line"):
            warnings.append("run_columns.first_line is empty - the trial / design tasks "
                            "have no response endpoints to read.")

    # 4. timeline keys the tasks rely on
    for key in ("baseline_day", "first_line_readout_day", "second_line_readout_day"):
        if key not in (tasks.get("timeline") or {}):
            warnings.append(f"timeline is missing '{key}' - a sensible default will be "
                            "used, but set it to match your model's events.")

    # 5. external fields still left as drafter stubs
    for field in ("drugs", "vpop_target", "clinical_trials", "refractory_target"):
        v = tasks.get(field)
        if v is None or (isinstance(v, str) and v.startswith(_STUB)) or v == {}:
            warnings.append(f"'{field}' is not filled in yet - the drafter leaves this "
                            "for you (external clinical data or your model's dose names).")

    # 6. spec: gsa_top and readout_targets sanity (structure benchmarks)
    for name in (spec.get("gsa_top") or []):
        if name not in params:
            warnings.append(f"spec gsa_top '{name}' is not a model parameter - the "
                            "sensitivity benchmark answer key expects real parameters.")
    if not spec.get("readout_targets"):
        warnings.append("spec.readout_targets is empty - the readout-mapping benchmark "
                        "needs the readout rule name(s).")

    return {"errors": errors, "warnings": warnings}


def format_report(result: dict) -> str:
    """Render the validate_project result as a readable block."""
    lines = []
    if result["errors"]:
        lines.append("ERRORS (must fix - these will not work):")
        lines += [f"  x {m}" for m in result["errors"]]
    if result["warnings"]:
        lines.append("WARNINGS (check these are intended):")
        lines += [f"  ! {m}" for m in result["warnings"]]
    if not lines:
        return "OK - config is consistent with the model."
    return "\n".join(lines)
