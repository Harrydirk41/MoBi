"""LLM-loop tools for the RA QSP virtual-trial task: inspect the trial + run-a-protocol.

The SimBiology analogue of osp_ddi_loop_tools. Where the DDI loop tunes interaction
parameters, this loop designs an in-silico TRIAL: the agent chooses a treatment
protocol (which drugs, and - for a second-line agent - when to switch it in) and
runs it across the virtual population; the model's own events compute the clinical
response.

  * ``ra_inspect``    (observe) - the disease/readout, the trial timeline, the drug
    formulary with mechanisms and dose names, the CALIBRATED reference arms with
    their known response rates (to validate the harness), and the held-out
    objective (predict the second-line response in MTX-inadequate responders).
  * ``ra_run_trial``  (act) - apply a protocol, run the Vpop, and return the
    model-computed first-line and second-line response RATES.

The disease model and virtual population are GIVEN and fixed; the unknown is the
PROTOCOL (which therapy, given how). Registered per-task, so the handlers close
over the trial context. The held-out target is NOT exposed to the agent - it
reasons its way to the right therapy and the example script scores the final
prediction against the held-out truth.
"""

from __future__ import annotations

from typing import Any

from ..engines.simbiology import SimBiologyEngine
from ..engines import osp_ra_trial
from .registry import Tool, ToolRegistry, ToolResult


def register_ra_trial_loop_tools(registry: ToolRegistry, config, ctx: dict) -> None:
    """ctx: {sb, vpop, calibrated_arms, objective, limit, stop_time,
    baseline_day, readout_day}."""
    sb: SimBiologyEngine = ctx["sb"]
    vpop: str = ctx["vpop"]
    calibrated_arms: list = ctx.get("calibrated_arms") or []
    objective: str = ctx.get("objective") or "Predict the second-line response."
    limit: int = int(ctx.get("limit") or 50)
    stop_time: float = float(ctx.get("stop_time") or 700.0)
    baseline_day: float = float(ctx.get("baseline_day") or 200.0)
    readout_day: float = float(ctx.get("readout_day") or 284.0)

    def _formulary() -> list[dict]:
        out = []
        for code, d in osp_ra_trial.DRUG_CATALOG.items():
            out.append({
                "code": code, "drug": d["drug"], "modality": d["modality"],
                "mechanism": d["mechanism"], "dose_names": d["doses"],
            })
        return out

    # -- observe -------------------------------------------------------- #
    def inspect(args: dict, session) -> ToolResult:
        return ToolResult.success(
            "RA QSP virtual-trial task: choose a treatment protocol and run it "
            "across the virtual population; the model's events compute the clinical "
            "response (ACR20/50/70, DAS28 remission).",
            objective=objective,
            disease="rheumatoid arthritis; readout DAS28-CRP; "
                    "ACR_Perc = 100*(DAS28_BL - DAS28_CRP)/DAS28_BL",
            virtual_population=f"{vpop} (each patient = one parameterized sim)",
            trial_timeline={
                "baseline_captured_day": 199,
                "treatment_start_day": 200,
                "first_line_readout_day": 284,
                "second_line_readout_day": 600,
                "note": "a dose name carries StartTime day 200 (suffix _t200); "
                        "append '@DAY' to a dose to switch it in later (sequential "
                        "second-line therapy)",
            },
            drug_formulary=_formulary(),
            calibrated_reference_arms=calibrated_arms,
            held_out="the second-line response rate (MTX_NonResp==1 subgroup, "
                     "day 600) is NOT given - predict it by choosing the therapy",
            edit_spec_help={
                "first_line": "[dose_name, ...] started at day 200 (e.g. "
                              "['MTX_15mg_Q1W_SC_t200'])",
                "second_line": "[dose_name, ...] switched in for MTX-inadequate "
                               "responders (e.g. ['TCZ8mgkg_Q4W_IV_t200'])",
                "switch_day": "day the second_line starts (e.g. 285, just after the "
                              "day-284 first-line readout); omit for concurrent",
            },
        )

    # -- act (run the protocol + read the model's response flags) -------- #
    def run_trial(args: dict, session) -> ToolResult:
        first = args.get("first_line") or []
        second = args.get("second_line") or []
        switch = args.get("switch_day")
        if not first and not second:
            return ToolResult.error(
                "no doses given - set first_line (and optionally second_line + "
                "switch_day). Call ra_inspect for the dose names.")
        spec = osp_ra_trial.build_dose_spec(
            first, second, float(switch) if switch is not None else None)
        r = sb.run_vpop(vpop, dose=spec, stop_time=stop_time,
                        baseline_day=baseline_day, readout_day=readout_day,
                        limit=limit)
        summary = osp_ra_trial.summarize_run(r)
        fl, sl = summary["first_line"], summary["second_line"]

        hist = session.get("ra_history") or []
        hist.append({"protocol": spec, "first_line": fl, "second_line": sl})
        session.put("ra_history", hist)
        session.put("ra_last_second_line", sl)
        session.put("ra_last_protocol", spec)

        ml = (r.get("matlab_log") or "").strip()
        warn = None
        if sl.get("n_MTX_IR", 0) == 0:
            warn = ("no MTX non-responders flagged - the second-line arm is empty. "
                    "Either no first-line MTX was given, or the sim ended before "
                    "day 600. Give MTX first-line and keep --stop-time past 600.")
        return ToolResult.success(
            f"protocol '{spec}' run on {fl.get('n')} patients: "
            f"first-line ACR20 {fl.get('ACR20')}%, remission {fl.get('remission')}%; "
            f"second-line (MTX-IR n={sl.get('n_MTX_IR')}) ACR20 {sl.get('ACR20')}%, "
            f"ACR50 {sl.get('ACR50')}%, ACR70 {sl.get('ACR70')}%"
            + (f"  [WARN] {warn}" if warn else ""),
            protocol=spec,
            first_line=fl,
            second_line=sl,
            das28=summary["das28"],
            warning=warn,
            matlab_tail=ml[-400:] if ml else "",
            iteration=len(hist),
        )

    # -- evaluate (commit the chosen protocol as the answer to be scored) ---- #
    def finalize(args: dict, session) -> ToolResult:
        first = args.get("first_line") or []
        second = args.get("second_line") or []
        switch = args.get("switch_day")
        spec = osp_ra_trial.build_dose_spec(
            first, second, float(switch) if switch is not None else None)
        hist = session.get("ra_history") or []
        match = next((h for h in reversed(hist) if h.get("protocol") == spec), None)
        if match is None:
            return ToolResult.error(
                f"protocol '{spec}' has not been run - call ra_run_trial with it "
                "first, then finalize the exact same protocol.")
        sl = match["second_line"] or {}
        if not sl.get("n_MTX_IR"):
            return ToolResult.error(
                f"protocol '{spec}' has an empty MTX-IR arm - finalize a protocol "
                "that gives MTX first-line and a second-line drug switched in, so "
                "there is a second-line response to report.")
        session.put("ra_final", match)
        return ToolResult.success(
            f"committed final answer: protocol '{spec}', predicted second-line "
            f"(MTX-IR n={sl.get('n_MTX_IR')}) ACR20 {sl.get('ACR20')}% "
            f"ACR50 {sl.get('ACR50')}% ACR70 {sl.get('ACR70')}%",
            protocol=spec, second_line=sl)

    registry.register(Tool(
        name="ra_inspect",
        description=(
            "OBSERVE the RA virtual-trial task: the disease and DAS28/ACR readout, "
            "the trial timeline (baseline day 199, first-line readout day 284, "
            "second-line readout day 600), the drug formulary (MTX + the biologics "
            "adalimumab/tocilizumab/secukinumab/anakinra, with mechanism and dose "
            "names), the calibrated reference arms with their KNOWN response rates "
            "(run these to confirm your harness reproduces them), and the held-out "
            "objective. Call this first."),
        input_schema={"type": "object", "properties": {}},
        handler=inspect, phase="observe"))

    registry.register(Tool(
        name="ra_run_trial",
        description=(
            "ACT: run a treatment protocol across the virtual population and read "
            "the model's own response flags. first_line doses start at day 200; "
            "second_line doses are switched in at switch_day (use a day just after "
            "the day-284 first-line readout, e.g. 285, so the first-line "
            "classification stays pure and the second-line therapy then acts on the "
            "MTX-inadequate responders through the day-600 readout). Returns the "
            "first-line ACR20/50/70/remission (all patients, day 284) and the "
            "second-line rates among MTX non-responders (day 600). A second-line "
            "given concurrently from day 200 conflates the arms - switch it in "
            "later."),
        input_schema={"type": "object", "properties": {
            "first_line": {"type": "array", "items": {"type": "string"},
                           "description": "dose names started at day 200"},
            "second_line": {"type": "array", "items": {"type": "string"},
                            "description": "dose names switched in for MTX-IR"},
            "switch_day": {"type": "number",
                           "description": "day the second_line starts (e.g. 285)"}}},
        handler=run_trial, phase="act"))

    registry.register(Tool(
        name="ra_finalize",
        description=(
            "COMMIT your final answer: the protocol you are predicting with (same "
            "first_line / second_line / switch_day you passed to ra_run_trial). Call "
            "this once, on the protocol you actually recommend, BEFORE finishing - "
            "it is the run that gets scored, so do not leave your answer as whatever "
            "control or dose-response check you happened to run last. The protocol "
            "must already have been run and have a non-empty MTX-IR arm."),
        input_schema={"type": "object", "properties": {
            "first_line": {"type": "array", "items": {"type": "string"}},
            "second_line": {"type": "array", "items": {"type": "string"}},
            "switch_day": {"type": "number"}}},
        handler=finalize, phase="evaluate"))
