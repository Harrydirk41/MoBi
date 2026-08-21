"""Model-agnostic loop tools for the virtual-trial task, driven by a QSPTaskConfig.

The agent chooses a treatment PROTOCOL (which drugs, and - for a second-line agent -
when to switch it in) and runs it across the virtual population; the model's own
events compute the clinical response. All vocabulary (disease, drug formulary,
timeline, readout column roles) is read off the config, so pointing at a different
QSP model needs only a different projects/<name>/tasks.json.

  * ``trial_inspect`` (observe) - the disease/readout, timeline, drug formulary, the
    calibrated reference arms, and the held-out objective.
  * ``trial_run``     (act)     - apply a protocol, run the Vpop, return the model's
    first-line and second-line response rates.
  * ``trial_finalize``(evaluate)- commit the chosen protocol to be scored.
"""

from __future__ import annotations

from ..engines.simbiology import SimBiologyEngine
from ..engines import qsp_tasks
from ..engines.qsp_config import QSPTaskConfig
from .registry import Tool, ToolRegistry, ToolResult


def register_qsp_trial_loop_tools(registry: ToolRegistry, config, ctx: dict) -> None:
    """ctx: {cfg, sb, vpop, limit, stop_time, baseline_day, readout_day,
    objective?, calibrated_arms?}."""
    cfg: QSPTaskConfig = ctx["cfg"]
    sb: SimBiologyEngine = ctx["sb"]
    vpop: str = ctx["vpop"]
    calibrated_arms: list = ctx.get("calibrated_arms") or cfg.calibrated_arms
    objective: str = ctx.get("objective") or cfg.trial_objective \
        or "Predict the second-line response."
    limit: int = int(ctx.get("limit") or 50)
    stop_time: float = float(ctx.get("stop_time") or 700.0)
    baseline_day: float = float(ctx.get("baseline_day")
                                or cfg.timeline.get("baseline_day", 200.0))
    readout_day: float = float(ctx.get("readout_day")
                               or cfg.timeline.get("first_line_readout_day", 284.0))

    def _formulary() -> list[dict]:
        return [{"code": code, "drug": d.get("drug"), "modality": d.get("modality"),
                 "mechanism": d.get("mechanism"), "dose_names": d.get("doses")}
                for code, d in cfg.drugs.items()]

    def inspect(args: dict, session) -> ToolResult:
        return ToolResult.success(
            f"{cfg.name} virtual-trial task: choose a treatment protocol and run it "
            "across the virtual population; the model's events compute the clinical "
            "response.",
            objective=objective,
            disease=cfg.disease,
            virtual_population=f"{vpop} (each patient = one parameterized sim)",
            trial_timeline={**{k: v for k, v in cfg.timeline.items()},
                            "note": "the subgroup classification and the first-line "
                                    "readout are the model's own events - work out how "
                                    "a protocol interacts with them."},
            drug_formulary=_formulary(),
            calibrated_reference_arms=calibrated_arms,
            held_out="the second-line response (the subgroup at the second readout) is "
                     "NOT given - you determine the therapy and predict it",
            edit_spec_help={
                "first_line": "[dose_name, ...] applied from treatment start",
                "second_line": "[dose_name, ...] applied to the same patients",
                "switch_day": "OPTIONAL day to override the second_line start (days)",
                "dose_scale": "OPTIONAL multiplier on the second_line dose amount "
                              "(1.0 = labeled dose)"})

    def run_trial(args: dict, session) -> ToolResult:
        first = args.get("first_line") or []
        second = args.get("second_line") or []
        switch = args.get("switch_day")
        scale = args.get("dose_scale")
        if not first and not second:
            return ToolResult.error(
                "no doses given - set first_line (and optionally second_line, "
                "switch_day, dose_scale). Call trial_inspect for the dose names.")
        spec = qsp_tasks.build_dose_spec(
            first, second, float(switch) if switch is not None else None,
            float(scale) if scale is not None else None)
        r = sb.run_vpop(vpop, dose=spec, stop_time=stop_time,
                        baseline_day=baseline_day, readout_day=readout_day, limit=limit)
        summary = cfg.summarize_run(r)
        fl, sl = summary["first_line"], summary["second_line"]

        hist = session.get("trial_history") or []
        hist.append({"protocol": spec, "first_line": fl, "second_line": sl})
        session.put("trial_history", hist)
        session.put("trial_last_second_line", sl)

        ml = (r.get("matlab_log") or "").strip()
        warn = None
        if sl.get("n_subgroup", 0) == 0:
            warn = ("no subgroup patients flagged - the second-line arm is empty. "
                    "Either the first-line arm resolved everyone, or the sim ended "
                    "before the second readout. Keep --stop-time past it.")
        return ToolResult.success(
            f"protocol '{spec}' run on {fl.get('n')} patients: "
            f"first-line {fl}; second-line (subgroup n={sl.get('n_subgroup')}) {sl}"
            + (f"  [WARN] {warn}" if warn else ""),
            protocol=spec, first_line=fl, second_line=sl, das28=summary["das28"],
            warning=warn, matlab_tail=ml[-400:] if ml else "", iteration=len(hist))

    def finalize(args: dict, session) -> ToolResult:
        first = args.get("first_line") or []
        second = args.get("second_line") or []
        switch = args.get("switch_day")
        scale = args.get("dose_scale")
        spec = qsp_tasks.build_dose_spec(
            first, second, float(switch) if switch is not None else None,
            float(scale) if scale is not None else None)
        hist = session.get("trial_history") or []
        match = next((h for h in reversed(hist) if h.get("protocol") == spec), None)
        if match is None:
            return ToolResult.error(
                f"protocol '{spec}' has not been run - call trial_run with it first, "
                "then finalize the exact same protocol.")
        sl = match["second_line"] or {}
        if not sl.get("n_subgroup"):
            return ToolResult.error(
                f"protocol '{spec}' has an empty subgroup arm - finalize a protocol "
                "that produces a second-line response to report.")
        session.put("trial_final", match)
        return ToolResult.success(
            f"committed final answer: protocol '{spec}', predicted second-line "
            f"(subgroup n={sl.get('n_subgroup')}) {sl}",
            protocol=spec, second_line=sl)

    _sched = {"type": "object", "properties": {
        "first_line": {"type": "array", "items": {"type": "string"}},
        "second_line": {"type": "array", "items": {"type": "string"}},
        "switch_day": {"type": "number"}, "dose_scale": {"type": "number"}}}

    registry.register(Tool(
        name="trial_inspect",
        description=("OBSERVE the virtual-trial task: the disease and readout, the "
                     "trial timeline, the drug formulary with mechanisms and dose "
                     "names, the calibrated reference arms with their KNOWN response "
                     "rates (run these to confirm your harness reproduces them), and "
                     "the held-out objective. Call this first."),
        input_schema={"type": "object", "properties": {}},
        handler=inspect, phase="observe"))
    registry.register(Tool(
        name="trial_run",
        description=("ACT: run a treatment protocol across the virtual population and "
                     "read the model's own response flags (first-line, all patients; "
                     "second-line, the subgroup). first_line/second_line take dose "
                     "names from the formulary; switch_day and dose_scale optionally "
                     "retime and rescale the second_line."),
        input_schema=_sched, handler=run_trial, phase="act"))
    registry.register(Tool(
        name="trial_finalize",
        description=("COMMIT your final answer: the protocol you are predicting with "
                     "(same fields you passed to trial_run). Call once, on the "
                     "protocol you actually recommend, BEFORE finishing - it is the "
                     "run that gets scored. The protocol must already have been run "
                     "and have a non-empty subgroup arm."),
        input_schema=_sched, handler=finalize, phase="evaluate"))
