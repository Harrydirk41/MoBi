r"""Real-scenario PBPK modeling loop: the LLM decides, the optimizer fits.

Each round:
  1. osp_inspect        - read objective, known biology, priors, data, model
  2. DETERMINE the model - structure (methods/processes) + which parameters to
                           ESTIMATE (with bounds) vs FIX at literature values
  3. osp_optimize       - a numerical optimizer fits the chosen parameters
  4. CHECK              - GMFE, per-route bias, params_at_bound (identifiability),
                           systematic misfit (wrong structure?)
  5. REASON -> next round, or finish

This is the division of labor the project is built on: the LLM does the
judgment (structure, what to estimate), a numerical optimizer does the fitting.

    set ANTHROPIC_API_KEY=...
    set PKPD_PKSIM_CLI=C:\Program Files\Open Systems Pharmacology\PK-Sim 12.3\PKSim.CLI.exe

    python -m examples.run_llm_build ^
        --snapshot ..\OSP-PBPK-Model-Library\Alfentanil\benchmark\Alfentanil-Model.blanked.json ^
        --input    ..\OSP-PBPK-Model-Library\Alfentanil\json_input\Alfentanil-Model.input.json ^
        --target 1.6 --max-steps 6
"""

from __future__ import annotations

import argparse
import json
import os

from pkpd_agent.config import AgentConfig
from pkpd_agent.engines.osp_cli import OSPCli
from pkpd_agent.llm import LLMPolicy
from pkpd_agent.loop import DecisionLoop
from pkpd_agent.state import Decision, Finish, ModelingSession, Observation
from pkpd_agent.tools.registry import ToolRegistry
from pkpd_agent.tools.osp_loop_tools import register_osp_loop_tools


def _system_prompt(target: float) -> str:
    return (
        "You are a PBPK modeler using the OSP PK-Sim engine with a numerical "
        "optimizer. The whole-body physiological structure is fixed; you decide "
        "the DRUG model and let the optimizer fit the numbers.\n\n"
        "Each round:\n"
        "1. Call osp_inspect (task: objective, known biology, priors, data) and "
        "osp_options (the authoritative ACTION SPACE: every editable parameter "
        "with its role/range, the legal calculation methods, the molecules the "
        "model expresses, and the mechanisms you may ADD). Use osp_options as the "
        "source of truth for names and legal choices - do not assume OSP details. "
        "Each parameter carries a value_status: 'given' = a measured input, TRUST "
        "its value; 'placeholder' = a naive benchmark default that is NOT a "
        "measurement - the shown number is meaningless, you must DETERMINE the "
        "parameter (from the data, or from established physchem knowledge if you "
        "have it); 'structural' = a real constant present as-is (e.g. molecular "
        "weight). Never trust a placeholder's starting value as if it were data.\n"
        "2. DETERMINE the model:\n"
        "   - Structure: pick the distribution & permeability calculation methods "
        "and which processes are active, guided by the known biology (e.g. "
        "CYP3A4 metabolism, negligible renal clearance, no active transport). If "
        "the distribution / Vd is off, call osp_sweep_methods FIRST: it tries EVERY "
        "partition x permeability method and re-fits your physchem under each, then "
        "adopts the best - this is the reliable, deterministic way to choose the "
        "method. NEVER distort a measured physchem parameter (e.g. drive "
        "lipophilicity far from its literature logP, or drop GFR fraction well "
        "below 1) to force Vd before you have swept the methods; a wrong Vd is "
        "usually a wrong distribution METHOD, not a wrong measured value. You "
        "may disable a process (processes:{Molecule:false}) or ADD a mechanism "
        "(add_processes:[{type,molecule,parameters}]) to an expressed molecule if "
        "the data show a systematic misfit a single enzyme cannot explain.\n"
        "   - Decide which parameters to ESTIMATE vs FIX, following the "
        "MEASURED-QUANTITY PRINCIPLE and a STAGED strategy:\n"
        "     * Measured physical CONSTANTS (MW, pKa, reference pH) are NEVER "
        "estimated - the optimizer will reject them. Always fix them.\n"
        "     * USE WHAT YOU KNOW. If a parameter has a well-established value - "
        "from the given literature OR your own pharmacological knowledge (a "
        "reported logP/logD for lipophilicity, a measured pKa, a known fraction "
        "unbound) - FIX it to that value, or bound it TIGHTLY around it. Do NOT "
        "float a parameter you actually know across a wide range 'to let the "
        "optimizer decide' - a wide bound on a known quantity invites the "
        "optimizer to park it at a physically wrong value that merely compensates "
        "for other errors (e.g. lipophilicity driven to ~0 for a drug whose logP "
        "is ~2). The optimizer has no physical prior; you do - encode it.\n"
        "     * STAGE 1: estimate only the genuinely uncertain, NOT reliably "
        "measurable parameters (intrinsic clearance, permeabilities) - the minimal "
        "identifiable set (2-3). Fix everything else, including MEASURED-BUT-"
        "REFINABLE quantities (fraction unbound, solubility) and any parameter you "
        "have a literature/textbook value for, at that value.\n"
        "     * STAGE 2 (only if Stage 1 leaves a SYSTEMATIC misfit that a "
        "measured quantity could legitimately explain): free that measured "
        "quantity and re-optimize. It will be automatically constrained to its "
        "MEASURED RANGE - do not try to push it outside. If a good fit needs a "
        "measured value outside its measured range, your STRUCTURE is wrong, not "
        "the measurement.\n"
        "   - Prefer the simplest model that fits: do not free a measured quantity "
        "unless the data demand it. Use the per-parameter 'sensitivity' returned "
        "by osp_optimize - a low-sensitivity parameter is weakly identified, so "
        "fix it rather than reporting an uncertain fitted value.\n"
        "   - Give each estimated parameter a plausible [lo, hi] bound.\n"
        "3. Call osp_optimize with {estimate:{param:[lo,hi]}, fix:{...}, "
        "structure:{...}}. The optimizer fits the estimated parameters.\n"
        "4. CHECK the result:\n"
        "   - GMFE and per-route bias - is the fit acceptable?\n"
        "   - 'recommendations' - the optimizer returns concrete identifiability "
        "ACTIONS. ACT ON THEM: if a parameter has near-zero sensitivity, the data "
        "cannot pin it and its value is an artifact - FIX it to a known/literature "
        "value (yours or the given anchor) and refit the smaller set. If two "
        "parameters are collinear and ONE has a literature anchor, fix that one "
        "and estimate the other. If SEVERAL clearances are collinear and NONE has "
        "an anchor (per-enzyme CLspec/kcat that plasma data cannot split), do NOT "
        "fix one at a guess - that corrupts the identifiable TOTAL. Instead use "
        "link_scale=[{members:[..], bounds:[lo,hi]}] to fit their shared magnitude "
        "(the total clearance) while holding their ratio; the split stays "
        "unresolved but the total - which the data DO constrain - is recovered. A "
        "smaller, well-conditioned fit of the CORRECT structure beats a larger fit "
        "that trades parameters off.\n"
        "   - params_at_bound - a parameter pinned to a bound is unidentifiable "
        "or signals a WRONG STRUCTURE. Widen the bound, fix it, or change the "
        "structure.\n"
        "   - A systematic per-route misfit that fitting cannot remove (e.g. oral "
        "always off) points to a missing/incorrect MECHANISM - change the "
        "structure, not just the numbers.\n"
        "   - STRUCTURE by biology, not by a marginal GMFE gain. Do not drop a "
        "named elimination pathway (e.g. lump UGT + CYP into one term) just "
        "because a simpler model scores slightly better on this fit - a marginal "
        "improvement often means the richer model was fit poorly (too many "
        "collinear free parameters), not that the pathway is absent. Fix the "
        "conditioning (per the recommendations) and refit the biologically "
        "correct structure before abandoning it.\n"
        f"5. REASON and iterate until overall GMFE <= {target} with plausible, "
        "identifiable parameters, then finish with a summary: the structure, "
        "which parameters you estimated vs fixed, their fitted values, and the "
        "GMFE.\n\n"
        "ENGINE ERRORS ARE HARD FAILURES, NOT RESULTS. If a tool returns an error "
        "(e.g. 'PK-Sim run failed', 'optimization failed', 'no .pksim5 produced') "
        "the model did NOT run - do not accept it and do not route around it by "
        "abandoning optimization. Read the error, REVISE the offending edit (a "
        "wrong parameter/method/process name, an illegal structure), and RE-RUN. "
        "Your deliverable is a model whose parameters were fitted by a SUCCESSFUL "
        "osp_optimize; do not finish on an un-optimized or failed run.\n\n"
        "Do NOT hand-tune parameter values - that is the optimizer's job. Your "
        "job is the modeling decisions: structure, and what to estimate vs fix."
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--snapshot", required=True)
    ap.add_argument("--input", required=True)
    ap.add_argument("--pksim", default=None)
    ap.add_argument("--target", type=float, default=1.6)
    ap.add_argument("--max-steps", type=int, default=6)
    ap.add_argument("--model", default=None)
    ap.add_argument("--effort", default=None)
    ap.add_argument("--report", default=None,
                    help="write an evaluation report here (.html; a .pdf is also "
                         "written if matplotlib is installed)")
    ap.add_argument("--reference", default=None,
                    help="reference snapshot for the ground-truth comparison in "
                         "the report (e.g. the original json/<Compound>-Model.json)")
    ap.add_argument("--answer-edits", default=None,
                    help="answer_key edit spec for the parameter comparison")
    args = ap.parse_args()

    cfg = AgentConfig(mock=False, max_steps=args.max_steps)
    if args.model:
        cfg.model = args.model
    if args.effort:
        cfg.effort = args.effort
    if not cfg.anthropic_key_present():
        print("ANTHROPIC_API_KEY not set.")
        return

    cli = OSPCli(pksim_cli_path=args.pksim or cfg.pksim_cli_path,
                 timeout_s=cfg.pksim_timeout_s)
    if not cli.pksim_cli_path or not os.path.exists(cli.pksim_cli_path):
        print(f"PKSim.CLI not found at {cli.pksim_cli_path!r}.")
        return

    with open(args.input, encoding="utf-8") as fh:
        inp = json.load(fh)
    observed = inp["given_data"]["clinical_observed_data"]

    registry = ToolRegistry()
    register_osp_loop_tools(registry, cfg, {
        "cli": cli, "snapshot_path": args.snapshot,
        "observed": observed, "input": inp})

    goal = (f"{inp.get('objective','Build the PBPK model.')}\n\n"
            f"Target: overall GMFE <= {args.target}. Start by calling osp_inspect, "
            "then determine the model and call osp_optimize.")
    policy = LLMPolicy(cfg, registry, _system_prompt(args.target))
    loop = DecisionLoop(config=cfg, registry=registry, policy=policy)

    print(f"== LLM PBPK BUILD loop on {os.path.basename(args.snapshot)} "
          f"(target GMFE {args.target}, max {args.max_steps} steps) ==")
    print("(osp_optimize runs the model many times to fit - several minutes per "
          "call; progress streams below)\n")

    def show(ev):
        if isinstance(ev, Decision):
            if ev.text:
                print(f"\n[reason] {ev.text[:2500]}")
            for c in ev.calls:
                a = {k: v for k, v in c.arguments.items()}
                print(f"  -> {c.name} {json.dumps(a, ensure_ascii=False)[:400]}")
                if c.name in ("osp_optimize", "osp_try_model"):
                    print("     ...running PK-Sim, please wait...", flush=True)
        elif isinstance(ev, Observation):
            c = ev.content
            print(f"  <- {ev.tool}: {c.get('message','')}", flush=True)
            if c.get("optimized"):
                print(f"       optimized: {c['optimized']}")
            for r, m in (c.get("by_route") or {}).items():
                print(f"       {r}: GMFE {m.get('gmfe')} bias {m.get('bias')}")
            if c.get("params_at_bound"):
                print(f"       AT BOUND: {c['params_at_bound']}")
            if c.get("parameter_flags"):
                print(f"       PLAUSIBILITY: {c['parameter_flags']}")
        elif isinstance(ev, Finish):
            print(f"\n=== SUMMARY ===\n{ev.text}")

    session = loop.run(goal, ModelingSession(goal=goal), on_event=show)
    print(f"\nbest GMFE reached: {session.get('osp_best_gmfe')}")
    print(f"best model: {json.dumps(session.get('osp_best_edits'), ensure_ascii=False)}")

    if args.report:
        from pkpd_agent import report
        answer = None
        if args.answer_edits and os.path.exists(args.answer_edits):
            with open(args.answer_edits, encoding="utf-8") as fh:
                answer = json.load(fh)
        print("\nassembling report (re-running best + reference models)...", flush=True)
        d = report.assemble(session, cfg, cli, inp, args.snapshot,
                            session.get("osp_best_edits") or {},
                            ref_snapshot_path=args.reference, answer_edits=answer)
        # always write into the case's report/ folder (case = parent of json_input)
        case_dir = os.path.dirname(os.path.dirname(os.path.abspath(args.input)))
        report_dir = os.path.join(case_dir, "report")
        os.makedirs(report_dir, exist_ok=True)
        name = os.path.basename(args.report) or "report.html"
        if not name.endswith(".html"):
            name += ".html"
        html_path = os.path.join(report_dir, name)
        report.write_html(d, html_path)
        print(f"wrote {html_path}")
        json_path = html_path[:-5] + ".json"
        report.write_json(d, json_path)             # machine-readable scoreboard payload
        print(f"wrote {json_path}")
        pdf_path = html_path[:-5] + ".pdf"
        if report.write_pdf(d, pdf_path):
            print(f"wrote {pdf_path}")
        else:
            print("(matplotlib not installed - open the .html and Print -> Save as PDF)")


if __name__ == "__main__":
    main()
