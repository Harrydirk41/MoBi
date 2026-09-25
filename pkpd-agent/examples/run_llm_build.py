r"""[LEGACY - superseded by the Claude Code sandbox flow; kept as a scripted baseline.]

The primary way to build a model is now a CONTINUOUS Claude Code agent in a sealed sandbox,
not this fixed decision loop. That path removes the scaffolding this script imposes (the
agent decides its own order of operations) and walls the answers off from a filesystem-
capable agent. See SANDBOX_CC.md; in short:

    pkpd-bench seal   --model <Name> --out ../sandbox
    pkpd-bench verify --sandbox ../sandbox/<Name>
    # point a coding agent at ../sandbox/<Name>/workspace and have it follow PROMPT.md
    pkpd-bench judge  --sandbox ../sandbox/<Name>   # held-out grade

This script remains as a deterministic, no-CC baseline (and for the report/scoreboard
wiring). Prefer the sandbox flow for new runs.

---

Real-scenario PBPK modeling loop: the LLM decides, the optimizer fits.

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
import sys

# Force UTF-8 console output: the model's reasoning and PK-Sim units contain
# non-ASCII characters (µmol/l, alpha/beta metabolites), which crash a Windows
# console whose stdout defaults to ascii/cp1252 (UnicodeEncodeError). Reconfigure
# to UTF-8 with replacement so a stray glyph never aborts a multi-hour run.
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(encoding="utf-8", errors="replace")
    except (AttributeError, ValueError):
        pass

from pkpd_agent.config import AgentConfig
from pkpd_agent.engines.osp_cli import OSPCli
from pkpd_agent.llm import LLMPolicy
from pkpd_agent.loop import DecisionLoop
from pkpd_agent.state import Decision, Finish, ModelingSession, Observation
from pkpd_agent.tools.registry import ToolRegistry
from pkpd_agent.tools.osp_loop_tools import register_osp_loop_tools


def _system_prompt(target: float, self_extract: bool = False, web: bool = False) -> str:
    source = ("From the handed report and data, and from web search / web fetch "
              "when you need a value or the mechanism they do not give"
              if web else "ONLY from the handed report and data - you have no "
              "web access")
    web_rule = (
        "YOU HAVE WEB SEARCH (web_search / web_fetch) - USE IT, do not answer from "
        "memory. For any numeric prior NOT in the handed report - especially "
        "lipophilicity/logP, plasma clearance (or hepatic extraction), blood:plasma "
        "ratio, fraction metabolized, and in-vitro CLint/Km - actually CALL "
        "web_search, take the value from a real source, and record it via "
        "osp_record_givens with that citation (provenance 'given'). If the clearing "
        "enzyme or a transporter role is not explicit in the report, confirm it by "
        "web_search too. NEVER state a number or a mechanism as fact from your own "
        "memory when you could look it up - and never describe a search you did not "
        "run. If a published PBPK model exists for this compound, fetch it and build "
        "on its structure, citing it.\n\n") if web else ""
    step0 = (
        "0. SELF-EXTRACT the givens FIRST. osp_inspect will show NO pre-digested "
        "literature_physicochemical - you build the context yourself from the raw "
        "materials you were handed. Call osp_read_report (background, the "
        "physicochemical/in-vitro table, and the metabolizing enzyme named in the "
        "prose), and osp_list_studies / osp_read_study for the raw curves. Then "
        "call osp_record_givens with the values you extracted (lipophilicity/logP, "
        "fraction unbound, pKa, solubility, molecular weight, in-vitro Km/Vmax, and "
        "the clearing enzyme), each with its source citation and provenance "
        "('given' = read from the report, 'judged' = inferred). Your recorded "
        "givens drive the fix-vs-fit split, so do this BEFORE osp_options / "
        f"osp_optimize. Extract {source}; the report's chosen methods and fitted "
        "values are redacted.\n") if self_extract else ""
    return (
        "You are a PBPK modeler using the OSP PK-Sim engine with a numerical "
        "optimizer. The whole-body physiological structure is fixed; you decide "
        "the DRUG model and let the optimizer fit the numbers.\n\n"
        "NARRATE YOUR REASONING: before every tool call, write one or two plain "
        "sentences saying what you are about to do and WHY (the hypothesis or the "
        "evidence driving it). This text is shown to the user as your reasoning, so "
        "never emit a tool call with no accompanying sentence.\n\n"
        "ORIENT BEFORE YOU BUILD. Do not rush to the sweep. First, in your own "
        "words: (a) restate the task and the key constraints you must honour (dose "
        "range, routes, the known biology); (b) SURVEY the reference_library index "
        "from osp_inspect, then actually OPEN the 2-3 analogues closest in "
        "chemistry/clearance with osp_read_reference and READ each one's report "
        "(its full modeling write-up, not just the summary) - the index alone is "
        "not enough; say what distribution method and which processes each used and "
        "which you will lean on and why; (c) list the parameters you can "
        "FIX from givens vs must FIT. Show this reasoning before proposing a "
        "structure - it is the understanding step, and skipping it is the main "
        "failure mode.\n\n"
        "READ THE DATA before you fit: call osp_read_study on at least the IV curves "
        "(and one PO curve if present) to SEE the profile shape - how many "
        "distribution phases, the terminal slope, any dose-nonlinearity. That shape "
        "decides how many processes/compartments the model needs. Do not go from the "
        "study LIST straight to fitting.\n\n"
        + web_rule +
        "Each round:\n"
        + step0 +
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
        "   If osp_inspect returns a reference_library (LIBRARY-ASSISTED mode), USE "
        "it: it holds finished models of OTHER, analogous compounds (leave-one-out). "
        "Pick your STRUCTURE - distribution/permeability method and which "
        "enzyme/transporter processes - by analogy to the most chemically/"
        "mechanistically similar reference model, then FIT the parameters to THIS "
        "compound's own data (do not copy their numbers blindly).\n"
        "2. DETERMINE the model:\n"
        "   - Structure: pick the distribution & permeability calculation methods "
        "and which processes are active, guided by the known biology (e.g. "
        "CYP3A4 metabolism, negligible renal clearance, no active transport). If "
        "the distribution / Vd is off, choose the method with osp_sweep_methods, but "
        "REASON FIRST and screen a SHORTLIST rather than brute-forcing all five. "
        "From the compound's chemistry (a lipophilic neutral, a strong base, an "
        "acid distribute differently) and the reference-library analogues (a close "
        "analogue's distribution method is a strong prior), pick the 2-3 plausible "
        "partition methods and pass them as partition_methods=[...]; the sweep "
        "FAST-SCREENS only those, re-fits your physchem under each (so no method is "
        "judged at frozen physchem), and refines the best by the DATA - no answer "
        "is used. Sweep all five only when you have no basis to narrow. If the "
        "result says NO screened method fit adequately (screen_adequate=false), the "
        "lever is elsewhere: free a measured physchem WITHIN its measured range, or "
        "revisit the mechanism (add/disable a process), then re-sweep the shortlist "
        "- do NOT accept a poor winner. Otherwise, once it has chosen the method you "
        "REFINE with osp_optimize - you do NOT re-sweep the same grid (a repeat "
        "identical sweep is skipped anyway). "
        "It also AUTO-WIDENS a free (non-given) physchem parameter that "
        "rails to your bound: if the input gave NO measured logP, the effective "
        "Lipophilicity is a fit-target, not the measured logP - so do NOT cap it "
        "at a measured value you were never given (a hydrophilic drug can still "
        "need a high effective lipophilicity for distribution). NEVER distort a "
        "GIVEN measured physchem parameter (or drop GFR fraction well below 1) to "
        "force Vd before you have swept the methods; a wrong Vd is usually a wrong "
        "distribution METHOD, not a wrong measured value. You "
        "may disable a process (processes:{Molecule:false}) or ADD a mechanism "
        "(add_processes:[{type,molecule,parameters}]) to an expressed molecule if "
        "the data show a systematic misfit a single enzyme cannot explain.\n"
        "   - ABSORPTION (oral drugs): the Weibull 'Dissolution time (50% dissolved)' "
        "and 'Dissolution shape' are placeholders you must FIT for any ORAL "
        "formulation - they set the absorption rate/shape and are not measured. With "
        "more than one formulation (tablet vs capsule, IR vs ER) fit each separately "
        "with the qualified name 'Dissolution time (50% dissolved)@<FormulationName>'. "
        "For IV-only drugs, ignore them.\n"
        "   - Decide which parameters to ESTIMATE vs FIX, following the "
        "MEASURED-QUANTITY PRINCIPLE and a STAGED strategy:\n"
        "     * Measured physical CONSTANTS (MW, pKa, reference pH) are NEVER "
        "estimated - the optimizer will reject them. Always fix them.\n"
        "     * If the INPUT GIVES a physchem value (lipophilicity/logP, pKa), that "
        "is a MEASUREMENT - FIX it, never put it in 'estimate' (a given lipophilicity "
        "is moved to 'fix' automatically). Only fit lipophilicity as an effective "
        "parameter when the input gives none.\n"
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
        "   - STAGED FIT: if osp_inspect.method_guidance gives 'route_staging' "
        "(both IV and PO data present), fit in STAGES - first on the IV data set "
        "the distribution method and fit distribution + systemic clearance; then "
        "on the PO data HOLD those fixed and fit ONLY absorption parameters. A "
        "single joint fit lets absorption and clearance compensate and hide a "
        "structural gap. Respect method_guidance.saturable_identifiability: do "
        "not fit Michaelis-Menten/saturable constants without >=2 dose levels.\n"
        "4. CHECK the result:\n"
        "   - GMFE and per-route bias - is the fit acceptable?\n"
        "   - pk_parameter_gmfe (Cmax/tmax/AUC/t1/2 fold-errors): a good pointwise "
        "GMFE with a poor Cmax or t1/2 fold means the PEAK or TERMINAL SLOPE is "
        "off - a structural clue (absorption/distribution vs elimination), not a "
        "numbers-only issue.\n"
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
    ap.add_argument("--library", default=None, choices=["all", "same-type"],
                    help="LIBRARY-ASSISTED mode: inject the OTHER compounds' finished models "
                         "(leave-one-out) as a reference library. 'same-type' keeps analogues "
                         "sharing an enzyme/transporter with the target.")
    ap.add_argument("--self-extract", action="store_true",
                    help="agent builds its own givens from the handed report+data "
                         "(pbpk-realworld tree) instead of a pre-digested list")
    ap.add_argument("--realworld", default=None,
                    help="pbpk-realworld root (default: ../pbpk-realworld next to the library)")
    ap.add_argument("--library-full", action="store_true",
                    help="include the reference models' fitted VALUES (default: structure + "
                         "parameter names only, so the agent still fits the numbers itself)")
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

    # HELD-OUT split: the agent fits on the BUILDING studies only; the VERIFICATION
    # studies are held out and used only to grade prediction (report.assemble). Compute
    # the split from the blanked snapshot (its OutputMappings/classifications, which the
    # agent legitimately has) and hand the tools a building-only view of the data.
    from pkpd_agent.engines import osp_split as _osp_split
    with open(args.snapshot, encoding="utf-8") as _sf:
        _split = _osp_split.split_studies(json.load(_sf), observed)
    build_set = set(_split.get("building") or [])
    if _split.get("verification"):
        build_obs = [o for o in observed if o["dataset"] in build_set]
        inp_agent = dict(inp)
        inp_agent["given_data"] = {**inp["given_data"], "clinical_observed_data": build_obs}
        print(f"== held-out split [{_split['method']}]: agent fits {len(build_obs)} "
              f"building datasets; {len(_split['verification'])} verification datasets held "
              f"out ({len(_split['held_out_studies'])} studies) for grading ==")
    else:
        build_obs, inp_agent = observed, inp
        print(f"== held-out split: none ({_split.get('method')}) - graded on all data ==")

    # LIBRARY-ASSISTED mode: inject the OTHER compounds' finished models (leave-one-out).
    if args.library:
        from pkpd_agent.engines import reference_library as _rl
        inp_agent = dict(inp_agent)
        lib = _rl.library_for_snapshot(args.snapshot,
                                       same_type=(args.library == "same-type"),
                                       full=args.library_full)
        inp_agent["reference_library"] = lib
        print(f"== library-assisted [{lib['mode']}]: {len(lib['models'])} reference "
              f"model(s) available to the agent (leave-one-out) ==")

    # SELF-EXTRACT mode: strip the pre-digested givens and hand the agent the raw
    # report + data so it builds its own context data lake (osp_read_report /
    # osp_record_givens). Confined to the handed materials - no web = no leak.
    ctx_report = ctx_data = None
    if args.self_extract:
        comp = os.path.basename(os.path.dirname(os.path.dirname(args.snapshot)))
        rw = args.realworld or os.path.join(os.path.dirname(args.snapshot),
                                            "..", "..", "..", "pbpk-realworld")
        tdir = os.path.join(rw, comp, "test")
        import glob as _glob
        reps = _glob.glob(os.path.join(tdir, "report", "*.md"))
        ddir = os.path.join(tdir, "data")
        if reps and os.path.isdir(ddir):
            ctx_report, ctx_data = reps[0], ddir
            inp_agent = dict(inp_agent)
            inp_agent["given_data"] = {k: v for k, v in inp_agent["given_data"].items()
                                       if k != "literature_physicochemical"}
            print(f"== self-extract: agent builds givens from {comp}/test/"
                  "{report,data} (pre-digested list withheld) ==")
        else:
            print(f"== self-extract requested but no realworld tree for {comp}; "
                  "falling back to pre-digested givens ==")

    registry = ToolRegistry()
    register_osp_loop_tools(registry, cfg, {
        "cli": cli, "snapshot_path": args.snapshot,
        "observed": build_obs, "input": inp_agent,
        "self_extract": bool(ctx_report)})
    if ctx_report:
        from pkpd_agent.tools.context_tools import register_context_tools
        register_context_tools(registry, cfg, {
            "report_path": ctx_report, "data_dir": ctx_data, "input": inp_agent})

    goal = (f"{inp.get('objective','Build the PBPK model.')}\n\n"
            f"Target: overall GMFE <= {args.target}. Start by calling osp_inspect, "
            "then determine the model and call osp_optimize.")
    policy = LLMPolicy(cfg, registry, _system_prompt(args.target, self_extract=bool(ctx_report)))
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
