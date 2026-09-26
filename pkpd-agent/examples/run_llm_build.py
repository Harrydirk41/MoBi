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


def _system_prompt(target: float, self_extract: bool = False, web: bool = False,
                   minimal: bool = False) -> str:
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
    # Load the editable workflow policy (the "skill.md"): strip the design-note
    # HTML comments, then substitute the conditional/dynamic parts. The playbook is
    # the single source of truth; a missing file fails loudly rather than silently
    # running a stale copy.
    try:
        import re as _re
        import pkpd_agent as _pkg
        from pathlib import Path as _Path
        fname = "modeling_playbook_minimal.md" if minimal else "modeling_playbook.md"
        pb = (_Path(_pkg.__file__).resolve().parent / "prompts"
              / fname).read_text(encoding="utf-8")
        pb = _re.sub(r"(?s)<!--.*?-->", "", pb).strip()
        return (pb.replace("{{WEB_RULE}}", web_rule)
                  .replace("{{STEP0}}", step0)
                  .replace("{{TARGET}}", f"{target}"))
    except (OSError, ImportError) as e:
        raise RuntimeError(
            "modeling playbook not found (pkpd_agent/prompts/modeling_playbook.md) — it is the agent's workflow "
            "policy and must ship with the package.") from e


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
