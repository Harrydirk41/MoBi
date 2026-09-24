# PBPK Copilot (local web app)

A browser UI over the PBPK modeling agent. You pick a compound, press **Build
with agent**, and the real agent loop (inspect → sweep → optimize → grade)
streams into the page live. It runs on your machine because it needs PKSim.CLI
and an Anthropic key; the browser only talks to this local server.

## Run

```powershell
cd pkpd-agent
pip install -e ".[copilot,real,llm]"
$env:ANTHROPIC_API_KEY = "sk-ant-..."
$env:PKPD_PKSIM_CLI    = "C:\Program Files\Open Systems Pharmacology\PK-Sim 12.3\PKSim.CLI.exe"
python -m copilot.server
```

Then open <http://127.0.0.1:8765>.

## What it does

- **Left** — the compounds that have a benchmark snapshot, their status, and the
  run controls: LLM (Sonnet/Opus), max steps, and whether the agent **extracts
  the givens itself** from the handed report (on by default).
- **Centre** — first a *"What you're handed"* panel (the digitized
  concentration–time curves on a semi-log plot + the report with the answer
  redacted), then the live agent thread. Each tool call (`osp_read_report`,
  `osp_record_givens`, `osp_sweep_methods`, `osp_optimize`, …) becomes a step
  card; the optimizer's per-eval progress streams into it; each result shows its
  GMFE.
- **Right** — the best GMFE so far, the model structure the agent adopts, and a
  **Context the agent sees** selector: the finished projects it may read as
  leave-one-out analogues (default all, target always withheld; built-status
  dots; all/none). Unchecking narrows the reference library; unchecking all runs
  de-novo.

## Two honest defaults

- **Self-extract (context data lake).** By default the pre-digested givens are
  withheld and the agent builds them itself: it reads the handed report and raw
  data (`osp_read_report` / `osp_list_studies` / `osp_read_study`) and records a
  provenance-tagged parameter table (`osp_record_givens`), which then drives the
  fix-vs-fit split. Extraction is confined to the handed materials — **no web
  access** — so the agent can never pull the target's own published model off the
  internet.
- **Leave-one-out.** The context selector controls which *other* projects the
  agent reads; the selected compound's own reference (methods, fitted values,
  held-out data) is never readable, so held-out grading stays honest.

The raw report+data tree the agent reads lives in `../pbpk-realworld/` (generate
it with `python -m examples.build_realworld_projects --write`).

## Interactivity

- **Collapsed detail** — the noisy sweep grid and per-eval optimizer log are
  hidden inside each step; the header shows the live best GMFE and the outcome
  line. Click a step to expand its detail.
- **Toggle curves** — click a legend entry in the "What you're handed" plot to
  show/hide that study's curve.
- **Steer the run** — a composer note ("focus on oral absorption", "try
  CYP3A4 + P-gp") is passed to the agent as a modeler note on the goal.
- **Live givens table** — as the agent self-extracts, its provenance-tagged
  parameter table (given / judged / fitted) fills the right rail.
- **Deliverable** — each finished run writes a rerunnable handoff folder
  (`<compound>/deliverable/`: `adopted_model.json`, `parameter_table.csv`,
  `observed_curves.csv`, `README.md`); **Download deliverable** zips it.
- **New task from here** — follow-up chips reconfigure a control (more steps,
  de-novo, all context, Opus) and re-run in one click.

The held-out grade and the report (`.html` / `.json` / `.pdf`) are written to the
compound's `report/` folder, same as the CLI.

## Notes

- One run at a time (the server captures optimizer stdout, which is process-global).
- Set `PBPK_COPILOT_PORT` / `PBPK_COPILOT_HOST` to change where it binds.
- This is a thin UI over the same engine the CLI and scoreboard use — no separate
  modeling logic.
