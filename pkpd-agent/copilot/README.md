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

## Less scaffolding, more real modeling

Three changes make the task closer to what a modeler actually faces, and the
evaluation closer to the Kuepfer 2016 best-practice workflow:

- **Structure-blind (hard) mode** — a *Structure-blind* toggle runs against the
  `*-Model.hard_blanked.json` snapshot (the process skeleton and the metabolizing
  enzyme identity are stripped) plus the hard input (a
  `candidate_clearance_molecules` pool instead of the given enzyme). The agent
  then has to **discover the mechanism** — which enzyme clears the drug, which
  processes to add — not just fit values into a handed structure. Combine it with
  self-extract for the least-scaffolded task (structure *and* physchem withheld).
- **Staged IV→PO guidance** — `osp_inspect` now returns `method_guidance`: when
  both IV and PO data are present it tells the agent to fit in stages (IV sets
  distribution + clearance, then PO holds those and fits only absorption), and it
  flags when saturable kinetics are unidentifiable (fewer than two dose levels).
  Guidance from the data shape, never the answer.
- **PK-parameter goodness-of-fit** — beyond GMFE, every fit reports Cmax / tmax /
  AUC / t½ fold-errors (`pk_parameter_gmfe`), so a good pointwise GMFE with a bad
  peak or terminal slope surfaces as the structural clue it is. Shown on the
  editor result line and returned to the agent after each `osp_optimize`.

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
- **Open a support-context report** — click a project name in the "Context the
  agent sees" list to pop open that reference model's **full** report (its
  methods, fitted values, GMFE) and curves. This is the already-built modeling
  the agent leans on by analogy; the target's own report is redacted in the
  main panel and only openable here with an explicit "you're the modeler" note.
- **Model topology** — the viewer draws the whole-body PBPK structure: the
  fixed organ network + arterial/venous/portal blood flow in grey, and the
  drug-specific model the modeler chose/fit in teal (distribution & permeability
  methods, the metabolizing enzyme on the liver, renal GFR on the kidney, oral
  absorption on the gut). Fixed vs decided at a glance.
- **Overlay the reference fit** — inside that viewer, **▶ Run the model** runs
  the finished reference snapshot (forward simulation only, no fitting — fast)
  and overlays its simulated curves (lines) on the observed data (open points),
  per-study toggle + hover. Result cached on disk. Needs PKSim.CLI; the report's
  static goodness-of-fit figures render regardless.
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

## Reaching it from another computer

- **Same network:** start with `PBPK_COPILOT_HOST=0.0.0.0`; the console prints the
  LAN URLs. Open one from the other machine (allow the port in the firewall).
- **Any network:** run a tunnel to the port —
  `cloudflared tunnel --url http://localhost:8765` (public https URL, no account)
  or Tailscale (private mesh). The compute, PK-Sim and your API key stay on the host.
- **Protect it first.** There is no login by default. Before exposing it beyond a
  trusted LAN, set `PBPK_COPILOT_TOKEN=<something>`; then only URLs opened as
  `…/?token=<something>` get in (the server sets a cookie, so fetch and the SSE
  stream carry it afterwards). The startup banner prints the ready-to-use URL.

## Running in parallel (faster)

PK-Sim calls each use their own temp dir, so fits are independent and can run
concurrently. Two knobs (both default to serial, so behavior is unchanged unless
you opt in), sized near the machine's core count:

- `PKPD_JOBS=N` — run a build's method **sweep** with N fits at once (~N× faster
  sweep, the slowest phase of a build).
- `PKPD_COPILOT_JOBS=N` — allow N **agent runs** at once; the "Autonomous run all"
  scoreboard also has a *parallel* selector. Each run's stdout is routed
  per-thread so their live streams don't interleave.

Keep `PKPD_JOBS × PKPD_COPILOT_JOBS` around the core count — each PK-Sim process
uses a core and a few hundred MB, so over-subscribing thrashes.

## Stopping a run

**Stop** (single build), the scoreboard's **Stop**, and each row's **✕** now
cancel the run *on the server*, not just in the browser: the button POSTs to
`/api/cancel/{run_id}` and the agent loop checks that flag at every step/tool
boundary and finishes cleanly, freeing the worker slot. It can't interrupt a
PK-Sim call already in flight (that subprocess runs to its own timeout), but no
new step or fit is started, so a stop takes effect within one tool call. Closing
the page only stops *watching*; use a Stop button (or `POST /api/cancel/<id>`) to
stop the work itself.

## Clearing the cache

The only thing cached on disk is the reference-fit overlay (`▶ Run the model`),
one `.reffit.json` per compound under `OSP-PBPK-Model-Library/*/benchmark/`. The
agent build is never cached, and `/api/try` + Explore run live in throwaway temp
dirs. To clear it:

```powershell
python -m copilot.server --clear-cache     # deletes the .reffit.json files, reports the count
```

The in-memory copy of that cache clears whenever you restart the server. (Raw
equivalent: delete `OSP-PBPK-Model-Library\*\benchmark\.reffit.json`.)

## Notes

- Set `PKPD_COPILOT_JOBS`/`PKPD_JOBS` to 1 (default) for one run at a time.
- Set `PBPK_COPILOT_PORT` / `PBPK_COPILOT_HOST` to change where it binds.
- This is a thin UI over the same engine the CLI and scoreboard use — no separate
  modeling logic.
