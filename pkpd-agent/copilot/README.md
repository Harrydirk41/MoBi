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
  run controls: LLM (Sonnet/Opus), max steps, and the reference library toggle
  (off / same-type / all).
- **Centre** — the live agent thread. Each tool call (`osp_sweep_methods`,
  `osp_optimize`, …) becomes a step card; the optimizer's per-eval progress
  streams into it; each result shows its GMFE.
- **Right** — the best GMFE so far, the model structure the agent adopts, and the
  leave-one-out reference library (the target's own model is never included).

The held-out grade and the report (`.html` / `.json` / `.pdf`) are written to the
compound's `report/` folder, same as the CLI.

## Notes

- One run at a time (the server captures optimizer stdout, which is process-global).
- Set `PBPK_COPILOT_PORT` / `PBPK_COPILOT_HOST` to change where it binds.
- This is a thin UI over the same engine the CLI and scoreboard use — no separate
  modeling logic.
