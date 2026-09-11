# Continuous Claude Code agent in a sealed sandbox

This is the substrate for letting a general agent (Claude Code) **build a PBPK model
autonomously** - deciding what to do, in what order, until it judges the model done - instead
of running through the scripted `run_llm_build` decision loop. It removes the scaffolding
that was doing part of the modeling judgment for the agent, and it does so **without letting
the agent see the answers**.

## The one non-negotiable: two isolated worlds

A filesystem-capable agent can `cat` the reference model, the evaluation report, the
`answer_key/`, and the held-out concentrations. So the design is two trees, and only one is
ever given to the agent:

```
sandbox/<Model>/
  workspace/          <- answer-free. The agent runs HERE. Give the CC container ONLY this.
    model.blanked.json    whole-body snapshot, fitted values + observed data stripped
    task.input.json       objective, biology, physchem priors, BUILDING clinical data
    TASK.md               what to do
    .osp_session.json     the agent's own best-so-far/sweep state (created as it works)
  judge/              <- the answers. NEVER mounted into the agent's container.
    reference.json        the OSP reference model (fitted)
    heldout.observed.json the held-out concentration curves (what it's graded on predicting)
    split.json            building vs verification split
    forbidden.json        tokens the verifier scans for (reference fitted numbers, held-out concs)
    manifest.json         what was sealed/redacted/stripped
```

The wall is only as real as the deployment: run the CC session in a container/VM whose
filesystem contains **the `pkpd-agent` package and `workspace/` only** - not
`OSP-PBPK-Model-Library`, not `judge/`. On the same open filesystem the seal is advisory.

## What the seal removes (and why the model still runs)

The library's blanked snapshots were blanked at the Compound/Formulation level only. The
sandbox additionally makes the sealed copy answer-free:

- **Fitted values baked into `Simulations` overrides** (e.g. the fitted gut-wall permeability
  in every gut Neighborhood) and into fitted `ExpressionProfiles` reference concentrations
  are **redacted** - simulation overrides dropped (the building-block value governs),
  everything else neutralised and re-tagged.
- **`ObservedData`, `ParameterIdentifications`, and per-simulation `OutputMappings`** are
  **stripped**. The harness scores against the building data passed to the tools, not these
  blocks (and `osp_cli` already blanks `ParameterIdentifications` when it prunes), so removing
  them does not change how the model runs - it only removes the observed concentrations
  (building + held-out) and the saved fit targets from the file.

> A held-out study's **name/design** (author, route, dose, cohort) is left readable - a
> modeler predicting a named verification study legitimately knows it. Only its
> **concentration curve** is the answer, and that is gone.

## Run it

One CLI, `pkpd-bench`, does the whole flow. Install it once so it's on PATH and runs from
any directory (no `-m examples`, no cwd coupling):

```bash
pip install -e .        # gives you the `pkpd-bench` command
```

```bash
# 1. seal a case (or --all)
pkpd-bench seal --model Triazolam --out ../sandbox

# 2. verify the seal - exits non-zero on any leak (use it as a CI gate)
pkpd-bench verify --sandbox ../sandbox/Triazolam
pkpd-bench verify --library          # structural sweep over every library blanked snapshot

# 3. point the agent at the workspace. In its (isolated) environment it drives the tools:
cd ../sandbox/Triazolam/workspace
pkpd-bench inspect  --workspace .
pkpd-bench options  --workspace .
pkpd-bench sweep    --workspace . --estimate '{"Lipophilicity":[0,5]}'
pkpd-bench optimize --workspace . --estimate '{"Intrinsic clearance":[0.01,100]}' --structure '{"add_processes":[...]}'
pkpd-bench best     --workspace .    # best model so far

# 4. AFTER the agent stops, grade held-out prediction, and re-check for leaks:
pkpd-bench judge  --sandbox ../sandbox/Triazolam
pkpd-bench verify --sandbox ../sandbox/Triazolam --transcript session.log
```

Every subcommand prints one JSON object on stdout (progress goes to stderr), and state
persists in `workspace/.osp_session.json`, so the agent can stop and resume. If the console
script isn't on PATH, `python -m pkpd_agent.bench.cli <cmd> ...` is identical.
`PKPD_PKSIM_CLI` (or `--pksim`) points at `PKSim.CLI.exe`.

## What the agent is - and isn't - being tested on

Given inputs are trusted; the optimizer does the fitting; the agent supplies the **judgment**:
fit-vs-fix clearance (the IVIVE decision), which parameters to free without collinearity, and
which structural mechanism a systematic misfit calls for. `sweep` is a deterministic method
search (not withheld knowledge) - the agent runs it, the objective picks the winner.

## Known follow-ups (need a live PK-Sim run)

- The library's blanked snapshots themselves still carry the `Simulations`/`ExpressionProfile`
  fitted values (the sandbox redacts them per-seal). Extending the generator's blanking to
  those blocks is the permanent single-source fix; `verify_no_cheat --library` currently
  reports the 13 models that leak. Fixing them needs a PK-Sim re-build to confirm the models
  still run.
- The agent's tools expose Compound/Formulation/process parameters, not `ExpressionProfile`
  reference concentrations. Models whose fitted answer lives only there (Alprazolam/Midazolam
  ethnicity-scaled CYP3A4) need tool support or explicit exclusion.
- Wrapping `osp_agent_cli` as an MCP server (native CC tool integration) is a thin adapter
  over the same subcommands, if stdio JSON tools are preferred over Bash calls.
