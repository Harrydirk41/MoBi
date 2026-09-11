---
name: pkpd-modeling
description: >
  Build and grade pharmacology models as an autonomous agent in this repo — PBPK
  (whole-body PK), QSP (systems/immune networks), and popPK/PD (population). Use when asked
  to build, fit, benchmark, or grade a model here, or to drive the sealed-sandbox benchmark.
  The clean entry point is the `pkpd-bench` CLI; this skill is the input/output contract.
---

# Modeling-agent surface

This repo evaluates whether an agent can build pharmacology models to expert standard
**without seeing the answers**. All domain value is in tools + a held-out grader, so any
coding agent drives the same clean surface. Prefer these commands over the ~80 research
scripts in `examples/` (those are the lab notebook — see `examples/README.md`).

## PBPK — the primary, fully-clean CLI: `pkpd-bench`

Install once (`pip install -e .`) so it's on PATH and runs from any directory. Each
subcommand prints ONE JSON object on stdout (progress → stderr); state persists in
`workspace/.osp_session.json`.

**Setup & grading (human/orchestrator side, needs the answers):**
```
pkpd-bench seal   --model <Name> --out ../sandbox     # build answer-free workspace + judge tree
pkpd-bench verify --sandbox ../sandbox/<Name>         # assert no answer leaks (exit≠0 on leak)
pkpd-bench verify --library                           # structural leak sweep over all models
pkpd-bench judge  --sandbox ../sandbox/<Name>         # grade held-out prediction after the run
```

**Building the model (the agent, working INSIDE `<sandbox>/<Name>/workspace`):**
```
pkpd-bench inspect  --workspace .                     # objective, biology, priors, data, model
pkpd-bench options  --workspace .                     # editable params, legal methods, mechanisms
pkpd-bench sweep    --workspace . --estimate '{"Lipophilicity":[0,5]}'   # choose distribution method
pkpd-bench optimize --workspace . --estimate '{"<param>":[lo,hi]}' --structure '{"add_processes":[...]}'
pkpd-bench try      --workspace . --edits '{"parameters":{...}}'
pkpd-bench best     --workspace .                     # best model so far
```
(No PATH entry? `python -m pkpd_agent.bench.cli <cmd> ...` is identical.)

**The modeling job** (what the agent decides; the optimizer does the fitting):
1. `inspect` + `options` — read the given inputs. Each parameter has a `value_status`:
   `given` = measured, trust it; `placeholder` = unknown, determine it (the shown number is
   meaningless).
2. Pick the clearance/metabolism processes from the known biology.
3. If distribution/Vd is off, run `sweep` ONCE (it tries every partition × permeability method
   and re-fits your physchem under each, then adopts the best).
4. Choose which parameters to **fit vs fix**. Watch the IVIVE trap: an in-vitro Vmax being
   given does NOT mean clearance is set — the in-vivo kcat usually still needs fitting. Keep
   physically-needed parameters (e.g. oral absorption) even if they look redundant in-sample.
5. `optimize`, then read `recommendations` / `params_at_bound` / per-route bias; act on
   identifiability findings. Iterate until the fit is good AND parameters are identifiable, then stop.

You are graded afterwards on held-out studies (not in the workspace). Build physiology, don't
overfit the building data, and never look for the held-out data or reference values.

## The shared harness (same discipline for every domain)

- **Sealed sandbox** — the workspace has only the blanked model + building data + task; the
  reference, held-out curves, and fitted answers live in a separate `judge/` tree.
- **Held-out grading** — the agent's model and the reference are re-run on hidden studies;
  PASS = agent GMFE ≤ 1.25× the reference's.
- **No-cheat verifier** — scans the workspace (and optionally the agent transcript) for any
  leaked answer. `SANDBOX_CC.md` is the full runbook.

## QSP — systems/immune-network models

The whole pipeline runs from one command (needs MATLAB + SimBiology; the pure build stage runs
without them):
```
python -m examples.run_qsp_full_pipeline --model ra --live --out report.md   # + --sbproj/--modeldir/--paper-vpop
```
It builds the immune network from biology, calibrates it to baseline, checks stability, builds
a virtual population, and predicts a real trial (RA / DAS28 / ACR) — graded against the
validated published model AND the trial, with the second-line arm held out. Smaller stages:
`run_qsp_build_cell`, `run_qsp_build_network`, `run_qsp_calibrate`, `run_qsp_end_to_end`.

## popPK/PD — population models

The analysis engines are in place (NCA, population PK via pharmpy, NLME via nlmixr2, curve
fitting); the generic agent driver is `python -m examples.run_llm ...`. A sealed-sandbox
benchmark like PBPK's is the planned next step, not yet built.

## When unsure which command

PBPK model build/grade → `pkpd-bench`. QSP network/trial → `run_qsp_full_pipeline`. Anything
else (building benchmarks, extracting data, one-off analyses) → `examples/README.md`.
