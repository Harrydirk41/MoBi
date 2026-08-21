# Project: Vantage RA QSP Model

This is a **project descriptor** — everything specific to *this* model/paper lives here as
data + prose. The code under `pkpd_agent/` is model-agnostic; point it at this folder to run
the Stage-1 evaluation suite on this model. Nothing about rheumatoid arthritis is hardcoded
in the engines.

## The model
- **Paper:** Bedathuru et al., *A multiscale, mechanistic model of Rheumatoid Arthritis to
  enable decision making in late-stage drug development*, npj Syst Biol Appl (2024) 10:126.
- **Model file:** `RA-QSP-Model/Vantage RA QSP Model v1.0.sbproj` (SimBiology). 59 species,
  ~524 parameters, ~101 reactions. Virtual population of 300 patients.
- **Structure dump:** `network.json` (via `sb_network_json.m` / `dump_network.py`, or from an
  SBML export via `sbml_to_network_json.py`). This is the sole structural input.
- **Biological cast (derived, not listed here):** 9 cells (FLS, Endothelial, Macrophages,
  BCells, PlasmaCells, CTL, Th1, Th17, Treg) + 17 mediators (TNFa, IL6, IL1b, IL17, IFNg,
  IL10, IL12, IL23, TGFb, GMCSF, VEGF, BAFF, MCP1, MIP3, RANTES, CAM, AutoAb).

## Data files
| file | what it is | used by |
|---|---|---|
| `network.json` | model structure (species, rules, params) | all Stage-1 benchmarks |
| `RA-QSP-Model/Vpop1.xlsx` | 300-patient parameter sets | downstream (trial/Vpop/validation) |
| `41540_2024_454_MOESM2_ESM.xlsx` | 130 documented params + units + refs | parameter benchmark |
| `41540_2024_454_MOESM1_ESM.xlsx` | clinical trial reference data | trial/validation scoring |
| `original_paper.pdf` (Fig 9) | global sensitivity top-20 | sensitivity benchmark (`spec.json:gsa_top`) |
| `spec.json` | Stage-1 structure spec (readout/drug patterns, aliases, gsa_top) | Stage-1 benchmarks |
| `tasks.json` | downstream (2–7) task data (drugs, params, targets, columns) | 2–7 loop tools + runners |

## Clinical endpoints
- **DAS28-CRP** (disease severity) and **ACR20/50/70** (response). In the model, DAS28-CRP is a
  Hill-sum of the 9 cell densities (Treg negative); ACR is derived from the DAS28 change.
- Trials: MTX (Strand 1999), ADA (OPTIMA), TCZ (ROSE), TCZ-refractory (Emery/RADIATE).

## Stage-1 benchmarks (all model-agnostic code, keyed off `network.json` + `spec.json`)
scope · topology · signs · readout-mapping · parameters · sensitivity. See
`RA_QSP_AGENT_TASK.md` for the measured ladder and the biology-determined-vs-model-committed
finding.

## Downstream (2–7) tasks — all RA specifics live in `tasks.json`
The engine is model-agnostic: `engines/qsp_tasks.py` (vocab-free scorers, spec builders,
numerical routines) + `engines/qsp_config.py` (`QSPTaskConfig`, loaded from `tasks.json`) +
`tools/qsp_{trial,fit,vpop,design,validate}_loop_tools.py` (build every prompt/answer-key off
the config). Runners: `run_llm_qsp_{trial,fit,vpop_gen,design,validate,full}.py --model ra`.
`tasks.json` holds the drug formulary, PD `fit_params`, `vpop_drivers`+target, `design_targets`,
the clinical trials, the `run_columns` semantic map, and the trial timeline — no drug, cytokine,
or trial name appears in the 2–7 code.
1. **Trial design** — predict the second-line response of MTX-inadequate responders escalated
   to a biologic; choose drug + sequencing.
2. **Calibration** — fit a PD parameter (e.g. `KD_TCZ`) so the model reproduces a real trial.
3. **Vpop generation** — sample disease drivers, reweight to a target DAS28 distribution.
4. **Drug design** — design a new anti-cytokine biologic (pathway + efficacy).
5. **Validation** — reproduce the held-out TCZ-in-dual-IR validation vs RADIATE.

## The per-model MATLAB boundary
The 2–7 tasks drive SimBiology through `examples/matlab/sb_run_vpop.m` /
`sb_sample_vpop.m`, which apply doses and emit the `run_columns` CSV. Those `.m` scripts are the
model-specific readout boundary (the run-time analogue of `network.json` for Stage-1): a new QSP
model ships its own readout script + `tasks.json`, and the Python above is unchanged.

## Porting to a new model — what the author actually prepares
The goal: the author worries only about *their* stuff (the model + their real clinical
numbers); everything else is derived.

| artifact | how it's produced | author effort |
|---|---|---|
| `network.json` | `dump_network.py` / `sbml_to_network_json.py` | none (auto) |
| `spec.json` | `--infer`, or `run_llm_extract` (LLM reads any naming, regressed vs RA) | review a draft |
| `tasks.json` role assignments (`readout_states`, `vpop_drivers`, `design_targets`, `fit_params`) | `run_llm_draft_tasks` (LLM drafts from `network.json`, regressed vs this RA `tasks.json`) | review a draft |
| `tasks.json` external fields (`clinical_trials`, `refractory_target`, `vpop_target`, dose `drugs`, `timeline`) | author supplies — real trial data + `.sbproj` dose names | **the author's own data** |
| `.m` readout script | reuse `sb_run_vpop.m` as-is if the model is same-shaped (state names come from `readout_states`); rewrite only for a different trial structure | none / rewrite |

So a same-shaped QSP model ports to: dump `network.json` → run the two extractors to
draft `spec.json` + `tasks.json` → paste the real clinical numbers + dose names → run.

### For a traditional modeler (no JSON, no schema)
Naming is never required to match this project — the extractors read whatever the model
calls things, and `run_columns` maps roles to the modeler's own column names. A GUI-native
modeler never touches JSON; they write a few plain sentences ("psoriasis model, severity
PASI, active band 6–20, drug secukinumab dose SEC_300mg, match UNCOVER-2 PASI75 wk12 = 77%,
calibrate KD_SEC ref 1e-10 M") and use one of two entry points:
- `run_llm_init` — linear: extract structure → draft roles → fill `tasks.json` from the
  description → `project_validate` reports leftovers in plain English → write the folder.
- `run_llm_onboard` — an AGENT drives it: inspects the model, builds the config, reads its
  own validation report, fixes the errors it can infer (a `vpop_driver` that near-misses a
  real parameter, etc.), and saves only when clean, telling the modeler which clinical
  numbers they still owe. `onboard_save` refuses while any validation ERROR remains.

`engines/project_validate.py` is the safety net either way: it checks every config name
against the real model and turns the silent foot-guns (a role that is not a real parameter,
an unfilled stub) into plain-English errors/warnings. The modeler edits nothing by hand.

## The one input that cannot be derived
`spec.json:gsa_top` — the global-sensitivity ranking comes from a figure (Fig 9), read by a
vision step or by hand. Everything else in `spec.json` is auto-inferable from `network.json`
(`--infer`) or read by the LLM structure extractor.
