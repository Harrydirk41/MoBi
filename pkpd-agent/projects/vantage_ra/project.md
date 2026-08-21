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

## Clinical endpoints
- **DAS28-CRP** (disease severity) and **ACR20/50/70** (response). In the model, DAS28-CRP is a
  Hill-sum of the 9 cell densities (Treg negative); ACR is derived from the DAS28 change.
- Trials: MTX (Strand 1999), ADA (OPTIMA), TCZ (ROSE), TCZ-refractory (Emery/RADIATE).

## Stage-1 benchmarks (all model-agnostic code, keyed off `network.json` + `spec.json`)
scope · topology · signs · readout-mapping · parameters · sensitivity. See
`RA_QSP_AGENT_TASK.md` for the measured ladder and the biology-determined-vs-model-committed
finding.

## Downstream (2–7) task objectives (config in `tasks.json`)
1. **Trial design** — predict the second-line response of MTX-inadequate responders escalated
   to a biologic; choose drug + sequencing.
2. **Calibration** — fit a PD parameter (e.g. `KD_TCZ`) so the model reproduces a real trial.
3. **Vpop generation** — sample disease drivers, reweight to a target DAS28 distribution.
4. **Drug design** — design a new anti-cytokine biologic (pathway + efficacy).
5. **Validation** — reproduce the held-out TCZ-in-dual-IR validation vs RADIATE.

## The one input that cannot be derived
`spec.json:gsa_top` — the global-sensitivity ranking comes from a figure (Fig 9), read by a
vision step or by hand. Everything else in `spec.json` is auto-inferable from `network.json`
(`--infer`) or read by the LLM structure extractor.
