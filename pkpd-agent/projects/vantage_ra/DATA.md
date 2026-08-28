# Vantage RA — project data layout

This project separates the four kinds of input a QSP study consumes, so each fact
lives where it belongs. The loader (`qsp_config.load_tasks`) merges them back into one
config at load time; every split file is optional, so a project that keeps everything
in `tasks.json` still loads unchanged.

```
vantage_ra/
├── model/
│   ├── spec.json         # model-structure description (readout targets/aliases/GSA)
│   └── parameters.json   # parameter catalog: nominal values + plausible PHYSIOLOGICAL ranges
├── data/
│   ├── calibration.json  # TRAINING trials the model is fit to + baseline-severity target
│   └── validation.json   # HELD-OUT trial, used only to qualify predictions
├── scenarios.json        # protocols to SIMULATE (treat-to-target sequence, target probes)
├── tasks.json            # TASK SETUP: readout mapping, timeline, endpoints, objectives
└── DATA.md               # this file
```

## What goes where — and why

| File | Holds | Changes when… |
|---|---|---|
| `model/spec.json` | model-structure description: readout targets, node aliases, top GSA parameters | the **model** changes |
| `model/parameters.json` | each parameter's nominal + plausible physiological range | the **model** changes (not the clinical data) |
| `data/calibration.json` | aggregate outcomes of the trials the model is **fit** to; the clinical baseline-severity distribution | you use **different training trials** |
| `data/validation.json` | the **held-out** trial's aggregate outcomes | you use a **different qualification trial** |
| `scenarios.json` | protocols/regimens to run for **prediction** and target exploration | you ask a **different what-if** |
| `tasks.json` | which CSV column plays which role, the timeline, the endpoints, and the objectives/targets each analysis matches | the **analysis setup** changes |

**The key distinction the layout encodes:**

- **Parameter ranges are a MODEL property, not data.** The plausible range of a rate
  constant is intrinsic to the mechanism (sourced from the paper's supplementary
  material); swapping the clinical data does not change it — so it lives under `model/`.
  *Which* subset is fit or varied, and how far the search bounds are narrowed for a
  given run, is a task decision and stays in `tasks.json`.
- **Calibration vs. held-out data are kept physically separate**, mirroring the QSP
  qualification principle: never validate on data you fit to.
- **`scenarios.json` is not data** — it is a specification of experiments to simulate.
  It carries no observed outcomes; the outcomes come out of the model.

## Data character

All entries under `data/` are **aggregate published-trial summaries** — response rates
(ACR20/50/70, remission) and the baseline-severity distribution — **not individual
patient records** (those are never published). This is exactly why the virtual
population must be *reverse-engineered* to match the aggregate: the per-patient data
does not exist to fit against.

## Provenance

- Parameter ranges: paper ESM (supplementary material).
- Calibration trials: Strand 1999 (MTX), OPTIMA/Kavanaugh (ADA), ROSE/Yazici 2012 (TCZ).
- Held-out trial: RADIATE (Emery 2008, TCZ in anti-TNF inadequate responders).
- Model: Bedathuru et al., *npj Systems Biology and Applications*, 2024.
