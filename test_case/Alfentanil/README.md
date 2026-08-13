# Alfentanil benchmark case

Drop the source files for the Alfentanil masked-report benchmark here:

- `*.md`   — the evaluation report (e.g. `Alfentanil_evaluation_report.md`)
- `*.pdf`  — the report PDF, if you have it
- `*.json` — the OSP snapshot(s):
    - `AlfentanilModel.json`      — adult base model (fitted; 2 ParameterIdentifications)
    - `AlfentanilPediatrics.json` — pediatric extrapolation (no re-fitting)

## What this case is for

Two benchmark problems live in this folder:

1. **Adult build** (`AlfentanilModel.json`) — build + fit a PBPK model from the
   compound and 19 adult studies. Harder: it includes an actual estimation step.
2. **Pediatric extrapolation** (`AlfentanilPediatrics.json`) — predict the 3
   pediatric profiles from the adult model + child physiology, with **no
   parameter fitting**. Cleaner "recover the concentration profile" test.

The masking recipe: keep the description, data, and objective; mask the modeling
strategy, the fitted parameter values, and the results. Then check whether the
agent recovers a concentration profile / model close to the ground-truth snapshot.

Extract any snapshot with:

    python -m examples.extract_snapshot test_case/Alfentanil/AlfentanilModel.json
