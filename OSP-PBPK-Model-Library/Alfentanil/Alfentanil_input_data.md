# Alfentanil — Modeling Task (Agent-Facing)

This is the **problem statement** for the alfentanil PBPK modeling benchmark. It
describes the objective, the data you are given, what you must produce, and how
you will be graded. It deliberately contains **no modeling strategy and no
reference results** — recovering a sensible model is the task.

> **Leakage policy.** The reference model's structure, its fitted parameter
> values, the distribution/permeability method choices, and the achieved fit
> metrics are held in a separate `answer_key/` folder used only by the grader.
> Do **not** load `answer_key/`, the raw `json/` snapshots, or
> `Alfentanil_evaluation_report.md` when attempting the task — they contain the
> answers. Work only from `json_input/*.input.json` and this file.

The machine-readable task lives in `json_input/` (schema `osp-agent-task/v2`);
this note is the human-readable companion.

---

## 1. Objective

Alfentanil is a potent, fast- and short-acting synthetic opioid used for
anesthesia. There are two tasks (same drug, different question):

- **Adult build** (`Alfentanil-Model.input.json`) — Build a whole-body PBPK model
  that reproduces observed plasma concentration–time profiles after **IV and
  oral** dosing across 8 clinical studies (weight-based 0.015–0.075 mg/kg, plus
  absolute 1 mg IV / 4 mg PO). Choose the structure and any parameters you cannot
  read off the given data, and justify them.

- **Pediatric extrapolation** (`Alfentanil-Pediatrics.input.json`) — **Predict**
  plasma profiles in children (~0.3–14 y) by scaling a PBPK model with
  age-appropriate physiology and enzyme **ontogeny**. This is a prediction task:
  **do not estimate parameters** from the pediatric data — reuse adult drug
  properties and let physiology do the work.

### Public background (given, from the literature — not model output)

- Metabolized primarily by hepatic **CYP3A4**; CYP3A4 abundance matures with age.
- Renal excretion of unchanged drug is **< 1%** (negligible).
- **Not** a P-glycoprotein substrate.
- Clinically dosed IV; some drug-interaction studies also report oral profiles.

---

## 2. What you are given

Contained in each `*.input.json` under `given_data`:

### 2.1 Literature physicochemical / in-vitro data

| Parameter | Value | Unit | Source |
|---|---|---|---|
| Molecular weight | 416.52 | g/mol | DrugBank DB00802 |
| pKa (base) | 6.5 | — | Jansson 2008 |
| Solubility (at pH 6.5) | 992 | mg/L | Baneyx 2014 |
| logD (pH 7.4) | 2.1–2.2 (reported) | — | Baneyx 2014 / Jansson 2008 |
| Fraction unbound in plasma | 8.6–12 (reported range) | % | Gertz 2010 / Edginton 2008 / Almond 2016 |
| Plasma protein partner | α1-acid glycoprotein | — | — |

These are **inputs**, not answers. Some are single measured values; some are
reported *ranges* — where a range is given, you decide the value to use (fix it,
or estimate within the range) and justify it.

### 2.2 Clinical observed data (the target to reproduce)

All observed data are **plasma** concentrations, time in **hours**, concentration
in **mg/L** (converted from the source µg/L or mg/L), with arithmetic SD where the
source reported it. Each record carries its study, route, dose, and matrix.

**Adult build — 19 datasets, 8 studies, IV + PO:**

| Study | Route | Dose | # |
|---|---|---|---|
| Ferrier 1985 | IV | 0.05 mg/kg | 2 |
| Kharasch 1997 | IV | 0.02 mg/kg | 1 |
| Kharasch 2004 | IV / PO | 0.015 / 0.06 mg/kg | 2 |
| Kharasch 2011 | IV / PO | 0.015 / 0.075 mg/kg | 2 |
| Kharasch 2011b | IV | 1 mg | 4 |
| Kharasch 2011b | PO | 4 mg | 4 |
| Kharasch 2012 | IV / PO | 0.02 / 0.043 mg/kg | 2 |
| Meistelman 1987 | IV | 0.02 mg/kg | 1 |
| Phimmasone 2001 | IV | 0.015 mg/kg | 1 |

**Pediatric extrapolation — 3 datasets, IV only:** den Hollander 1992 (1.7 mg),
Meistelman 1987 (20 µg/kg, mean child), Goresky 1987 (50 µg/kg). Plus **38
individual children** (age 0.2–14 y, weight 4.8–51 kg) and **1 population** as the
demographic/physiology input for the extrapolation.

### 2.3 Study designs & demographics

- **`study_designs`** — dosing protocols: application type (IV bolus / oral),
  dose + unit, single vs. multiple, start time.
- **`demographics`** — per modeled subject: species, PK-Sim population, gender,
  age, weight; plus population blocks. (Adults: typical European individuals;
  pediatrics: the real children from the source publications.)

---

## 3. The unknowns (you decide)

Not everything is measurable from the given data. You must choose:

- a **distribution / partition-coefficient** approach and a **cellular
  permeability** approach;
- an **absorption** model (for the oral arms);
- values for compound parameters that can't be derived from the literature data,
  which you may need to **estimate** so the model reproduces the observed curves;
- how to represent **elimination** consistent with the known biology (CYP3A4
  metabolism; negligible renal).

The reference model's specific choices are withheld. Two reasonable modelers may
land on different-but-valid models — that's expected; the rubric (Section 5)
rewards fit **and** physical plausibility, not matching a specific answer.

---

## 4. What you must produce

Return an object with the shape in `what_you_must_produce.submission` (or add it
back into the input JSON under a top-level `submission` key):

```
submission
├── structural_model      # distribution model, absorption model,
│                         #   elimination pathway(s), engine, notes
├── parameters[]          # every value you fixed or estimated beyond the
│                         #   given data: {parameter, value, unit,
│                         #   fixed_or_estimated, rationale}
├── predicted_profiles[]  # per dataset: {dataset, time_h[], pred_conc_mg_L[]}
│                         #   ALIGNED to the observed dataset names & time grid
└── self_assessment       # your own GMFE / %within-2-fold + plausibility notes
```

`predicted_profiles` must use the **same dataset names and the same `time_h`
grid** as `clinical_observed_data`, so predictions can be paired point-for-point
with observations.

---

## 5. How you will be graded

Three independent dimensions (also in `evaluation_rubric`). You never see the
reference answers — aim for a model that fits **and** makes physical sense.

1. **Data fit** — predicted vs. observed plasma concentrations.
   - **GMFE** = `exp(mean(|ln(pred/obs)|))` over paired points — target ≤ 2 (good ≤ 1.5).
   - **% within 2-fold** of observed — the higher the better.
   - Reported overall, by route (IV vs PO), and by study.

2. **Parameter plausibility** — every parameter you set is physically sensible:
   - fraction unbound in (0, 1]; lipophilicity (logP/logD) ≈ [−2, 7];
   - clearance ≤ relevant organ blood flow (hepatic ≈ 1.5 L/min in an adult);
   - permeabilities > 0 and physiological; volume of distribution ≈ 0.05–50 L/kg;
   - all rate constants and doses positive; units consistent.

3. **Output plausibility** — simulated profiles behave like real PK:
   - concentrations ≥ 0 everywhere; no numerical blow-up;
   - IV bolus declines monotonically from t₀; oral rises to Cmax then declines;
   - terminal half-life consistent (within ~2×) across doses unless nonlinearity
     is justified; dose-normalized AUC roughly constant across doses, or the
     nonlinearity is explained; Cmax and AUC increase with dose.

A self-assessment against these three is required in your submission.

---

## 6. Grading (for the grader — not part of the agent's input)

Submissions are scored by `grade_submission.py` in two layers that mirror the
rubric:

- **Numerical (deterministic, stdlib):** pairs `predicted_profiles` against the
  observed data (interpolating onto the observed time grid) and computes **GMFE**
  and **% within 2-fold** overall / by route / by study; runs the rule-based
  **parameter** bound checks and **output** plausibility checks. Optionally reports
  auxiliary *closeness to the reference* from `answer_key/` (not a primary score —
  many valid models differ from the reference).
- **Agentic (Claude judge, optional):** reads the numerical scorecard plus the
  agent's structural choices and rationales and renders the **physical-reasoning**
  verdict the numbers can't — mechanistic soundness, whether each flag is a real
  problem, a 0–5 score per dimension, an overall verdict, and actionable feedback.
  Runs when `anthropic` is installed and `ANTHROPIC_API_KEY` is set (or `--reason`);
  otherwise the numerical scorecard is produced alone.

```bash
# 1. (demo) fabricate a submission just to exercise the pipeline
python make_demo_submission.py json_input/Alfentanil-Model.input.json --out demo_submission.json

# 2. grade it (add --reason for the Claude physical-reasoning layer)
python grade_submission.py \
    --input json_input/Alfentanil-Model.input.json \
    --submission demo_submission.json \
    --key answer_key/Alfentanil-Model.answer_key.json
# writes scorecards/<submission>.scorecard.json

python grade_submission.py --selftest   # checks the GMFE math
```

## 7. Files

```
Alfentanil/
├── Alfentanil_input_data.md            # THIS problem statement (agent-facing)
├── build_clean_input.py                # snapshot -> input + answer key
├── grade_submission.py                 # numerical + agentic grader
├── make_demo_submission.py             # synthetic submission for pipeline testing
├── problem_cards.json                  # curated PUBLIC context (objective, literature)
├── json_input/                         # AGENT INPUT — start here, no leaks
│   ├── Alfentanil-Model.input.json
│   └── Alfentanil-Pediatrics.input.json
├── answer_key/                         # GRADING ONLY — do not show the agent
│   ├── Alfentanil-Model.answer_key.json
│   └── Alfentanil-Pediatrics.answer_key.json
├── scorecards/                         # grader output (generated)
├── json/                               # raw snapshots (full truth) — grader/source only
└── Alfentanil_evaluation_report.md     # the original OSP report — contains answers
```

### Regenerating the input/key

```bash
python build_clean_input.py json/Alfentanil-Model.json
python build_clean_input.py json/Alfentanil-Pediatrics.json
# writes json_input/<stem>.input.json (agent) and answer_key/<stem>.answer_key.json (grader)
```

`problem_cards.json` holds only public, literature-sourced context. Everything the
builder finds that is model-derived (fitted values, method choices, processes) is
routed to `answer_key/`, never into `json_input/`.
