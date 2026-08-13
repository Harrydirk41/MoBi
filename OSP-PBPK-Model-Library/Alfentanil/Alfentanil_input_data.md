# Alfentanil — Objective & Input Data

This note describes the **inputs** for the alfentanil PBPK modeling case: the
objective, the clinical data to be reproduced, and the physicochemical/in‑vitro
information you are given. It is written so an agent (or a person) can attempt to
**rebuild the model from scratch** without seeing the original modeling strategy
or the fitted results.

Everything here is derived from the two OSP snapshots in `json/`
(`Alfentanil-Model.json`, `Alfentanil-Pediatrics.json`) and the published
evaluation report (`Alfentanil_evaluation_report.md`). The machine-readable,
already-cleaned inputs live in `json_input/` and are produced by
`build_clean_input.py` (see [How the clean input was built](#how-the-clean-input-was-built)).

---

## 1. Objective

Alfentanil is a potent, fast/short-acting synthetic opioid used for anesthesia.
Two modeling tasks share the same compound:

- **Adult build** (`Alfentanil-Model`) — Build a whole-body PBPK model that
  reproduces the observed plasma concentration–time profiles after **IV and oral**
  dosing across 8 clinical studies (weight-based doses 0.015–0.075 mg/kg, plus
  absolute 1 mg IV and 4 mg PO). A handful of ADME parameters are unknown and must
  be **estimated** so that simulations match the observed data.

- **Pediatric extrapolation** (`Alfentanil-Pediatrics`) — Take the adult model and
  **predict** plasma profiles in children (~0.3–14 y) using PK‑Sim's built-in
  physiology and enzyme **ontogeny**. This is a pure prediction task: **no
  parameters are re-estimated.**

**Success criterion (either task):** simulated plasma concentrations fall close to
the observed data — the OSP report used geometric mean fold error (GMFE) with a
target of roughly ≤ 2, achieving GMFE 1.26 (IV), 1.45 (PO), 1.32 (all).

### Known biology (given background, not something to discover)

- Alfentanil is metabolized **solely by CYP3A4**; < 1% is excreted unchanged in urine.
- It is **not a P‑gp substrate**.
- Absorption is fully explained by **passive** permeation (no active uptake needed).
- Plasma protein binding partner: **α1‑acid glycoprotein**.

---

## 2. What is given vs. what must be estimated

PBPK separates *system* knowledge (physiology, taken from the database) from
*drug* knowledge (compound properties). For the drug, some properties are measured
in the lab / taken from literature (**given**), while a few that cannot be measured
reliably are **estimated** to fit the clinical data.

### 2.1 Given physicochemical / in‑vitro inputs

| Parameter | Value | Unit | Source |
|---|---|---|---|
| Molecular weight | 416.52 | g/mol | DrugBank DB00802 |
| pKa (base) | 6.5 | — | Jansson 2008 |
| Solubility (at reference pH 6.5) | 992 | mg/L | Baneyx 2014 |
| GFR fraction (renal) | 0.06 | — | Hanke 2018 |
| logD (reference) | ~2.1–2.2 | — | Baneyx 2014 / Jansson 2008 |
| fu (fraction unbound, reported range) | 8.6–12 | % | Gertz 2010 / Edginton 2008 / Almond 2016 |

> Note on fu: the literature reports a range (8.6–12%). In the reference model the
> final unbound fraction was refined during fitting, so in the clean input it is
> listed under *parameters to identify* (below) rather than as a fixed given.

### 2.2 Parameters to identify (values withheld — these are the unknowns)

These are the parameters the reference model estimated with PK‑Sim's Parameter
Identification. In the clean input their **values are hidden** so the task is
genuine; only their names/units are provided so you know what to solve for.

- Lipophilicity (Log Units)
- Fraction unbound (plasma, reference value)
- Specific intestinal permeability (transcellular) — cm/min
- Permeability (specific organ permeability) — cm/min
- CYP3A4 intrinsic clearance — l/min
- *(report also lists basolateral mucosa permeability, tuned to raise gut-wall metabolism)*

**Withheld as the answer key** (not in the clean input): the fitted values of the
above, the distribution/permeability **method choices**, and the simulated output.
See `Alfentanil_evaluation_report.md` §2.3.4 and §3.1 if you want to grade against them.

---

## 3. Clinical data (the target to reproduce)

All observed data are **plasma** concentrations in peripheral venous blood, time in
hours, concentration converted to **mg/L** (native units in the snapshot are µg/L or
mg/L). Arithmetic standard deviations are included where the source reported them.

### 3.1 Adult build data — 19 datasets, 8 studies, IV + PO

| Study | Route | Dose | # datasets |
|---|---|---|---|
| Ferrier 1985 | IV | 0.05 mg/kg | 2 (healthy + cirrhosis) |
| Kharasch 1997 | IV | 0.02 mg/kg | 1 |
| Kharasch 2004 | IV / PO | 0.015 / 0.06 mg/kg | 2 |
| Kharasch 2011 | IV / PO | 0.015 / 0.075 mg/kg | 2 |
| Kharasch 2011b | IV | 1 mg | 4 (±grapefruit, seq/simul) |
| Kharasch 2011b | PO | 4 mg | 4 (±grapefruit, seq/simul) |
| Kharasch 2012 | IV / PO | 0.02 / 0.043 mg/kg | 2 |
| Meistelman 1987 | IV | 0.02 mg/kg | 1 (adults) |
| Phimmasone 2001 | IV | 0.015 mg/kg | 1 |

Doses span weight-based (0.015–0.075 mg/kg) and absolute (1 mg IV, 4 mg PO). The
grapefruit-juice arms of Kharasch 2011b are control profiles from a DDI study;
grapefruit affects gut (not hepatic) CYP3A4, so IV arms are effectively unaffected.

### 3.2 Pediatric extrapolation data — 3 datasets, IV only

| Study | Route | Dose | Population |
|---|---|---|---|
| den Hollander 1992 | IV | 1.7 mg (individual I3) | children |
| Meistelman 1987 | IV | 20 µg/kg | mean child (n=8) |
| Goresky 1987 | IV | 50 µg/kg | typical individual |

The pediatric snapshot additionally defines **38 individual children** (age
**0.2–14 y**, weight **4.8–51 kg**) and **1 population** ("Meistelman 1987 Children"),
which supply the demographic/physiology inputs for the extrapolation.

---

## 4. Study designs & demographics

Beyond the raw curves, the clean input also carries:

- **`study_designs`** — dosing protocols: application type (IV bolus / oral),
  dose + unit, single vs. multiple, start time, formulation (solution).
- **`demographics`** — for each modeled individual: species, PK‑Sim population
  (e.g. `European_ICRP_2002`), gender, age, weight; and population blocks with the
  number of virtual individuals. Adults are typical European individuals (age 30);
  pediatrics are the real children listed in the source publications.

These let a simulator set up the correct dose and virtual subject without guessing.

---

## 5. The two snapshots at a glance

| | `Alfentanil-Model` (adult) | `Alfentanil-Pediatrics` |
|---|---|---|
| Role | base model — **built & fitted** | **extrapolation** to children |
| Individuals | 2 adults (age 30) | 38 children (0.2–14 y) |
| Populations | 0 | 1 |
| ParameterIdentifications | **2** (fitting happened here) | **0** (no re-fitting) |
| Observed datasets | 19 (IV + PO) | 3 (IV) |
| Params to identify (clean input) | 5 | 0 |

The classic PBPK workflow: **fit in adults where data is rich, then extrapolate to
children via physiology + ontogeny.** The pediatric file uses the published Hanke
2018 drug parameters directly; the adult file carries locally re-optimized values
(PK‑Sim "Parameter Identification 4", 2019) — so the drug parameters are not
byte-identical between the two, which is expected given the separate provenance.

---

## 6. Files

```
Alfentanil/
├── Alfentanil_evaluation_report.md   # full OSP report (contains the answer key)
├── Alfentanil_evaluation_report.pdf
├── Alfentanil_input_data.md          # this file
├── build_clean_input.py              # snapshot -> clean input JSON
├── json/                             # raw OSP snapshots (input + strategy + results)
│   ├── Alfentanil-Model.json
│   └── Alfentanil-Pediatrics.json
└── json_input/                       # CLEAN inputs (this is what you feed the agent)
    ├── Alfentanil-Model.input.json
    └── Alfentanil-Pediatrics.input.json
```

### Clean-input JSON schema (`osp-clean-input/v1`)

```
compound, objective, source_snapshot
data_overview            # counts, routes, study list
compound_identity        # name, small-molecule flag, protein partner, pKa
given_physicochemical[]  # parameter, value, unit, source   (measured/published)
parameters_to_identify[] # parameter, unit, note            (VALUES WITHHELD)
clinical_observed_data[] # dataset, study, route, dose, matrix,
                         #   time_h[], conc_mg_L[], conc_native[], sd_mg_L[]
study_designs[]          # application_type, dose, dosing_interval, start_time
demographics[]           # species, population, gender, age, weight (+ populations)
```

---

## How the clean input was built

`build_clean_input.py` reads an OSP snapshot and keeps only the input half,
routing anything that would leak the modeling strategy or results into the hidden
answer key. Concretely it:

1. tidies every observed profile → **time in h**, **concentration in mg/L** (from
   µg/L or molar via molecular weight), plus native values and arithmetic SD;
2. splits compound parameters by `ValueOrigin`: `Publication`/measured → **given**,
   `ParameterIdentification` → **parameters to identify** (name/unit only, value
   withheld);
3. extracts study designs (protocols) and demographics (individuals/populations);
4. **drops** fitted values, distribution/permeability method choices, the
   parameter-identification configuration, and simulation results.

Regenerate at any time (stdlib only, no OSP install needed):

```bash
python build_clean_input.py json/Alfentanil-Model.json
python build_clean_input.py json/Alfentanil-Pediatrics.json
# optional: --objective "..."   --outdir json_input
```
