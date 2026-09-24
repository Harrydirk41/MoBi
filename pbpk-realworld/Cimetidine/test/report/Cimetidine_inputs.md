# Building and evaluation of a PBPK model for cimetidine in healthy adults


> You are handed this report and the clinical data files. Extract the

> physicochemical and in-vitro inputs yourself, decide the model

> structure, and fit the parameters. The reference model's chosen

> methods and fitted values are withheld.


# Background

Cimetidine is a histamine H2 receptor antagonist that inhibits stomach acid production. It is mainly used as an antacid for the treatment of gastric and duodenal ulcers, Zollinger-Ellison syndrome and esophageal reflux.

The herein presented model was developed and published by Hanke et al. ([Hanke 2020](#5-references)) and adjusted later on to PK-Sim V10 by refitting CYP3A4 K<sub>i</sub> and MATE1 K<sub>i</sub>.

Cimetidine is mainly excreted unchanged via the kidneys (40–80% of the dose) with a high renal clearance of 400 ml/min. Metabolism is reported to account for 25– 40% of of the total elimination of cimetidine, with less than 2% of the dose excreted unchanged with the bile. Cimetidine inhibits several transporters and CYP enzymes and it is recommended by the FDA as strong inhibitor of OCT2/MATE and as weak inhibitor of CYP3A4 and CYP2D6 for the use in clinical DDI studies and drug labeling.

The cimetidine model was established using 27 clinical studies, covering a dosing range from 100 to 800 mg. The final model applies active uptake of cimetidine into the liver by OCT1,
uptake into the kidney by OAT3 and secretion from the kidney into the urine by MATE1, as well
as an unspecific hepatic clearance and passive renal glomerular filtration.

The herein presented model building and evaluation report evaluates the performance of the PBPK model for cimetidine in (healthy) adults.

## 2.2 Data<a id="data"></a>

### 2.2.1 In vitro / physico-chemical Data

A literature search was performed to collect available information on physiochemical properties of cimetidine. The obtained information from literature is summarized in the table below. 

| **Parameter**   | **Unit** | **Value**       | Source                                                       | **Description**                                 |
| :-------------- | -------- | --------------- | ------------------------------------------------------------ | ----------------------------------------------- |
| MW              | g/mol    | 252.34    | [Wishart 2006](#5-references) | Molecular weight                                |
| pK<sub>a</sub>1 | 6.93  | (base)          | [Avdeef 2001](#5-references)                       | Acid dissociation constant                      |
| pK<sub>a</sub>2 | 13.38  | (acid)          | [Wishart 2006](#5-references)                    | Acid dissociation constant                      |
| Solubility (pH) | mg/L     | 24.00 (6.8) | [Avdeef 2001](#5-references) | Water solubility                               |
| f<sub>u</sub> |   %      | 78.00 | [Taylor 1978](#5-references) | Fraction unbound in plasma                      |
| B/P ratio              |         | 0.98 | [Somogyi 1983](#5-references) | Blood to plasma ratio        |
| OCT1 K<sub>m</sub> | μmol/l | 2600 | [Umehara 2007](#5-references) | Michaelis-Menten constant |
| OAT3 K<sub>m</sub> | µmol/l | 149     | [Tahara 2005](#5-references) | Michaelis-Menten constant |
| MATE1 K<sub>m</sub> | µmol/l | 8.0 | [Ohta 2005](#5-references) | Michaelis-Menten constant |
| OCT1 K<sub>i</sub> | µmol/l | 104     | [Ito 2012](#5-references) | Inhibition constant for competitive inhibition |
| OCT2 K<sub>i</sub> | µmol/l | 124 | [Ito 2012](#5-references) | Inhibition constant for competitive inhibition |
| MATE1 K<sub>i</sub> (refitted in PK-Sim V10) | µmol/l | 3.8 (0.65) | [Ito 2012](#5-references) | Inhibition constant for competitive inhibition |
| CYP3A4 K<sub>i</sub> (refitted in PK-Sim V10) | µmol/l | 268 (30.51266) | [Wrighton 1994](#5-references) | Inhibition constant for competitive inhibition |

### 2.2.2 Clinical Data

A literature search was performed to collect available clinical data on efavirenz in healthy adults.

#### 2.2.2.1 Model Building

The following studies were used for model building:

| Publication                       | Arm / Treatment / Information used for model building        |
| :-------------------------------- | :----------------------------------------------------------- |
| [Bodemar 1981](#5-references)     | Peptic ulcer patients receiving a single intravenous dose of 200 mg and oral doses of 200, 400 and 800 mg |
| [Morgan 1983](#5-references)      | Peptic ulcer patients receiving a single intravenous dose of 200 mg (5 min infusion) |
| [Bodemar 1979](#5-references)     | Healthy subjects receiving single oral doses of 200 and 400mg (tablet) |
| [Walkenstein 1978](#5-references) | Healthy subjects receiving a single oral dose of 300mg (solution) |
| [D'Angio 1986](#5-references)     | Healthy subjects receiving a single oral dose of 300mg (tablet) |

#### 2.2.2.2 Model verification

The following studies were used for model verification:

| Publication                       | Arm / Treatment / Information used for model verification    |
| :-------------------------------- | :----------------------------------------------------------- |
| [Grahnen 1979](#5-references)     | Healthy subjects receiving a single intravenous dose of 100 mg and a single oral dose of 400 mg (tablet) |
| [Larsson 1982](#5-references)     | Peptic ulcer patients receiving a single intravenous dose of 200 mg |
| [Mihaly  1984](#5-references)     | Peptic ulcer patients receiving a single intravenous and a single oral dose of 200 mg |
| [Morgan 1983](#5-references)      | Peptic ulcer patients receiving a single intravenous dose of 200 mg (30 min infusion) |
| [Lebert 1981](#5-references)      | Healthy subjects receiving a single intravenous dose of 300 mg (2 min infusion) |
| [Walkenstein 1978](#5-references) | Healthy subjects receiving a single intravenous dose of 300 mg (2 min infusion) and a single oral dose of 300 mg (tablet) |
| [Kanto 1981](#5-references)       | Healthy subjects receiving a single oral dose of 200 mg      |
| [Burland 1975](#5-references)     | Healthy subjects receiving single oral doses of 200 mg solution and capsule |
| [Bodemar 1979](#5-references)     | Peptic ulcer patients receiving a single oral dose of 200 mg (tablet) |
| [Bodemar 1981](#5-references)     | Peptic ulcer patients receiving single oral doses of 800 mg and multiple oral doses of 200 and 400 mg |
| [Barbhaiya 1995](#5-references)   | Healthy subjects receiving multiple oral doses of 300 mg (tablet) |
| [Somogyi 1981](#5-references)     | Healthy subjects receiving a single oral dose of 400 mg (tablet) |
| [Tiseo 1998](#5-references)       | Healthy subjects receiving multiple oral doses of 800 mg (tablet) |

#### 2.2.2.3 Model update due to PK-Sim V10 conversion

As a consequence of updating the cimetidine PBPK model to PK-Sim version 10, the CYP3A4 K<sub>i</sub> value needed to be readjusted. For this purpose, AUC ratios of the following clinical DDI studies were used to inform K<sub>i</sub> in an additional parameter identification:

| Publication                      | Interaction of cimetidine with:                              |
| :------------------------------------- | :------------------------------|
| [Kienlen 1993](#5-references)    | Alfentanil |
| [Abernethy 1983](#5-references)  | Alprazolam and triazolam |
| [Elliott 1984](#5-references)    | Midazolam |
| [Fee 1987](#5-references)        | Midazolam |
| [Greenblatt 1986](#5-references) | Intravenous and oral midazolam |
| [Martinez 1999](#5-references)   | Midazolam |
| [Salonen 1986](#5-references)    | Midazolam |
| [Pourbaix 1985](#5-references)   | Triazolam. NOTE: The interaction of cimetidine with alprazolam of this publication was not used for parameterization due to very long simulation duration! |
| [Cox 1986](#5-references)        | Triazolam |
| [Friedman 1988](#5-references)   | Triazolam |

Similarly, MATE1 K<sub>i</sub> value was adjusted to reproduce the observed inhibition effect on metformin PK (https://github.com/Open-Systems-Pharmacology/Cimetidine-Metformin-DDI).

## Literature in-vitro kinetics (given)

These are the measured in-vitro values from the literature - the fitted
catalytic rates are NOT given; you must determine them.

| molecule | process | parameter | value | unit |
| --- | --- | --- | --- | --- |
| MATE1 | ActiveTransportSpecific_MM | Vmax | 0.0 | µmol/l/min |
| OAT3 | ActiveTransportSpecific_MM | Vmax | 0.0 | µmol/l/min |
| OAT3 | ActiveTransportSpecific_MM | Km | 149.0 | µmol/l |
| OCT1 | ActiveTransportSpecific_MM | Vmax | 0.0 | µmol/l/min |
| OCT1 | ActiveTransportSpecific_MM | Km | 2600.0 | µmol/l |
