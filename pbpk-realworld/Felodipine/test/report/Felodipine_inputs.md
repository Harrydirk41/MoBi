# Building and evaluation of a PBPK model for Felodipine in healthy adults


> You are handed this report and the clinical data files. Extract the

> physicochemical and in-vitro inputs yourself, decide the model

> structure, and fit the parameters. The reference model's chosen

> methods and fitted values are withheld.


# Background

Felodipine is a calcium-channel blocker, indicated for angina pectoris and arterial hypertension. It is mostly metabolized by CYP3A4 making it a sensitive probe and victim drug for the investigation of CYP3A4 activity *in vivo*. It is a BCS class II compound. Felodipine shows substantial first pass metabolism resulting in a bioavailability of 15%. 

The model has been developed and evaluated by comparing observed data to simulations of a large number of clinical studies covering a dose range of 1.5 mg to 10 mg after intravenous and oral administrations. Furthermore, it has been evaluated within a CYP3A4 DDI modeling network as a victim drug. 

Model features include:

- metabolism by CYP3A4
- metabolism by an unknown enzyme *via* unspecific hepatic clearance
- a decrease in the permeability between the intracellular and interstitial space (model parameters `P (intracellular->interstitial)` and `P (interstitial->intracellular)`) in intestinal mucosa to optimize quantitatively the extent of gut wall metabolism

## 2.2 Data<a id="data"></a>

### 2.2.1 In vitro and physicochemical data

A literature search was performed to collect available information on physicochemical properties of felodipine. The obtained information from literature is summarized in the table below, and is used for model building.

| **Parameter**                         | **Unit**                   | **Value**        | Source                            | **Description**                                              |
| :------------------------------------ | -------------------------- | ---------------- | --------------------------------- | ------------------------------------------------------------ |
| MW                                    | g/mol                      | 384.254          | [DrugBank DB01023](#5-references) | Molecular weight                                             |
| pK<sub>a1</sub>                       |                            | n.a.             | [Alskar 2018](#5-references)      | Acid dissociation constant of conjugate acid                 |
| pK<sub>a1</sub>                       |                            | 5.07             | [Pandey 2013](#5-references)      | Acid dissociation constant of conjugate acid; compound type: acid |
| Solubility (pH)                       | mg/L                       | 14.3<br />(7.1)  | [Takano 2016](#5-references)      | Aqueous Solubility                                           |
|                                       |                            | 1                | [Scholz 2002](#5-references)      | Aqueous Solubility                                           |
|                                       |                            | 53               | [Söderlind 2010](#5-references)   | Solubility in fasted state simulated intestinal fluid I      |
|                                       |                            | 12               | [Söderlind 2010](#5-references)   | Solubility in fasted state simulated intestinal fluid II     |
|                                       |                            | 14               | [Söderlind 2010](#5-references)   | Solubility in fasted human intestinal fluid                  |
|                                       |                            | 15<br />(7.5)    | [Persson 2005](#5-references)     | Solubility in fasted human intestinal fluid                  |
|                                       |                            | 413<br />(6.1)   | [Persson 2005](#5-references)     | Solubility in fed human intestinal fluid                     |
|                                       |                            | 191              | [Persson 2005](#5-references)     | Solubility in fed state simulated intestinal fluid           |
|                                       |                            | 77<br />(6.35)   | [Scholz 2002](#5-references)      | Solubility in fasted state chyme                             |
|                                       |                            | 56<br />(4.93)   | [Scholz 2002](#5-references)      | Solubility in fed state chyme                                |
| fu                                    | %                          | 0.36             | [Soons 1993](#5-references)       | Fraction unbound in plasma                                   |
|                                       | %                          | 0.36             | [Ushimura 2010](#5-references)    | Fraction unbound in plasma                                   |
| V<sub>max</sub>, K<sub>m</sub> CYP3A  | pmol/mg/min,<br />µmol/L   | 1630<br />2.81   | [Walsky 2004](#5-references)      | CYP3A liver microsomes Michaelis-Menten kinetics             |
| V<sub>max</sub>, K<sub>m</sub> CYP3A  | pmol/mg/min,<br />µmol/L   | 240<br />6.9     | [Bu 2006](#5-references)          | CYP3A liver microsomes Michaelis-Menten kinetics             |   
| V<sub>max</sub>, K<sub>m</sub> CYP3A4 | pmol/mg/min,<br />µmol/L   | 36.8<br />0.938  | [Walsky 2004](#5-references)      | Recombinant CYP3A4 Michaelis-Menten kinetics                 |
| V<sub>max</sub>, K<sub>m</sub> CYP3A5 | pmol/mg/min,<br />µmol/L   | 24.2<br />1.41   | [Walsky 2004](#5-references)      | Recombinant CYP3A5 Michaelis-Menten kinetics                 |              

### 2.2.2 Clinical data

A literature search was performed to collect available clinical data on felodipine in adults. 

The following publications were found in adults for model building:

| Publication                            | Arm / Treatment / Information used for model building                                                                                                       |
| :------------------------------------- | :---------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [Lundahl 1997](#5-references)          | Plasma PK profiles in healthy subjects with single dose administrations of a felodipine 1.5 mg intravenous infusion                                         |
| [Edgar 1987](#5-references)            | Plasma PK profiles in healthy subjects with single dose administrations of a felodipine:<br />- 10 mg oral solution <br />- 10 mg extended release tablet <br />- 10 mg immediate release tablet   |
| [Blychert 1990](#5-references)         | Plasma PK profiles in healthy subjects with multiple dose administrations of a felodipine:<br />- 10 mg oral solution <br />- 10 mg extended release tablet |
| [Goosen 2004](#5-references)           | Plasma PK profiles in healthy subjects with single dose administrations of a felodipine 5 mg extended release tablet                                        |
| [Jalava 1997](#5-references)           | Plasma PK profiles in healthy subjects with single dose administrations of a felodipine 5 mg extended release tablet <br />- alone (control) <br />- with itraconazole (treatment)|
| [Bailey 1996](#5-references)           | Plasma PK profiles in healthy subjects with single dose administrations of a felodipine 10 mg extended release tablet                                       |
| [Gelal 2005](#5-references)            | Plasma PK profiles in healthy subjects with single dose administrations of a felodipine 10 mg extended release tablet                                       |
| [Bailey 2003](#5-references)           | Plasma PK profiles in healthy subjects with single dose administrations of a felodipine 10 mg extended release tablet                                       |
| [Bailey 1993](#5-references)           | Plasma PK profiles in healthy subjects with single dose administrations of a felodipine 5 mg immediate release tablet                                       |
| [Edgar 1992](#5-references)            | Plasma PK profiles in healthy subjects with single dose administrations of a felodipine 5 mg immediate release tablet                                       |
| [Blychert 1990](#5-references)         | Plasma PK profiles in healthy subjects with multiple dose administrations of a felodipine:<br />- 10 mg oral solution <br />- 10 mg extended release tablet |
| [Lundahl 1998](#5-references)          | Plasma PK profiles in healthy subjects with multiple dose administrations of a felodipine 10 mg extended release tablet                                     |
| [Bailey 1995](#5-references)           | Plasma PK profiles in healthy subjects with multiple dose administrations of a felodipine 10 mg extended release tablet                                     |
| [Aberg 1997](#5-references)            | Plasma PK profiles in healthy subjects with multiple dose administrations of a felodipine 10 mg extended release tablet                                     |

The following dosing scenarios were simulated and compared to respective data for model verification:

| Scenario                                                     | Data reference                       |
| ------------------------------------------------------------ | ------------------------------------ |
| po 5 mg single dose (extended release tablet)                | [Dresser 2000](#5-references)        |
| po 10 mg single dose (extended release tablet)               | [Dresser 2017](#5-references)        |
|                                                              | [Dresser 2002](#5-references)        |
|                                                              | [Bailey 2000](#5-references)         |
|                                                              | [Bailey 1998](#5-references)         |
|                                                              | [Lundahl 1997](#5-references)        |
|                                                              | [Madsen 1996](#5-references)         |
| po 2.5 / 5 mg once daily (extended release tablet)           | [Dresser 2000](#5-references)        |
| po 10 mg once daily (immediate release tablet)               | [Blychert 1990](#5-references)       |
| po 10 mg twice daily (immediate release tablet)              | [Blychert 1990](#5-references)       |
|                                                              | [Smith 1987](#5-references)          |
| po 10 mg three times daily (immediate release tablet)        | [Bratel 1989](#5-references)         |

## Literature in-vitro kinetics (given)

These are the measured in-vitro values from the literature - the fitted
catalytic rates are NOT given; you must determine them.

| molecule | process | parameter | value | unit |
| --- | --- | --- | --- | --- |
| CYP3A4 | MetabolizationLiverMicrosomes_MM | In vitro Vmax for liver microsomes | 1630.0 | pmol/min/mg mic. protein |
| CYP3A4 | MetabolizationLiverMicrosomes_MM | Km | 2.81 | µmol/l |
