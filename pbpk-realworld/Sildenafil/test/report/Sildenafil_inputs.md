# Building and evaluation of a PBPK model for Sildenafil in healthy adults


> You are handed this report and the clinical data files. Extract the

> physicochemical and in-vitro inputs yourself, decide the model

> structure, and fit the parameters. The reference model's chosen

> methods and fitted values are withheld.


# Background

Sildenafil is a cGMP-specific phosphodiesterase 5 inhibitor, indicated for erectile dysfunction and pulmonary arterial hypertension. It is mostly metabolized by CYP3A4 making it a sensitive probe and victim drug for the investigation of CYP3A4 activity *in vivo*. Other CYPs are involved in sildenafil metabolism: CYP2C9 and CYP2C19. It is a BCS class II compound. Sildenafil shows substantial first pass metabolism resulting in a bioavailability of 40%. 

The model has been developed and evaluated by comparing observed data to simulations of a large number of clinical studies covering a dose range of 20 mg to 100 mg after intravenous and oral administrations. Furthermore, it has been evaluated within a CYP3A4 DDI modeling network as a victim drug. 

Model features include:

- metabolism by CYP3A4
- metabolism by CYP2C9
- metabolism by CYP2C19
- a decrease in the permeability between the intracellular and interstitial space (model parameters `P (intracellular->interstitial)` and `P (interstitial->intracellular)`) in intestinal mucosa to optimize quantitatively the extent of gut wall metabolism

## 2.2 Data<a id="data"></a>

### 2.2.1 In vitro and physicochemical data

A literature search was performed to collect available information on physicochemical properties of sildenafil. The obtained information from literature is summarized in the table below, and is used for model building.

| **Parameter**                         | **Unit**                   | **Value**        | Source                            | **Description**                                              |
| :------------------------------------ | -------------------------- | ---------------- | --------------------------------- | ------------------------------------------------------------ |
| MW                                    | g/mol                      | 474.576          | [DrugBank DB00203](#5-references) | Molecular weight                                             |
| pK<sub>a1</sub>                       |                            | 5.97             | [Salerno 2021](#5-references)     | Acid dissociation constant of conjugate acid; compound type: basic |
| pK<sub>a1</sub>                       |                            | 6.78             | [Gobry 2000](#5-references)       | Acid dissociation constant of conjugate acid; compound type: ampholyte |
| pK<sub>a2</sub>                       |                            | 9.12             | [Gobry 2000](#5-references)       | Acid dissociation constant of conjugate acid; compound type: ampholyte |
| Solubility (pH)                       | mg/mL                      | 0.025<br />(7.1) | [Takano 2016](#5-references)      | Aqueous Solubility                                           |
|                                       |                            | 3.5              | [Salerno 2021](#5-references)     | Aqueous Solubility                                           |
|                                       |                            | 3.5              | [DrugBank DB00203](#5-references) | Aqueous Solubility                                           |
|                                       |                            | 4.1              | [Jung 2011](#5-references)        | Aqueous Solubility                                           |
|                                       |                            | 3.965<br />(3)   | [Wang 2008](#5-references)        | Aqueous Solubility                                           |
|                                       |                            | 7.077<br />(4)   | [Wang 2008](#5-references)        | Aqueous Solubility                                           |
|                                       |                            | 2.068<br />(5)   | [Wang 2008](#5-references)        | Aqueous Solubility                                           |
|                                       |                            | 0.114<br />(6)   | [Wang 2008](#5-references)        | Aqueous Solubility                                           |
|                                       |                            | 0.025<br />(7)   | [Wang 2008](#5-references)        | Aqueous Solubility                                           |
|                                       |                            | 0.027<br />(8)   | [Wang 2008](#5-references)        | Aqueous Solubility                                           |
|                                       |                            | 0.04<br />(9)    | [Wang 2008](#5-references)        | Aqueous Solubility                                           |
|                                       |                            | 0.103<br />(10)  | [Wang 2008](#5-references)        | Aqueous Solubility                                           |
|                                       |                            | 0.322<br />(11)  | [Wang 2008](#5-references)        | Aqueous Solubility                                           |
| fu                                    | %                          | 4                | [Walker 1999](#5-references)      | Fraction unbound in plasma (α1-acid glycoprotein)            |
|                                       | %                          | 4.3              | [Muirhead 2002b](#5-references)    | Fraction unbound in plasma (α1-acid glycoprotein)            |
|                                       | %                          | 2.7              | [Muirhead 2002b](#5-references)    | Fraction unbound in plasma (α1-acid glycoprotein)            |
|                                       | %                          | 3.46             | [Muirhead 2002b](#5-references)    | Fraction unbound in plasma (α1-acid glycoprotein)            |
| V<sub>max</sub>, K<sub>m</sub> CYP3A4 | pmol/min/pmol P450,<br />µmol/L | 78.6<br />4.34   | [Takano 2016](#5-references)      | Recombinant CYP3A4 Michaelis-Menten kinetics                 |
| V<sub>max</sub>, K<sub>m</sub> CYP3A4 | relative units,<br />µmol/L| 1.9<br />23.10   | [Warrington 2000](#5-references)  | Recombinant CYP3A4 Michaelis-Menten kinetics                 |   
| V<sub>max</sub>, K<sub>m</sub> CYP2C9 | relative units,<br />µmol/L| 0.2<br />9.60    | [Warrington 2000](#5-references)  | Recombinant CYP3A4 Michaelis-Menten kinetics                 |
| V<sub>max</sub>, K<sub>m</sub> CYP2C19| relative units,<br />µmol/L| 0.02<br />23.10  | [Warrington 2000](#5-references)  | Recombinant CYP3A4 Michaelis-Menten kinetics                 |              

### 2.2.2 Clinical data

A literature search was performed to collect available clinical data on sildenafil in adults. 

The following publications were found in adults for model building:

| Publication                            | Arm / Treatment / Information used for model building                                                                                                |
| :------------------------------------- | :--------------------------------------------------------------------------------------------------------------------------------------------------- |
| [Muirhead 2002a](#5-references)         | Plasma PK profiles in healthy subjects with single dose administrations of a sildenafil 25 mg intravenous infusion                                   |
| [Nichols 2002](#5-references)          | Plasma PK profiles in healthy subjects with single dose administrations of a sildenafil:<br />- 50 mg intravenous infusion <br />- 100mg oral tablet |
| [FDA 2009](#5-references) | Plasma PK profiles in healthy subjects with single dose administrations of a sildenafil:<br />- 20 mg intravenous infusion <br />- 40 mg intravenous infusion <br />- 80 mg intravenous infusion|
| [Walker 1999](#5-references)           | Plasma PK profiles in healthy subjects with single dose administrations of a sildenafil 50 mg oral solution                                          |
| [Spence 2008](#5-references)           | Plasma PK profiles in healthy subjects with single dose administrations of a sildenafil 20 mg tablet                                                 |
| [Lee 2021](#5-references)              | Plasma PK profiles in healthy subjects with single dose administrations of a sildenafil:<br />- 25 mg tablet (in the absence of itraconazole)<br />- 25 mg tablet (in the absence of clarithromycin)|
| [Abdelkawy 2016](#5-references)        | Plasma PK profiles in healthy subjects with single dose administrations of a sildenafil 50 mg tablet                                                 |
| [Gillen 2017](#5-references)           | Plasma PK profiles in healthy subjects with single dose administrations of a sildenafil 50 mg tablet (Panel 1)                                       |
| [Jetter 2002](#5-references)           | Plasma PK profiles in healthy subjects with single dose administrations of a sildenafil 50 mg tablet                                                 |
| [Murtadha 2021](#5-references)         | Plasma PK profiles in healthy subjects with single dose administrations of a sildenafil 50 mg tablet (non-smoker group)                              |
| [Wilner 2002](#5-references)           | Plasma PK profiles in healthy subjects with single dose administrations of a sildenafil 50 mg tablet (study I)                                       |

The following dosing scenarios were simulated and compared to respective data for model verification:

| Scenario                                                     | Data reference                       |
| ------------------------------------------------------------ | ------------------------------------ |
| po SD 50mg                                                   | [Al-Ghazawi 2010](#5-references)     |
|                                                              | [Hedaya 2006](#5-references)         |
|                                                              | [Wilner 2002](#5-references)         |
|                                                              | [Gillen 2017](#5-references)         |
| po SD 100mg                                                  | [Muirhead 2000](#5-references)       |
| po MD 20/80 mg                                               | [Burgess 2008](#5-references)        |
| po MD 20 mg                                                  | [Gotzkowsky 2013](#5-references)     |

## Literature in-vitro kinetics (given)

These are the measured in-vitro values from the literature - the fitted
catalytic rates are NOT given; you must determine them.

| molecule | process | parameter | value | unit |
| --- | --- | --- | --- | --- |
| CYP2C19 | rCYP450_MM | In vitro Vmax/recombinant enzyme | 0.01009105 | pmol/min/pmol rec. enzyme |
| CYP2C19 | rCYP450_MM | Km | 23.1 | µmol/l |
| CYP2C9 | rCYP450_MM | In vitro Vmax/recombinant enzyme | 0.02 | pmol/min/pmol rec. enzyme |
| CYP2C9 | rCYP450_MM | Km | 9.6 | µmol/l |
| CYP3A4 | rCYP450_MM | In vitro Vmax/recombinant enzyme | 0.16918611 | pmol/min/pmol rec. enzyme |
| CYP3A4 | rCYP450_MM | Km | 23.1 | µmol/l |
