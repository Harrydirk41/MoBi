# Building and evaluation of a PBPK model for tizanidine in adults


> You are handed this report and the clinical data files. Extract the

> physicochemical and in-vitro inputs yourself, decide the model

> structure, and fit the parameters. The reference model's chosen

> methods and fitted values are withheld.


# Background

The presented PBPK model of tizanidine has been developed to be used in a PBPK Drug-Drug-Interactions (DDI) network with tizanidine as a substrate of CYP1A2.

Tizanidine is a centrally acting skeletal muscle relaxant generally used for the symptomatic treatment of acute painful muscle spasms and chronic spasticity resulting from diverse neurologic disorders ([Granfors 2004](#5-references)).

**Absorption**: After oral administration of tizanidine-HCl it is absorbed fast and completely with peak plasma concentrations reached within 1 hour. Administration of tizanidine together with food increases plasma concentrations. The influence of food on the concentration-time profile is also dependent on the formulation. A capsule given with food leads to a prolonged Tmax, with a Cmax slightly lower than when the drug is given without food. In contrast, giving the tablet with food leads to higher peak concentrations while Tmax remains unchanged.

**Distribution**: Approximately 30% is bound to plasma proteins. The concentration-time profile elicits a monophasic shape.

**Metabolism**: Over 95% of a dose is metabolized. The main enzyme system involved in the metabolism of tizanidine is CYP1A2.

**Excretion**: Only a minor part of the dose is recovered unchanged in the urine <5%.

## 2.2 Data<a id="data"></a>

### 2.2.1 In vitro and physico-chemical data

A literature search was performed to collect available information on physico-chemical properties of tizanidine ([Table 1](#table-1)).

| **Parameter**                   | **Unit**  | **Value**        | Source                            | **Description**                 |
| :------------------------------ | --------- | ---------------- | --------------------------------- | ------------------------------- |
| MW<sup>+</sup>                  | g/mol     | 253.711          | [DrugBank DB00697](#5-references) | Molecular weight                |
| pK<sub>a,base</sub><sup>+</sup> |           | 7.49             | [DrugBank DB00697](#5-references) | Acidic dissociation constant    |
| Solubility (pH)<sup>+</sup>     | mg/mL     | 0.133<br />(7.4) | [DrugBank DB00697](#5-references) | Aqueous Solubility              |
| fu<sup>+</sup>                  | %         | 70               | [SmPC tizanidine](#5-references)  | Fraction unbound in plasma      |
| Intrinsic CL                    | ml/min/kg | 17               | [Granfors 2004](#5-references)    | Predicted from microsomal assay |

**Table 1:**<a name="table-1"></a> Physico-chemical and *in-vitro* metabolization properties of tizanidine extracted from literature. *<sup>+</sup>: Value used in final model*

### 2.2.2 Clinical data

A literature search was performed to collect available clinical data on tizanidine ([Table 2](#table-2)).

| **Source**           | **Dose [mg]/**  **Schedule \*** | **Pop.**     | Age [yrs] (mean or range) | Weight [kg] (mean or range) | **Sex** | **N** | **Form.** | Fasted or Fed | **Comment**                       |
| -------------------- | ------------------------------- | ------------ | ------- | ----- | --------- | --------------------------------- | --------------------------------- | --------------------------------- | -------------------- |
| [Momo 2010](#5-references)<sup>+</sup> | 2               | HV                | 29                      | 70                        | m    | 12   | Tablet         | Fed       |              |
| [Granfors 2004](#5-references)<sup>+</sup> | 4                               | HV       | 21-31                     | 65-83                       | m       | 10    | Tablet           | Fasted        |                                          |
| [Schellenberger 1999](#5-references) | 4 t.i.d.                        | HV       | 19-37                     | 70-97                       | m       | 12    | Tablet           | Fasted        |                                          |
| [Henney 2007](#5-references)<sup>+</sup> | 4                               | HV       | 26                        | 71                          | m12/f6  | 18    | Tablet /Capsule  | Fed           |                                          |
| [Backman 2008](#5-references)<sup>+</sup> | 4                               | HV       | 23                        | 78                          | m       | 38    | Tablet           | Fasted        | Male/female-non-smokers and male-smokers |
| [Backman 2006](#5-references)<sup>+</sup> | 4                               | HV       | 21                        | 71                          | m6/f4   | 10    | Tablet           | Fasted        | Only control group used                  |
| [Shah 2006](#5-references)           | 8                               | HV       | 18-52                     | 46-102                      | m54/f42 | 96    | Tablet / Capsule | Fed/   Fasted |                        |
| [Henney 2008](#5-references)         | 6                               | HV       | 18-39                     | NA                          | m19/f9  | 27    | Capsule          | Fasted        | |
| [Tse 1987](#5-references)            | 4 t.i.d.                        | HV       | 21-48                     | 57-86                       | m       | 6     | Tablet           | Fasted        | |
| [Al-Ghazawi 2013](#5-references)<sup>+</sup> | 4                               | HV       | 28                        | 75                          | m       | 36    | Tablet           | Fasted        | |

**Table 2:**<a name="table-2"></a> Literature sources of clinical concentration data of tizanidine used for model development and validation.  *\*: single dose unless otherwise specified; EM: extensive metabolizers;<sup>+</sup>: Data used for final parameter identification*
