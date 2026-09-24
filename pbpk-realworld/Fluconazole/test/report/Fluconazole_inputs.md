# Building and evaluation of a PBPK model for Fluconazole in healthy adults


> You are handed this report and the clinical data files. Extract the

> physicochemical and in-vitro inputs yourself, decide the model

> structure, and fit the parameters. The reference model's chosen

> methods and fitted values are withheld.


# Background

Fluconazole is an antifungal medication for the treatment of a number of different fungal infections. It is metabolized by UGT2B7 and renally excreted by glomerular filtration. About 70% of fluconazole is excreted unchanged in urine during 120 hours.  Fluconazole is a moderate inhibitor of CYP3A4 and can therefore serve as a perpetrator for CYP3A4 in DDI studies. 

The herein presented PBPK model of fluconazole has been developed using published pharmacokinetic clinical data by Palkama et al. ([Palkama 1998](#5-references)), Ripa et al. ([Ripa 1993](#5-references)), Shiba et al. ([Shiba 1990](#5-references)), Ahonen et al. ([Ahonen 1997](#5-references)), Brammer et al. ([Brammer 1990](#5-references)), Brammer et al. ([Brammer 1991](#5-references)), Thorpe et al. ([Thorpe 1990](#5-references) and Varhe et al. ([Varhe 1996](#5-references)). 
The model has then been evaluated by comparing observed data to simulations of both intravenously and orally administered fluconazole covering a dose range of 25 mg to 400 mg. Both single dose and multiple dose administration are represented in the development and evaluation data sets. 

The presented model includes the following features:

- metabolism by UGT2B7,
- inhibition of CYP3A4,
- renal clearance by glomerular filtration,
- oral absorption with dissolution rate assigned to the Weibull function for low, medium and high doses.

## 2.2 Data used<a id="data-used"></a>

### 2.2.1 In vitro and physicochemical data

A literature search was performed to collect available information on physicochemical properties of fluconazole. The obtained information from literature is summarized in the table below, and is used for model building.

| **Parameter**           | **Unit** | **Value** | Source                               | **Description**                                              |
| :---------------------- | -------- | --------- | ------------------------------------ | ------------------------------------------------------------ |
| MW                      | g/mol    | 306.28    | [Debruyne 1993](#5-references)       | Molecular weight                                             |
| pK<sub>a</sub>          |          | 2.03      | [Charoo 2014](#5-references)         | acid dissociation constant of conjugate acid; compound type: weak base |
| Solubility (pH)         | mg/mL    | 6.9 (7.4) | [Charoo 2014](#5-references)         | Aqueous Solubility                                           |
| fu                      | %        | 89        | Pfizer                               | Fraction unbound in plasma                                   |
|                         | %        | 88        | [Christofoletti 2016](#5-references) | Fraction unbound in plasma                                   |
| CL<sub>int</sub> UGT2B7 | 1/min    | 0.005     | [Watt 2018](#5-references)           | First order intrinsic clearance UGT2B7                       |
| P<sub>eff</sub>         | cm/s     | 0.000213  | GastroPlus                           | Specific intestinal permeability                             |
| Ki CYP3A4               | µmol/L   | 13.1     | [Sakaeda 2005](#5-references)           | Non-competitive inhibition of CYP3A4                         |
| GFR fraction            |          | 0.17      | [Watt 2018](#5-references)           | glomerular filtration fraction                               |

### 2.2.2 Clinical data

A literature search was performed to collect available clinical data on fluconazole in adults. 

The following publications were found in adults for model building:

| Publication                   | Arm / Treatment / Information used for model building        |
| :---------------------------- | :----------------------------------------------------------- |
| [Ahonen 1997](#5-references)  | Plasma PK profiles in healthy subjects after single dose oral administration of 400 mg fluconazole |
| [Brammer 1990](#5-references) | Plasma PK profiles in healthy subjects after multiple oral administration of 200, 300 and 400 mg fluconazole |
| [Brammer 1991](#5-references) | Plasma PK profiles in healthy subjects after single dose oral administration of 50 mg fluconazole |
| [Shiba 1990](#5-references)   | Plasma PK profiles and urine data in healthy subjects after single iv and oral administration of 25 and 50 mg fluconazole and single oral dose of 100 mg fluconazole |
| [Palkama 1998](#5-references) | Plasma PK profiles in healthy subjects after single dose iv and oral administration of 400 mg fluconazole |
| [Ripa 1993](#5-references)    | Plasma PK profiles in healthy subjects after single dose iv administration of 100 mg fluconazole |
| [Thorpe 1990](#5-references)  | Plasma PK profiles in healthy subjects after single dose oral administration of 100 mg fluconazole |
| [Varhe 1996](#5-references)   | Plasma PK profiles in healthy subjects after multiple oral administration of 100 and 200 mg fluconazole |

The following table shows the data from the excretion studies ([Shiba 1990](#5-references)) used for model building:

| Observer                                                     | Value |
| ------------------------------------------------------------ | ----- |
| Fraction excreted  to urine of unchanged fluconazole after iv administration 25 mg | 74%   |
| Fraction excreted  to urine of unchanged fluconazole after iv administration 25 mg | 72%   |
| Fraction excreted  to urine of unchanged fluconazole after oral administration 25 mg | 69%   |
| Fraction excreted  to urine of unchanged fluconazole after oral administration 50 mg | 67%   |
| Fraction excreted  to urine of unchanged fluconazole after oral administration 100 mg | 75%   |

The following dosing scenarios were simulated and compared to respective data for model verification:

| Scenario                                                     | Data reference                       |
| ------------------------------------------------------------ | ------------------------------------ |
| iv 400 mg (60 min)                          | [Ahonen 1997](#5-references) |
| iv 50 mg (2 min)                                 | [Brammer 1990](#5-references) |
| iv 25 mg (1 min) | [Shiba 1990](#5-references) |
| iv 800 mg (240 min) once daily for 14 days | [Sobue 2004](#5-references) |
| iv 100 mg (20 min)                             | [Yeates 1994](#5-references) |
| po 400 mg                              | [Palkama 1998](#5-references) |
| po 100 mg | [Shiba 1990](#5-references) |
| po 50 mg once daily for 4 days          | [Varhe 1996](#5-references) |
