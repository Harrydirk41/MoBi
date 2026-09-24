# Building and evaluation of a PBPK model for montelukast in adults


> You are handed this report and the clinical data files. Extract the

> physicochemical and in-vitro inputs yourself, decide the model

> structure, and fit the parameters. The reference model's chosen

> methods and fitted values are withheld.


# Background

The presented model building and evaluation report evaluates the performance of a physiology-based  pharmacokinetic (PBPK) model for montelukast in adults.

Montelukast is a selective and orally active leukotriene receptor antagonist that inhibits the cysteinyl leukotriene (CysLT) receptor 1, used in the maintenance treatment of asthma. Montelukast is mainly metabolized by CYP2C8 (72%) ([Marzolini 2017](#5-references)).  Montelukast is a strongly lipophilic drug. The final lipophilicity was estimated to be lower than the reported values, as with lipophilicity values above 3-4 log units the drug already reached maximal permeability levels. The final montelukast model applies metabolism by CYP2C8, and to a minor extend involved clearance by the enzymes CYP3A4/5 (16% ), CYP2C9 (12%) and glomerular filtration ([Marzolini 2017, Filppula 2011, Zhou 2017](#5-references)) and adequately described the pharmacokinetics of montelukast in adults.

The montelukast model is a whole-body PBPK model, allowing for dynamic translation between individuals. The montelukast report demonstrates the level of confidence in the montelukast PBPK model built with the OSP suite with regard to reliable predictions of montelukast pharmacokinetics (PK) in adults during model-informed drug development.

## 2.2 Data used<a id="data-used"></a>

### 2.2.1 In vitro / physicochemical data

A literature search was performed to collect available information on physicochemical properties of montelukast. The obtained information from literature is summarized in the table below and is used for model building.

| **Parameter**   | **Unit**    | **Value (reference)** | **Description**                                  |
| :-------------- | ----------- | ----------------------------------- | ------------------------------------------------ |
| MW              | g/mol       | 586.2 ([Marzolini 2017](#5-references)) | Molecular weight                                 |
| pKa             |             | 4.4 ([Marzolini 2017](#5-references)) | Acid dissociation constant                   |
| Solubility (pH) | mg/mL       | 8.2E-06 (7) ([Drugbank](#5-references)) | Solubility                                       |
| fu              |             | 0.0018 ([Marzolini 2017](#5-references)) | Fraction unbound                                 |
| fe**   |          | <0.002 ([Marzolini 2017](#5-references)) | fraction of dose excreted unchanged in urine |
| CYP3A4-CLint | µl/min/pmol | 1.8 ([Marzolini 2017](#5-references)) | Cytochrome-P450 3A4 mediated intrinsic clearance |
| CYP3A5-CLint | µl/min/pmol | 1.8 ([Marzolini 2017](#5-references)) | Cytochrome-P450 3A5 mediated intrinsic clearance |
| CYP2C8-CLint | µl/min/pmol | 3.6 ([Marzolini 2017](#5-references)) | Cytochrome-P450 2C8 mediated intrinsic clearance |
| CYP2C9-CLint | µl/min/pmol | 0.48 ([Marzolini 2017](#5-references)) | Cytochrome-P450 2C9 mediated intrinsic clearance |

** fe was matched by modeling unchanged renal excretion in PK-Sim as glomerular filtration (GF)

### 2.2.2 Clinical data

A literature search was performed to collect available clinical data on montelukast PK in adults. 

The following publications were found in adults for model building and evaluation:

| Publication                       | Study description                                            |
| :-------------------------------- | :----------------------------------------------------------- |
| [Cheng 1996](#5-references) | Pharmacokinetics, bioavailability, and safety of montelukast sodium (MK-0476) in healthy males and females |
| [Fey 2014](#5-references) | Bioequivalence of two formulations of montelukast sodium 4 mg oral granules in healthy adults |
| [Knorr 2000](#5-references) | Montelukast adult (10-mg film-coated tablet) and pediatric (5-mg chewable tablet) dose selections |
| [Zhao 1997](#5-references) | Pharmacokinetics and bioavailability of montelukast sodium (MK-0476) in healthy young and elderly volunteers |
