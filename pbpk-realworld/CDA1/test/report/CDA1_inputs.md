# Building and evaluation of a PBPK model for antibody CDA1 in healthy adults


> You are handed this report and the clinical data files. Extract the

> physicochemical and in-vitro inputs yourself, decide the model

> structure, and fit the parameters. The reference model's chosen

> methods and fitted values are withheld.


# Background

CDA1 is a human monoclonal antibody (IgG1) against the toxin A of *Clostridium difficile*.

CDA1 shows a pharmacokinetic behavior which is typical for an antibody without endogenous target. The plasma concentration–time profiles after i.v. infusion of 5, 10 and 20 mg/kg CDA1 in healthy adults ([Taylor 2008](#5-references)) were used together with pharmacokinetic (PK) data from 5 other compounds to identify unknown parameters during the development of the generic large molecule physiologically based pharmacokinetic (PBPK) model in PK-Sim ([Niederalt 2018](#5-references)). 

The herein presented evaluation report evaluates the performance of the PBPK model for CDA1 in healthy adults for the PK data used for the development of the generic large molecule model in PK-Sim.

The presented CDA1 PBPK model as well as the respective evaluation plan and evaluation report are provided open-source (https://github.com/Open-Systems-Pharmacology/CDA1-Model).

## 2.2 Data<a id="methods-data"></a>

### 2.2.1 In vitro / physico-chemical Data <a id="invitro-and-physico-chemical-data"></a>

A literature search was performed to collect available information on physicochemical properties of CDA1. The obtained information from literature is summarized in the table below. 

| **Parameter** | **Unit** | **Value** | Source                       | **Description**                                              |
| :------------ | -------- | --------- | ---------------------------- | ------------------------------------------------------------ |
| MW            | g/mol    | 150000    | [Lobo 2004](#5-references)   | Molecular weight                                             |
| r             | nm       | 5.34      | [Taylor 1984](#5-references) | Hydrodynamic solute radius                                   |
| Kd (FcRn)     | µM       | 0.63      | [Zhou 2003](#5-references)   | Dissociation constant for binding of a human IgG1 antibody to human FcRn at pH 6 |

### 2.2.2 PK Data <a id="PK-data"></a>

Published clinical PK data on CDA1 in healthy adults were used.

| Publication                 | Description                                                  |
| :-------------------------- | :----------------------------------------------------------- |
| [Taylor2008](#5-references) | The plasma concentration–time profiles after single i.v. infusion of 5, 10 and 20 mg/kg CDA1 in healthy adults were used. The data for the dosages 0.3 and 1 mg/kg were not used since the PK data could not be read with sufficient accuracy from the published figure. |
