# Building and evaluation of a PBPK model for antibody MEDI-524-YTE in cynomolgus monkeys


> You are handed this report and the clinical data files. Extract the

> physicochemical and in-vitro inputs yourself, decide the model

> structure, and fit the parameters. The reference model's chosen

> methods and fitted values are withheld.


# Background

MEDI-524-YTE is variant of the humanized monoclonal IgG1 antibody MEDI-524 against the respiratory syncytial virus (RSV). The triple YTE mutation introduced into the Fc region led to an increased affinity to FcRn and consequently to an increased plasma half-life ([Dall’Acqua2006](#5-references)). 

The plasma concentration–time profile after intravenous application of a 30 mg/kg dose in cynomolgus monkeys ([Dall’Acqua2006](#5-references)) were used together with pharmacokinetic (PK) data from 5 other compounds to identify unknown parameters during the development of the generic large molecule physiologically based pharmacokinetic (PBPK) model in PK-Sim ([Niederalt 2018](#5-references)). 

The herein presented evaluation report evaluates the performance of the PBPK model for MEDI-524-YTE in cynomolgus monkeys for the PK data used for the development of the generic large molecule model in PK-Sim.

The presented MEDI-524-YTE PBPK model as well as the respective evaluation plan and evaluation report are provided open-source (https://github.com/Open-Systems-Pharmacology/MEDI524YTE-Model).

## 2.2 Data<a id="methods-data"></a>

### 2.2.1 In vitro / physico-chemical Data <a id="invitro-and-physico-chemical-data"></a>

A literature search was performed to collect available information on physicochemical properties of MEDI-524-YTE. The obtained information from literature is summarized in the table below. 

| **Parameter** | **Unit** | **Value** | Source                           | **Description**                                              |
| :------------ | -------- | --------- | -------------------------------- | ------------------------------------------------------------ |
| MW            | g/mol    | 150000    | [Lobo 2004](#5-references)       | Molecular weight                                             |
| r             | nm       | 5.34      | [Taylor 1984](#5-references)     | Hydrodynamic solute radius                                   |
| Kd (FcRn)     | µM       | 0.134     | [Dall'Acqua 2006](#5-references) | Dissociation constant for binding to cynomolgus monkey FcRn for the Fc variant MEDI-524-YTE (pH 6) |

### 2.2.2 PK Data <a id="PK-data"></a>

Published PK data on MEDI-524-YTE in cynomolgus monkeys were used.

| Publication                      | Description                                                  |
| :------------------------------- | :----------------------------------------------------------- |
| [Dall'Acqua 2006](#5-references) | The plasma concentration–time profiles after single i.v. infusion of 30 mg/kg MEDI-524-YTE in in cynomolgus monkeys were used. |
