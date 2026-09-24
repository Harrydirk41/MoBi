# Building and evaluation of a PBPK model for inulin in rats


> You are handed this report and the clinical data files. Extract the

> physicochemical and in-vitro inputs yourself, decide the model

> structure, and fit the parameters. The reference model's chosen

> methods and fitted values are withheld.


# Background

Inulin is a highly hydrophilic polysaccharide which does not distribute into cells and is cleared via glomerular filtration. 

Inulin has a considerably smaller solute radius than the proteins which had been used to develop the generic large molecule physiologically based pharmacokinetic (PBPK) model in PK-Sim ([Niederalt 2018](#5-references)).   

The herein presented evaluation report evaluates the performance of the PBPK model for inulin in rats using the large molecule model in PK-Sim.

The presented inulin PBPK model as well as the respective evaluation plan and evaluation report are provided open-source (https://github.com/Open-Systems-Pharmacology/Inulin-Model).

## 2.2 Data<a id="methods-data"></a>

### 2.2.1 In vitro / physico-chemical Data <a id="invitro-and-physico-chemical-data"></a>

A literature search was performed to collect available information on physicochemical properties of Inulin. The obtained information from literature is summarized in the table below. 

| **Parameter** | **Unit** | **Value** | Source                           | **Description**                                              |
| :------------ | -------- | --------- | -------------------------------- | ------------------------------------------------------------ |
| MW            | g/mol    | 5000-5500 | [Ohno 1978](#5-references)       | Molecular weight                                             |
| r             | nm       | 1.39      | [Ghandehari 1997](#5-references) | Hydrodynamic solute radius                                   |
| Kd (FcRn)     | µM       | 999,999   |                                  | Dissociation constant for binding to FcRn. High value representing no FcRn binding. |

### 2.2.2 PK Data <a id="PK-data"></a>

Published plasma and tissue PK data on inulin in rats were used.

| Publication                 | Description                                                  |
| :-------------------------- | :----------------------------------------------------------- |
| [Tsuji 1983](#5-references) | Plasma and tissue concentrations after i.v. application of 20 and 200 mg/kg  inulin in rats (for 200 mg/kg plasma only). |
