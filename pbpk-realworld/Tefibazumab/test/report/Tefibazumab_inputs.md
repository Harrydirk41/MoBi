# Building and evaluation of a PBPK model for tefibazumab in healthy adults


> You are handed this report and the clinical data files. Extract the

> physicochemical and in-vitro inputs yourself, decide the model

> structure, and fit the parameters. The reference model's chosen

> methods and fitted values are withheld.


# Background

Tefibazumab is a humanized monoclonal antibody (IgG1) against the clumping factor A (ClfA) of *Staphylococcus aureus*. Tefibazumab shows a pharmacokinetic (PK) behavior which is typical for an antibody without endogenous target.  

The herein presented evaluation report evaluates the performance of the physiologically based pharmacokinetic (PBPK) model for tefibazumab in healthy adults. 

The presented Tefibazumab PBPK model as well as the respective evaluation plan and evaluation report are provided open-source (https://github.com/Open-Systems-Pharmacology/Tefibazumab-Model).

## 2.2 Data<a id="methods-data"></a>

### 2.2.1 In vitro / physico-chemical Data <a id="invitro-and-physico-chemical-data"></a>

A literature search was performed to collect available information on physicochemical properties of tefibazumab. The obtained information from literature is summarized in the table below. 

| **Parameter** | **Unit** | **Value** | Source                       | **Description**                                              |
| :------------ | -------- | --------- | ---------------------------- | ------------------------------------------------------------ |
| MW            | g/mol    | 150000    | [Lobo 2004](#5-references)   | Molecular weight                                             |
| r             | nm       | 5.34      | [Taylor 1984](#5-references) | Hydrodynamic solute radius                                   |
| Kd (FcRn)     | µM       | 0.63      | [Zhou 2003](#5-references)   | Dissociation constant for binding of a human IgG1 antibody to human FcRn at pH 6 |

### 2.2.2 PK Data <a id="PK-data"></a>

Published clinical PK data on tefibazumab in healthy adults were used.

| Publication                  | Description                                                  |
| :--------------------------- | :----------------------------------------------------------- |
| [Reilly 2005](#5-references) | The plasma concentration–time profiles after single dose 15 min i.v. infusion of 2, 5, 10, or 20 mg/kg body weight in healthy adults were used. |
