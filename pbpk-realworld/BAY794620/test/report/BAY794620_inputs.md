# Building and evaluation of a PBPK model for BAY 79-4620 in mice


> You are handed this report and the clinical data files. Extract the

> physicochemical and in-vitro inputs yourself, decide the model

> structure, and fit the parameters. The reference model's chosen

> methods and fitted values are withheld.


# Background

BAY 79-4620 is an antibody–drug conjugate consisting of a human IgG1 mAb directed against human carbonic anhydrase IX (CA IX) conjugated to the toxophore monomethylauristatin E (MMAE) via a valine-citrulline based linker ([Petrul 2012](#5-references)).

For BAY 79-4620, tissue concentration-time profiles for a large number of different mice tissues were measured and used together with pharmacokinetic (PK) data from 5 other compounds to identify unknown parameters during the development of the generic large molecule physiologically based pharmacokinetic (PBPK) model in PK-Sim ([Niederalt 2018](#5-references)). 

The herein presented evaluation report evaluates the performance of a PBPK model for BAY 79-4620 in xenograft mice for the PK data used for the development of the generic large molecule model in PK-Sim. A standard PK-Sim model was used without an additional tumor organ and without target mediated drug disposition effects from CA IX binding in the tumor - in contrast to the model which has been used in Ref. ([Niederalt 2018](#5-references)). 

The presented BAY 79-4620 PBPK model as well as the respective evaluation plan and evaluation report are provided open-source (https://github.com/Open-Systems-Pharmacology/BAY794620-Model).

## 2.2 Data<a id="methods-data"></a>

### 2.2.1 In vitro / physico-chemical Data <a id="invitro-and-physico-chemical-data"></a>

A literature search was performed to collect available information on physicochemical properties of BAY 79-4620. The obtained information from literature is summarized in the table below. 

| **Parameter** | **Unit** | **Value** | Source                       | **Description**                                              |
| :------------ | -------- | --------- | ---------------------------- | ------------------------------------------------------------ |
| MW            | g/mol    | 152000    | Bayer in-house data          | Molecular weight                                             |
| r             | nm       | 5.34      | [Taylor 1984](#5-references) | Hydrodynamic solute radius                                   |
| Kd (FcRn)     | µM       | 0.082     | [Zhou 2003](#5-references)   | Dissociation constant for binding of a humane IgG1 antibody to murine FcRn at pH 6 |

### 2.2.2 PK Data <a id="PK-data"></a>

The biodistribution data from mice for BAY 79-4620 were Bayer AG in house data taken from two studies:

| Data                        | Description                                                  |
| :-------------------------- | :----------------------------------------------------------- |
| Whole-body autoradiography  | Female nude mice (NMRI nu/nu), bearing HT-29 human colon carcinoma xenografts, were dosed intravenously with 1.25 mg/kg body weight of 125I-labeled BAY 79-4620. The distribution of total radioactivity in organs and tissues was determined by quantitative whole-body autoradiography after sacrificing the mice (two per time) at various time points after administration. |
| Wet-tissue dissection study | Female nude mice (NMRI nu/nu), bearing HT-29 human colon carcinoma xenografts, were dosed intravenously with 2 µCi (approx. 500 ng) of 125I-labeled BAY 79-4620. The distribution of total radioactivity in organs and tissues was determined after sacrificing the mice (three per time) and dissection of the organs at various time points after administration by determination of radioactivity using a gamma-counter. The concentrations were reported as percentage of dose / g tissue. These values were converted to concentrations in ng/ml assuming a density of 1 g/ml for all tissues except for bone for which a density of 1.5 g/ml was assumed (as in Ref. [Baxter 1994](#5-references)). |
