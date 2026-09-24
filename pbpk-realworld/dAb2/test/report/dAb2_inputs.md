# Building and evaluation of a PBPK model for domain antibody dAb2 in mice


> You are handed this report and the clinical data files. Extract the

> physicochemical and in-vitro inputs yourself, decide the model

> structure, and fit the parameters. The reference model's chosen

> methods and fitted values are withheld.


# Background

The dAb2 domain antibody is a fusion protein consisting of a VH (heavy chain) and a Vk (light chain) antibody fragment without known binding affinity which was used to a develop a PBPK model  ([Sepp2015](#5-references)). 

Since the dAb2 is smaller than antibodies, the PK data (blood and tissue concentration–time profiles in mice)  ([Sepp2015](#5-references)) were also used  together with pharmacokinetic (PK) data from 5 other compounds to identify unknown parameters during the development of the generic large molecule physiologically based pharmacokinetic (PBPK) model in PK-Sim ([Niederalt 2018](#5-references)). 

The herein presented evaluation report evaluates the performance of the PBPK model for dAb2 in mice for the PK data used for the development of the generic large molecule model in PK-Sim.

The presented dAb2 PBPK model as well as the respective evaluation plan and evaluation report are provided open-source (https://github.com/Open-Systems-Pharmacology/dAb2-Model)

## 2.2 Data<a id="methods-data"></a>

### 2.2.1 In vitro / physico-chemical Data <a id="invitro-and-physico-chemical-data"></a>

A literature search was performed to collect available information on physicochemical properties of dAb2. The obtained information from literature is summarized in the table below. 

| **Parameter** | **Unit** | **Value** | Source                    | **Description**                                              |
| :------------ | -------- | --------- | ------------------------- | ------------------------------------------------------------ |
| MW            | g/mol    | 25600     | [Sepp2015](#5-references) | Molecular weight                                             |
| r             | nm       | 2.43      | calculated   from MW      | Hydrodynamic solute radius. Calculated by empirical equation given in [Niederalt2018](#5-references), supplemental material |
| Kd (FcRn)     | µM       | 999,999   |                           | high value representing no FcRn binding                      |

### 2.2.2 PK Data <a id="PK-data"></a>

Published plasma and tissue PK data on dAb2 in mice were used.

| Publication               | Description                                                  |
| :------------------------ | :----------------------------------------------------------- |
| [Sepp2015](#5-references) | Plasma and tissue PK data after an intravenous dose of dose of 10 mg/kg in mice. Tissue concentrations were analyzed using quantitative whole-body autoradiography. The concentrations were reported as percentage of injected dose / g tissue. These values were converted to concentrations in µg/ml assuming a density of 1 g/ml for all tissues except for bone for which a density of 1.5 g/ml was assumed (as in Ref. [Baxter 1994](#5-references)). Furthermore, a body weight of 29 g (i.e. a dose of  290 µg) was assumed for unit conversion of the experimental concentrations (body weight range reported: 26-33 g). |
