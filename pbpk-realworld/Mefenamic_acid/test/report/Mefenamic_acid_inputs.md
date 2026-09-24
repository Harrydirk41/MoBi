# Building and Evaluation of a PBPK Model for Mefenamic Acid in Adults


> You are handed this report and the clinical data files. Extract the

> physicochemical and in-vitro inputs yourself, decide the model

> structure, and fit the parameters. The reference model's chosen

> methods and fitted values are withheld.


# Background

Mefenamic acid is a nonsteroidal anti-inflammatory drug (NSAID). The mechanism of action of mefenamic acid, like that of other NSAIDs, is not completely understood but involves inhibition of cyclooxygenase (COX-1 and COX-2).

Mefenamic acid has been described to undergo metabolism by CYP2C9; it is also glucuronidated directly ([DrugBank DB00784](#5-references)).

Furthermore, mefenamic acid is known to be a potent inhibitor of uridine diphosphate-glucuronosyltransferase 1A9 (UGT1A9) and used in clinical drug-drug interaction (DDI) studies as a perpetrator to investigate the DDI potential of potential UGT1A9 substrates.

The presented model building and evaluation report evaluates the performance of a PBPK model for mefenamic acid in adults.

The objective is to establish a whole-body PBPK model for mefenamic acid featuring:

* a description of the systemic plasma concentration of mefenamic acid after oral administration.
* reversible UGT1A9 inhibition. 

The presented model building and evaluation report evaluates the performance of the PBPK model for mefenamic acid in (healthy) adults.

## 2.2 Data<a id="data"></a>

### 2.2.1 In vitro / physicochemical Data

A literature search was performed to collect available information on physicochemical properties of mefenamic acid. The obtained information from literature is summarized in the table below. 

| **Parameter**   | **Unit** | **Value** | Source                                                       | **Description**                                 |
| :-------------- | -------- | --------- | ------------------------------------------------------------ | ----------------------------------------------- |
| MW              | g/mol    | 241.29    | [DrugBank DB00784](#5-references)                           | Molecular weight                                |
| pK<sub>a</sub>  |          | 4.2       | [DrugBank DB00784](#5-references)                           | Acid dissociation constant                      |
| Solubility (pH) | mg/L     | 20 (7)    | [DrugBank DB00784](#5-references)                           | Aqueous Solubility                              |
| fu              | %        | 1.9       | [Goosen 2016](#5-references)                                 | Fraction unbound in plasma                      |

With regard to UGT1A9 inhibition, mefenamic acid inhibited propofol glucuronidation in recombinant UGT1A9 by a mixed-type mechanism, however close to a competitive type (*BAYER in-house*: [Jungmann 2019](#5-references)):

| **Parameter**    | **Unit** | **Value** | Source                          | **Description**                      |
| :--------------- | -------- | --------- | ------------------------------- | ------------------------------------ |
| K<sub>i</sub>    | µmol/L   | 0.30      | [Jungmann 2019](#5-references) | Inhibition constant                  |
| Alpha            |          | 71        | [Jungmann 2019](#5-references) | Alpha value in mixed-type inhibition |
| fu<sub>inc</sub> | %        | 1         | [Fricke 2020](#5-references)   | determined *in vitro* at 0.30 µmol/L of mefenamic acid |

### 2.2.2 Clinical Data

A literature search was performed to collect available clinical data on mefenamic acid in adults. 

The following publications were found for adults and, unless noted otherwise, used for model building and evaluation:

| Publication                                           | Study description                                            |
| :---------------------------------------------------- | :----------------------------------------------------------- |
| [Hamaguchi 1987](#5-references)                       | Treatment 2 - fasted with 200 mL of water - with an oral single dose of 250 mg, fasted |
| [Mahadik 2012](#5-references)                         | Reference (Ponstan capsule)  with an oral single dose of 250 mg, fasted |
| [Rouini 2005](#5-references)                          | Reference (Ponstan capsule) with an oral single dose of 250 mg, fasted |
| [Becker 2015](#5-references) <br />(*BAYER in-house*) | 500 mg oral dose, fed condition,<br />then 250 mg oral dose every 6 h (8 doses), fed conditions<br />***confidential data*** |
| [Goosen 2017](#5-references)                          | ***not used for model building (unclear study design)***<br />500 mg oral dose |
