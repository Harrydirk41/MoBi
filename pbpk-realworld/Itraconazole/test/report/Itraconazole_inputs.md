# Building and evaluation of a PBPK model for Itraconazole in healthy adults


> You are handed this report and the clinical data files. Extract the

> physicochemical and in-vitro inputs yourself, decide the model

> structure, and fit the parameters. The reference model's chosen

> methods and fitted values are withheld.


# Background

Itraconazole is a triazole agent prescribed for the treatment of fungal infections. It is predominantly metabolized by CYP3A4, resulting in the sequential formation of several metabolites, starting with the major metabolite hydroxy-itraconazole, followed by keto-itraconazole and N-desalkyl-itraconazole. All three metabolites are further metabolized by CYP3A4 and parent and metabolites are reported to competitively inhibit CYP3A4 ([Isoherranen 2004](#5-references)). Therefore, the metabolites inhibit their own formation and itraconazole inhibits further conversion of its metabolites by CYP3A4. Itraconazole has been proposed as one of the most appropriate CYP3A4 inhibitors for clinical DDI studies, to replace the currently no longer recommended CYP3A4 inhibitor drug ketoconazole.

The herein presented model represents an update of the itraconazole model publisdhed by Hanke *et al.* ([Hanke 2018](#5-references)). The model was originally established using various clinical studies, covering a dosing range from 100 to 200 mg in different formulations (solution and capsules), administered under fasted conditions or together with food. Although the plasma concentrations of keto-itraconazole and N-desalkyl-itraconazole are lower than those of itraconazole and hydroxy-itraconazole, N-desalkyl-itraconazole is reported to be a very potent inhibitor *in vitro*, and integration of the two further metabolites into the model with their inhibitory effects enabled the description the strong non-linearity and plasma accumulation of itraconazole. The model applies sequential metabolism of itraconazole to hydroxy-itraconazole to keto-itraconazole to N-desalkyl-itraconazole by CYP3A4, including competitive inhibition of CYP3A4 by the parent drug and all three metabolites, and glomerular filtration. Competitive inhibition of P-gp was included for itracaonazole.

In comparison to the published version by Hanke *et al.* 2018 ([Hanke 2018](#5-references)), the dissolution and solubility has been optimized and updated for the administration of itraconazole capsules in fasted state (by integrating additional data ([Jalava 1997](#5-references)) into the optimization routine).

## 2.2 Data<a id="data"></a>

### 2.2.1 In vitro and physicochemical data

A literature search was performed to collect available information on physicochemical properties of itraconazole and metabolites. The obtained information from literature is summarized in the table below and was used for model building. Note that not all parameters were used in the final model. A list of final model parameters is provided below in later sections. 

#### Itraconazole

| **Parameter**                         | **Unit**                   | **Value**      | Source                            | **Description**                                              |
| :------------------------------------ | -------------------------- | -------------- | --------------------------------- | ------------------------------------------------------------ |
| MW                                    | g/mol                      | 705.633        | [DrugBank DB01167](#5-references) | Molecular weight                                             |
| pK<sub>a,base</sub>                   |                            | 3.7            | [Heykants 1989](#5-references)    | acid dissociation constant of conjugate acid; compound type: base |
| Solubility (pH)                       | mg/L                       | 8.0<br />(6.5) | [Taupitz 2013](#5-references)     | Solubility in FaSSIF (fasted state simulated intestinal fluid) |
| fu                                    | %                          | 0.2            | [Heykants 1989](#5-references)    | Fraction unbound in plasma                                   |
|                                       | %                          | 0.2            | [Riccardi 2015](#5-references)    | Fraction unbound in plasma                                   |
|                                       | %                          | 1.58           | [Ishigam 2001](#5-references)     | Fraction unbound in plasma                                   |
|                                       | %                          | 3.6            | [Templeton 2008](#5-references)   | Fraction unbound in plasma                                   |
| V<sub>max</sub>, K<sub>m</sub> CYP3A4 | pmol/min/nmol,<br />nmol/L | 270<br />3.9   | [Isoherranen 2004](#5-references) | CYP3A4 supersomes Michaelis-Menten kinetics                  |
| K<sub>i</sub> CYP3A4                  | nmol/L                     | 1.3            | [Isoherranen 2004](#5-references) | CYP3A4 inhibition constant                                   |
| K<sub>i</sub> P-gp                    | nmol/L                     | 8.0            | [Shityakov 2014](#5-references)   | P-gp inhibition constant                                     |

#### Hydroxy-itraconazole

| **Parameter**                         | **Unit**                   | **Value**   | Source                               | **Description**                                              |
| :------------------------------------ | -------------------------- | ----------- | ------------------------------------ | ------------------------------------------------------------ |
| MW                                    | g/mol                      | 721.633     | [DrugBank DBMET00374](#5-references) | Molecular weight                                             |
| fu                                    | %                          | 0.5         | [Templeton 2008](#5-references)      | Fraction unbound in plasma                                   |
|                                       | %                          | 1.7         | [Riccardi 2015](#5-references)       | Fraction unbound in plasma                                   |
|                                       | %                          | 2.12        | [Chen 2016](#5-references)           | Fraction unbound in plasma                                   |
| V<sub>max</sub>, K<sub>m</sub> CYP3A4 | pmol/min/nmol,<br />nmol/L | 543<br />27 | [Isoherranen 2004](#5-references)    | CYP3A4 supersomes Michaelis-Menten kinetics                  |
| K<sub>i</sub> CYP3A4                  | nmol/L                     | 14.4        | [Isoherranen 2004](#5-references)    | CYP3A4 inhibition constant                                   |

#### Keto-itraconazole

| **Parameter**                         | **Unit**                   | **Value**     | Source                                | **Description**                                              |
| :------------------------------------ | -------------------------- | ------------- | ------------------------------------- | ------------------------------------------------------------ |
| MW                                    | g/mol                      | 719.617       | [PubChem CID 53865186](#5-references) | Molecular weight                                             |
| fu                                    | %                          | 1.0           | [Riccardi 2015](#5-references)        | Fraction unbound in plasma                                   |
|                                       | %                          | 5.3           | [Templeton 2008](#5-references)       | Fraction unbound in plasma                                   |
| V<sub>max</sub>, K<sub>m</sub> CYP3A4 | pmol/min/nmol,<br />nmol/L | 86.9<br />1.4 | [Isoherranen 2004](#5-references)     | CYP3A4 supersomes Michaelis-Menten kinetics                  |
| IC<sub>50</sub> CYP3A4<sup>†</sup>    | nmol/L                     | 7.0           | [Isoherranen 2004](#5-references)     | CYP3A4 inhibition constant for half maximal inhibitory concentration at constant substrate concentration |

<sup>†</sup> The IC<sub>50</sub> values was converted to a K<sub>i</sub> value via Cheng-Prusoff equation ([Chen 1973](#5-references)) with a substrate (*midazolam*) concentration of 1 µmol/L and an assumed substrate (*midazolam*) K<sub>m</sub> value of 2.73 µmol/L: **5.12 nmol/L**

#### N-Desalkyl-itraconazole

| **Parameter**                      | **Unit** | **Value** | Source                                | **Description**                                              |
| :--------------------------------- | -------- | --------- | ------------------------------------- | ------------------------------------------------------------ |
| MW                                 | g/mol    | 649.527   | [PubChem CID 53789808](#5-references) | Molecular weight                                             |
| fu                                 | %        | 1.1       | [Riccardi 2015](#5-references)        | Fraction unbound in plasma                                   |
|                                    | %        | 1.2       | [Templeton 2008](#5-references)       | Fraction unbound in plasma                                   |
| IC<sub>50</sub> CYP3A4<sup>†</sup> | nmol/L   | 0.44      | [Isoherranen 2004](#5-references)     | CYP3A4 inhibition constant for half maximal inhibitory concentration at constant substrate concentration |

<sup>†</sup> The IC<sub>50</sub> values was converted to a K<sub>i</sub> value via Cheng-Prusoff equation ([Chen 1973](#5-references)) with a substrate (*midazolam*) concentration of 1 µmol/L and an assumed substrate (*midazolam*) K<sub>m</sub> value of 2.73 µmol/L: **0.32 nmol/L**

### 2.2.2 Clinical data

A literature search was performed to collect available clinical data on itraconazole and its metabolites in adults. The itraconazole model was built and verified using various clinical studies, covering a dosing range of 100 to 200 mg with different formulations (solution *vs.* capsule), administered under fasting conditions or together with food.

The following dosing senarios were simulated and compared to respective data:

| Route | Formulation | Food state | Dose<br />[mg] | Dosing | PK Data                                                      | Used for<sup>†</sup> | Reference                          |
| ----- | ----------- | ---------- | -------------- | ------ | ------------------------------------------------------------ | -------------------- | ---------------------------------- |
| iv    | -           | -          | 100            | SD     | Itraconazole                                                 | mv                   | [Heykants 1989](#5-references)     |
|       |             |            | 200            | OD     | Itraconazole<br />Hydroxy-Itr.                               | mbb                  | [Mouton 2006](#5-references)       |
| po    | solution    | fasted     | 100            | SD     | Itraconazole                                                 | mbb                  | [Van de Velde 1996](#5-references) |
|       |             |            |                |        | Itraconazole<br />Hydroxy-Itr.                               | mbb                  | [Van Peer 1989](#5-references)    |
|       |             |            |                | OD     | Itraconazole<br />Hydroxy-Itr.<br />Keto-Itr.<br />N-Desalkyl-Itr. | mbb                  | [Templeton 2008](#5-references)    |
|       |             |            | 200            | OD     | Itraconazole<br />Hydroxy-Itr.                               | mbb                  | [Barone 1998a](#5-references)      |
|       |             | fed        | 100            | SD     | Itraconazole<br />Hydroxy-Itr.                               | mbe                  | [Van de Velde 1996](#5-references) |
|       |             |            |                |        | Itraconazole                                                 | mbe                  | [Heykants 1989](#5-references)     |
|       |             |            | 200            | SD     | Itraconazole<br />Hydroxy-Itr.                               | mbe                  | [Barone 1998b](#5-references)      |
|       |             |            |                | OD     | Itraconazole<br />Hydroxy-Itr.                               | mbe                  | [Barone 1998a](#5-references)      |
|       | capsule     | fasted     | 100            | SD     | Itraconazole                                                 | mbe                  | [Van Peer 1989](#5-references)     |
|       |             |            |                | BID    | Itraconazole                                                 | mv                   | [Kivistö 1997](#5-references)       |
|       |             |            | 200            | SD     | Itraconazole<br />Hydroxy-Itr.                               | mbe                  | [Barone 1993](#5-references)       |
|       |             |            |                |        | Itraconazole                                                 | mv                   | [Neuvonen 1996](#5-references)      |
|       |             |            |                | OD     | Itraconazole                                                 | mbe                  | [Jalava 1997](#5-references)       |
|       |             |            |                |        | Itraconazole                                                 | mbe                  | [Olkkola 1994](#5-references)      |
|       |             |            |                |        | Itraconazole                                                 | mv                   | [Varhe 1994](#5-references)         |
|       |             | fed        | 100            | SD     | Itraconazole                                                 | mbe                  | [Van Peer 1989](#5-references)     |
|       |             |            |                | OD     | Itraconazole                                                 | mbe                  | [Van Peer 1989](#5-references)    |
|       |             |            |                |        | Itraconazole                                                 | mv                   | [Hardin 1988](#5-references)       |
|       |             |            | 200            | SD     | Itraconazole<br />Hydroxy-Itr.                               | mbe                  | [Barone 1993](#5-references)       |
|       |             |            |                |        | Itraconazole<br />Hydroxy-Itr.                               | mbe                  | [Barone 1998b](#5-references)      |
|       |             |            |                |        | Itraconazole                                                 | mv                   | [Neuvonen 1996](#5-references)      |
|       |             |            | 200            | OD     | Itraconazole                                                 | mv                   | [Hardin 1988](#5-references)       |
|       |             |            |                | BID    | Itraconazole                                                 | mbe                  | [Barone 1993](#5-references)       |
|       |             |            |                |        | Itraconazole<br />Hydroxy-Itr.                               | mv                   | [Hardin 1988](#5-references)       |

<sup>†</sup> *mbb* model building basic: used to inform the basic model parameters (see [Section 2.3.5](#235-automated-parameter-identification)); *mbe* model building extended: used to inform solubility and, if applicable, formulation-depenendent parameters only (see [Section 2.3.5](#235-automated-parameter-identification)); *mv* model verification only

## Literature in-vitro kinetics (given)

These are the measured in-vitro values from the literature - the fitted
catalytic rates are NOT given; you must determine them.

| molecule | process | parameter | value | unit |
| --- | --- | --- | --- | --- |
| CYP3A4 | rCYP450_MM | In vitro Vmax/recombinant enzyme | 0.27 | pmol/min/pmol rec. enzyme |
| CYP3A4 | rCYP450_MM | Km | 2.0688492598 | nmol/l |
