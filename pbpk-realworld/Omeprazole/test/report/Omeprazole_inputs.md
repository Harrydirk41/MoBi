# Building and evaluation of a PBPK model for Omeprazole in adults


> You are handed this report and the clinical data files. Extract the

> physicochemical and in-vitro inputs yourself, decide the model

> structure, and fit the parameters. The reference model's chosen

> methods and fitted values are withheld.


# Background

The presented PBPK model of omeprazole has been developed to be used in a PBPK Drug-Drug-Interactions (DDI) network with omeprazole as a substrate of CYP2C19 and CYP3A4 and an inhibitor of CYP2C19.

Omeprazole is a proton pump inhibitor (PPI) for the treatment of gastric acid related diseases. Omeprazole is administered as a racemic mixture of its two enantiomers, S-omeprazole and R-omeprazole.

The following ADME characteristics for (R-/S-) omeprazole were taken from omeprazole 20 mg capsule [FDA SPC](#5-references):

**Absorption**: Absorption of (R-/S-) omeprazole is rapid, with peak plasma levels occurring approximately 1-2 hours after dose. Absorption of omeprazole takes place in the small intestine and is usually completed within 3-6 hours. The systemic availability (bioavailability) from a single oral dose of omeprazole is approximately 40%. After repeated once-daily administration, the bioavailability increases to about 60%.

**Distribution**: The apparent volume of distribution in healthy subjects is approximately 0.3 l/kg body weight. Omeprazole is 97% plasma protein bound.

**Metabolism**: (R-/S-) Omeprazole is completely metabolized by the cytochrome P450 system (CYP). The major part of its metabolism is dependent on the polymorphically expressed CYP2C19, responsible for the formation of hydroxyomeprazole, the major metabolite in plasma. The remaining part is dependent on another specific isoform, CYP3A4, responsible for the formation of omeprazole sulphone. As a consequence of high affinity of omeprazole to CYP2C19, there is a potential for competitive inhibition and metabolic drug-drug interactions. No metabolite has been found to be pharmacologically active.

**Elimination**: Almost 80% of an oral dose of omeprazole is excreted as metabolites in the urine, the remainder in the faeces, primarily originating from bile secretion.

**Linearity/non-linearity:** The AUC of omeprazole increases with repeated administration following a non-linear dose-AUC relationship. This time- and dose- dependency is due to a decrease of first pass metabolism and systemic clearance probably caused by an inhibition of the CYP2C19 enzyme by omeprazole and/or its metabolites (e.g. the sulphone).

**Stereoselectivity:** Because both enantiomers are substrates and inhibitors of CYP2C19, the PK profile of each enantiomer is expected to differ when the racemic mixture is administered as compared to as a single agent because of mutual inhibition between the enantiomers ([Äbelö 2000](#5-references)). The Ki for R-omeprazole is 2-5 higher compared to S-omeprazole ([Liu 2005](#5-references), [Wu 2014](#5-references)).

## 2.2 Data<a id="data"></a>

### 2.2.1 In vitro and physico-chemical data

A literature search was performed to collect available information on physico-chemical properties of S- and R-omeprazole and summarized in [Table 1](#table-1) and [Table 2](#table-2), respectively.

| **Parameter**                                              | **Unit** | **Value**     | Source                                          | **Description**                                              |
| :--------------------------------------------------------- | -------- | ------------- | ----------------------------------------------- | ------------------------------------------------------------ |
| MW<sup>+</sup>                                             | g/mol    | 345.42        | [DrugBank DB00338](#5-references)               | Molecular weight.                                            |
| pK<sub>a,acid</sub><sup>+</sup>                            |          | 9.29          | [DrugBank DB00338](#5-references)               | Acidic dissociation constant                                 |
| pK<sub>a,base</sub><sup>+</sup>                            |          | 4.77          | [DrugBank DB00338](#5-references)               | Basic dissociation constant                                  |
| Solubility (pH)<sup>+</sup>                                | mg/mL    | 0.36<br />(7) | [DrugBank DB00338](#5-references)               | Aqueous Solubility                                           |
| logD                                                       |          | 2.23          | [Ogilvie 2011](#5-references)                   | Distribution coefficient                                     |
| fu<sup>+</sup>                                             | %        | 3             | [Nexium prescribing information](#5-references) | Fraction unbound in plasma                                   |
| Kinact  CYP2C19<sup>+</sup>                                | l/h      | 5             | [Wu 2014](#5-references)                        | Kinact of time dependent inhibition on CYP2C19               |
| K<sub>i</sub> CYP2C19<sup>+</sup>                          | µM       | 0.3           | [Wu 2014](#5-references)                        | KI of time dependent inhibition on CYP2C19                   |
| K<sub>i</sub> CYP2C19 (competitive inhibition)<sup>+</sup> | µM       | 3.1           | [Liu 2005](#5-references)                       | The total ki  value reported by Liu was 3.4 µmol/L and corrected with an fu_mic of 0.92 |
| Renal Elimination<sup>+</sup>                              | l/h      | 0.037         | [Wu 2014](#5-references)                        | Assumed same as omeprazole                                   |

**Table 1:**<a name="table-1"></a> Physico-chemical and *in-vitro* metabolization properties of S-omeprazole extracted from literature. *<sup>+</sup>: Value used in final model*

| **Parameter**                                              | **Unit** | **Value**     | Source                            | **Description**                                              |
| :--------------------------------------------------------- | -------- | ------------- | --------------------------------- | ------------------------------------------------------------ |
| MW<sup>+</sup>                                             | g/mol    | 345.42        | [DrugBank DB00338](#5-references) | Molecular weight.                                            |
| pK<sub>a,acid</sub><sup>+</sup>                            |          | 9.29          | [DrugBank DB00338](#5-references) | Acidic dissociation constant                                 |
| pK<sub>a,base</sub><sup>+</sup>                            |          | 4.77          | [DrugBank DB00338](#5-references) | Basic dissociation constant                                  |
| Solubility (pH)<sup>+</sup>                                | mg/mL    | 0.36<br />(7) | [DrugBank DB00338](#5-references) | Aqueous Solubility                                           |
| logD                                                       |          | 2.23          | [Ogilvie 2011](#5-references)     | Distribution coefficient                                     |
| fu<sup>+</sup>                                             | %        | 4             | [Ogilvie 2011](#5-references)     | Fraction unbound in plasma; assumed same as omeprazole       |
| Kinact  CYP2C19<sup>+</sup>                                | l/h      | 4             | [Wu 2014](#5-references)          | Kinact of time dependent inhibition on CYP2C19               |
| K<sub>i</sub> CYP2C19<sup>+</sup>                          | µM       | 1.6           | [Wu 2014](#5-references)          | KI of time dependent inhibition on CYP2C19                   |
| K<sub>i</sub> CYP2C19 (competitive inhibition)<sup>+</sup> | µM       | 5.3           | [Liu 2005](#5-references)         | The total ki  value reported by Liu was 5.7 µmol/L and corrected with an fu_mic of 0.92 |
| Renal Elimination<sup>+</sup>                              | l/h      | 0.037         | [Wu 2014](#5-references)          | Assumed same as omeprazole                                   |

**Table 2:**<a name="table-2"></a> Physico-chemical and *in-vitro* metabolization properties of R-omeprazole extracted from literature. *<sup>+</sup>: Value used in final model*

### 2.2.2 Clinical data

A literature search was performed to collect available clinical data on omeprazole in adults. Data used for model development and validation for omeprazole, S-omeprazole, and R-omeprazole are listed in [Table 3](#table-3), [Table 4](#table-4), and [Table 5](#table-5), respectively.

| **Source**           | **Route** | **Dose [mg]/**  **Schedule \*** | **Pop.**     | **Sex** | **N** | **Form.** | **Comment**                       |
| -------------------- | --------- | ------------------------------- | ------------ | ------- | ----- | --------- | --------------------------------- |
| [Andersson 1990](#5-references)<sup>+</sup> | i.v.      | 40 - 80                         | HV                        | M       | 10    | solution                      |                    |
| [Andersson 1990](#5-references)              | p.o.      | 40 - 80                         | HV                        | M       | 10    | Oral solution                 |                    |
| [Andersson 1991](#5-references)              | p.o.      | 10 - 20 – 40 q.d.               | HV                        | M       | 12    | e.c. granules                 |                    |
| [Andersson 1991](#5-references)<sup>+</sup>  | i.v.      | 10 - 20 - 40                    | HV                        | M       | 12    | solution                      |                    |
| [Andersson 1998](#5-references)              | p.o.      | 20 q.d.                         | HV                        | M       | 12    | capsule                       | EM                 |
| [Andersson 1998](#5-references)              | p.o.      | 20 q.d.                         | HV                        | M       | 2     | capsule                       | PM                 |
| [Oosterhuis 1992](#5-references)<sup>+</sup> | i.v.      | 40 - 80                         | HV                        | M       | 8     | solution                      |                    |
| [Uno 2007](#5-references)<sup>+</sup>        | i.v.      | 20                              | HV japanese               | M - F   | 6     | solution                      | hmEM               |
| [Uno 2007](#5-references)<sup>+</sup>        | i.v.      | 20                              | HV japanese               | M - F   | 6     | solution                      | PM                 |
| [Uno 2007](#5-references)                    | p.o.      | 40                              | HV japanese               | M - F   | 6     | tablet                        | hmEM               |
| [Uno 2007](#5-references)                    | p.o.      | 40                              | HV japanese               | M - F   | 6     | tablet                        | PM                 |
| [Regårdh 1990](#5-references)<sup>+</sup>    | i.v.      | 10                              | HV                        | M       | 8     | solution                      |                    |
| [Regårdh 1990](#5-references)                | p.o.      | 20                              | HV                        | M       | 8     | Oral solution                 |                    |
| [Andersson 2000](#5-references)              | p.o.      | 15 q.d.                         | HV                        | -       | 4     | Oral solution                 | EM                 |
| [Andersson 2000](#5-references)              | p.o.      | 60 q.d.                         | HV                        | -       | 5     | Oral solution                 | PM                 |
| [Hassan-Alin 2005](#5-references)            | p.o.      | 20 – 40 q.d.                    | HV                        | -       | -     | Oral solution                 |                    |
| [Cho 2002](#5-references)                    | p.o.      | 20                              | HV asian                  | -       | -     | capsule                       | EM +/- moclobemide |
| [Cho 2002](#5-references)                    | p.o.      | 20                              | HV asian                  | -       | -     | capsule                       | PM +/- moclobemide |
| [Yasui-Furukori 2004](#5-references)         | p.o.      | 40                              | HV japanese               | M - F   | 6     | omepral                       |                    |
| [Yasui-Furukori 2004](#5-references)         | p.o.      | 40                              | HV japanese               | M - F   | 6     | omepral                       | PM +/- fluvoxamine |
| [Wu 2016](#5-references)      | p.o.      | 40 q.d.                         | HV caucasian              | M - F   | 15    | gastro-resistant hard capsule |                    |

**Table 3:**<a name="table-3"></a> Literature sources of clinical concentration data of omeprazole used for model development and validation. *e.c.: enteric coated; -: respective information was not provided in the literature source; \*:single dose unless otherwise specified; EM: extensive metabolizers; PM: poor metabolizers; <sup>+</sup>: Data used for final parameter identification*

| **Source**                                     | **Route** | **Dose [mg]/**  **Schedule \***                              | **Pop.** | **Sex** | **N** | **Form.**     | **Comment** |
| ---------------------------------------------- | --------- | ------------------------------------------------------------ | -------- | ------- | ----- | ------------- | ----------- |
| [Wilder-Smith 2005](#5-references)<sup>+</sup> | i.v.      | 40 q.d.                                                      | HV       | M - F   | 39    | solution      |             |
| [Wilder-Smith 2005](#5-references)<sup>+</sup> | p.o.      | 40 q.d.                                                      | HV       | M - F   | 37    | capsule       |             |
| [Hassan-Alin 2000](#5-references)<sup>+</sup>  | i.v.      | 20 - 40                                                      | HV       | -       | -     | solution      |             |
| [Andersson 2000](#5-references)                | p.o.      | 15 q.d.                                                      | HV       | -       | 4     | oral solution | EM          |
| [Andersson 2000](#5-references)<sup>+</sup>    | p.o.      | 60 q.d.                                                      | HV       | -       | 5     | oral solution | PM          |
| [Hassan-Alin 2005](#5-references)<sup>+</sup>  | p.o.      | 20 - 40 q.d.                                                 | HV       | -       | -     | oral solution |             |
| [FDA Nexium Review](#5-references)             | p.o.      | 40                                                           | HV       | M       | -     | -             | EM          |
| [FDA Nexium Review](#5-references)<sup>+</sup> | p.o.      | 40                                                           | HV       | M       | -     | -             | PM          |
| [Rohss 2007](#5-references)<sup>+</sup>        | i.v.      | 120mg(30min)+8mg/h - 120mg(2h)+8mg/h -  80mg(30min)+4mg/h - 80mg(30min)+8mg/h - 40mg(30min)+8mg/h | HV       | -       | 25    | solution      |             |

**Table 4:**<a name="table-4"></a> Literature sources of clinical concentration data of S-omeprazole used for model development and validation. *e.c.: enteric coated; -: respective information was not provided in the literature source; \*:single dose unless otherwise specified; EM: extensive metabolizers; PM: poor metabolizers; <sup>+</sup>: Data used for final parameter identification*

| **Source**                                    | **Route** | **Dose [mg]/**  **Schedule \*** | **Pop.** | **Sex** | **N** | **Form.**     | **Comment** |
| --------------------------------------------- | --------- | ------------------------------- | -------- | ------- | ----- | ------------- | ----------- |
| [Andersson 2000](#5-references)<sup>+</sup>   | p.o.      | 15 q.d.                         | HV       | -       | 4     | oral solution | EM          |
| [Andersson 2000](#5-references)<sup>+</sup>   | p.o.      | 60 q.d.                         | HV       | -       | 5     | oral solution | PM          |
| [Hassan-Alin 2005](#5-references)<sup>+</sup> | p.o.      | 20 - 40 q.d.                    | HV       | -       | -     | oral solution |             |

**Table 5:**<a name="table-5"></a> Literature sources of clinical concentration data of R-omeprazole used for model development and validation. *e.c.: enteric coated; -: respective information was not provided in the literature source; \*:single dose unless otherwise specified; EM: extensive metabolizers; PM: poor metabolizers; <sup>+</sup>: Data used for final parameter identification*
