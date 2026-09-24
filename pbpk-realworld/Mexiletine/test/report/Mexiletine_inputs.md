# Building and evaluation of a PBPK model for Mexiletine in adults


> You are handed this report and the clinical data files. Extract the

> physicochemical and in-vitro inputs yourself, decide the model

> structure, and fit the parameters. The reference model's chosen

> methods and fitted values are withheld.


# Background

The presented PBPK model of mexiletine has been developed to be used in a PBPK Drug-Drug-Interactions (DDI) network with mexiletine as a substrate and inhibitor of CYP1A2.

Mexiletine is a non-selective voltage-gated sodium channel blocker which belongs to the Class IB anti-arrhythmic group of medicines. It is used to treat arrhythmias within the heart, or seriously irregular heartbeats. The following ADME properties characterize mexiletine pharmacokinetics ([Mexiletine Drugs.com](#5-references), [SmPC Namuscla](#5-references)):

**Absorption**: Mexiletine is well absorbed (~90%) from the gastrointestinal tract. Its first-pass metabolism is low. Peak blood levels are reached in two to three hours.

**Distribution**: It is 50% to 60% bound to plasma protein, with a volume of distribution of 5 to 7 l/kg.

**Metabolism**: Mexiletine is mainly metabolized in the liver, the primary pathway being CYP2D6 metabolism, although it is also a substrate for CYP1A2. With involvement of CYP2D6, there can be either poor or extensive metabolizer phenotypes. The metabolic degradation proceeds via various pathways including aromatic and aliphatic hydroxylation, dealkylation, deamination and N-oxidation. Several of the resulting metabolites are submitted to further conjugation with glucuronic acid (phase II metabolism); among these are the major metabolites p-Hydroxymexiletine, hydroxy-methylMexiletine, and N-hydroxy-Mexiletine.

**Elimination**: In normal subjects, the plasma elimination half-life of mexiletine is approximately 10 to 12 hours. Approximately 10% is excreted unchanged by the kidney.

After i.v. administration, mexiletine shows linear pharmacokinetics in the dose range 167-200 mg (free base) and healthy volunteers and patients show similar profiles. p.o. data appear dose linear in the range of 83-500 mg as free base. The summary of product characteristics (SPC) for mexiletine ([Mexiletine, Drugs.com](#5-references)) reports that absorption rate of mexiletine is reduced in clinical situations such as acute myocardial infarction in which gastric emptying time is increased. For this reason, clinical data from patients after p.o. administration have not been considered during model development.

## 2.2 Data<a id="data"></a>

### 2.2.1 In vitro and physico-chemical data

A literature search was performed to collect available information on physico-chemical properties of mexiletine, see [Table 1](#table-1).

| **Parameter**                   | **Unit** | **Value**             | Source                            | **Description**                                              |
| :------------------------------ | -------- | --------------------- | --------------------------------- | ------------------------------------------------------------ |
| MW<sup>+</sup>                  | g/mol    | 179.26                | [DrugBank DB00379](#5-references) | Molecular weight                                             |
| pK<sub>a,base</sub><sup>+</sup> |          | 9.2                   | [DrugBank DB00379](#5-references) | Basic dissociation constant                                  |
| Solubility (pH)<sup>+</sup>     | mg/mL    | 0.54<br />(7)         | [DrugBank DB00379](#5-references) | Aqueous Solubility                                           |
| logD                            |          | 2.15 - 2.46           | [DrugBank DB00379](#5-references) | Distribution coefficient                                     |
| fu<sup>+</sup>                  | %        | 50                    | [DrugBank DB00379](#5-references) | Fraction unbound in plasma                                   |
| CYP1A2 CL                       | l/h      | 0.5 - 11              | [Labbé 2000](#5-references)       | Partial metabolic clearance of mexiletine to N-Hydroxymexiletine |
| CYP2D6 CL                       | l/h      | 12 - 13               | [Labbé 2000](#5-references)       | Difference in non-renal CL between CYP2D6 extensive and poor metabolizers |
| Unspecific liver CL             | l/h      | 12 - 24               | [Labbé 2000](#5-references)       | Non-renal CL – CYP2D6 CL – CYP1A2 CL                         |
| Renal elimination<sup>+</sup>   | l/h      | 1.8<sup>+</sup> - 2.1 | [Labbé 2000](#5-references)       | Renal clearance                                              |
| Ki CYP1A2<sup>+</sup>           | µmol/l   | 0.28                  | [Wei 1991](#5-references)         | CYP1A2 inhibition constant                                   |

**Table 1:**<a name="table-1"></a> Physico-chemical and *in-vitro* metabolization properties of mexiletine extracted from literature. *<sup>+</sup>: Value used in final model*

### 2.2.2 Clinical data

A literature search was performed to collect available clinical data on mexiletine, see [Table 2](#table-2).

| **Source**           | Route | **Dose [mg]/**  **Schedule \*** | **Pop.**     | **Sex** | **N** | **Form.** | **Comment**                       |
| -------------------- | ------------------------------- | ------------ | ------- | --------------------------------- | --------------------------------- | --------------------------------- | -------------------- |
| [Campbell 1978](#5-references)<sup>+</sup> | i.v.   | 200                                      | HV                           | m       | 5     | solution  |              |
| [Campbell 1978](#5-references)<sup>+</sup>    | p.o.  | 200                                      | HV                           | m       | 5     | -         |                              |
| [Begg 1982](#5-references)<sup>+</sup>     | p.o.  | 333.24                                   | HV                           | m/f     | 6     | tablet    | 6 IDs                        |
| [Labbé 2000](#5-references)                   | p.o.  | 83.31 b.i.d.                             | HV                           | m/f     | 1     | -         | EM/PM                        |
| [Campbell 1978](#5-references)<sup>+</sup>    | i.v.  | 200                                      | patients                     | -       | 10    | solution  |                              |
| [Pringle 1986](#5-references)<sup>+</sup>     | p.o.  | 83.31 - 166.62 - 249.9 -  333.24 - 499.9 | HV                           | m       | 12    | capsule   |                              |
| [Kusumoto 1998](#5-references)<sup>+</sup>    | p.o.  | 166.62                                   | HV                           | m       | 9     | capsule   |                              |
| [Pentikäinen 1984](#5-references)<sup>+</sup> | i.v.  | 166.62                                   | acute myocardial  infarction | -       | 18    | solution  | acute myocardial  infarction |
| [Joeres 1987](#5-references)<sup>+</sup>      | p.o.  | 200                                      | HV                           | -       | 1     | -         |                              |

**Table 2:**<a name="table-2"></a> Literature sources of clinical concentration data of mexiletine used for model development and validation. *\*: single dose unless otherwise specified; EM: extensive metabolizers; PM: poor metabolizers; <sup>+</sup>: Data used for final parameter identification*
