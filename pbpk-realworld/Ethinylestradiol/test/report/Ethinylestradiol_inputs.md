# Building and evaluation of a PBPK model for Ethinylestradiol in adults


> You are handed this report and the clinical data files. Extract the

> physicochemical and in-vitro inputs yourself, decide the model

> structure, and fit the parameters. The reference model's chosen

> methods and fitted values are withheld.


# Background

The presented PBPK model of ethinylestradiol (EE) has been developed to be used in a PBPK Drug-Drug-Interactions (DDI) network with ethinylestradiol as perpetrator of CYP1A2.

Ethinylestradiol is an estrogen medication which is used widely as a birth control pills in combination with progestins. The following ADME properties characterize ethinylestradiol ([SmPC Namuscla](#5-references), [FDA. QUARTETTE](#5-references)):

**Absorption**: ethinylestradiol is rapidly and completely absorbed from the gut but it undergoes some first pass metabolism in the gut wall (mediated by a.o. CYP3A4 ([Wiesinger 2015](#5-references), [Wang 2004](#5-references))). After oral administration, an initial peak occurs in plasma at 2 to 3 hours, with a secondary peak at about 12 hours after dosing; the second peak is interpreted as evidence for extensive enterohepatic circulation of ethinylestradiol.

**Distribution**: ethinylestradiol is rapidly distributed throughout most body tissues with the largest concentration found in adipose tissue. It distributes into breast milk, with low concentrations. More than 80% of ethinylestradiol in serum is conjugated as sulphate and almost all the conjugated form is bound to albumin.

**Metabolism**: ethinylestradiol is metabolized in the liver. Hydroxylation appears to be the main metabolic pathway. 60% of a dose is excreted in the urine and 40% in the faeces. 

**Excretion**: About 30% is excreted in the urine and bile as the glucuronide or sulphate conjugate. The rate of metabolism of ethinylestradiol is affected by several factors, including enzyme-inducing agents, antibiotics, and cigarette smoking. The elimination half-life of ethinylestradiol ranges from 5 to 16 hours.

After i.v. administration, ethinylestradiol displays approximately linear dose relationship in the dose range 30-100 µg. A wide variability is present in the terminal part of the dose-normalized concentrations.

After p.o. single dose, ethinylestradiol shows linear dose relationship in the dose range 30-3000 µg. Secondary peaks can be observed in individual data, compatible with enterohepatic re-circulation. However, mean data do not display such feature as a result of such peak being averaged out. Therefore, enterohepatic re-circulation was not taken into account in the model.

## 2.2 Data<a id="data"></a>

### 2.2.1 In vitro and physico-chemical data

A literature search was performed to collect available information on physico-chemical properties of ethinylestradiol, see [Table 1](#table-1).

| **Parameter**                   | **Unit**          | **Value**        | Source                            | **Description**                                |
| :------------------------------ | ----------------- | ---------------- | --------------------------------- | ---------------------------------------------- |
| MW<sup>+</sup>                  | g/mol             | 296.4            | [DrugBank DB00977](#5-references) | Molecular weight                               |
| pK<sub>a,acid</sub><sup>+</sup> |                   | 10.33            | [DrugBank DB00977](#5-references) | Acidic dissociation constant                   |
| Solubility (pH)<sup>+</sup>     | mg/mL             | 6.77e-3<br />(7) | [DrugBank DB00977](#5-references) | Aqueous Solubility                             |
| logD                            |                   | 3.63 - 3.9       | [DrugBank DB00977](#5-references) | Distribution coefficient                       |
| fu<sup>+</sup>                  | %                 | 3                | [DrugBank DB00977](#5-references) | Fraction unbound in plasma                     |
| CYP1A2 CL<sup>+</sup>           | µl/min/pmol       | 0.51             | [Ezuruike 2018](#5-references)    | Clearance by CYP1A2                            |
| CYP2C8 CL<sup>+</sup>           | µl/min/pmol       | 0.13             | [Ezuruike 2018](#5-references)    | Clearance by CYP2C8                            |
| CYP2C9 CL<sup>+</sup>           | µl/min/pmol       | 0.51             | [Ezuruike 2018](#5-references)    | Clearance by CYP2C9                            |
| CYP3A4 CL<sup>+</sup>           | µl/min/pmol       | 0.5              | [Ezuruike 2018](#5-references)    | Clearance by CYP3A4                            |
| Km UGT1A1<sup>+</sup>           | µmol/l            | 19.22            | [Ezuruike 2018](#5-references)    | UGT1A1 saturation constant                     |
| Vmax UGT1A1<sup>+</sup>         | pmol/min/mg prot. | 408.5            | [Ezuruike 2018](#5-references)    | Maximal metabolization rate by UGT1A1          |
| Renal Elimination<sup>+</sup>   | l/h               | 2.079            | [Stanczyk 2013](#5-references)    | Renal clearance                                |
| Clint HLM<sup>+</sup>           | µL/min/mg prot.   | 118.83           | [Ezuruike 2018](#5-references)    | Intrinsic clearance in Human Liver Microsomes |
| Ki CYP1A2                       | µmol/l            | 10.6             | [Karjalainen 2008](#5-references) | CYP1A2 inhibition constant                     |

**Table 1:**<a name="table-1"></a> Physico-chemical and *in-vitro* metabolization properties of ethinylestradiol extracted from literature. *<sup>+</sup>: Value used in final model*

### 2.2.2 Clinical data

A literature search was performed to collect available clinical data on ethinylestradiol, see [Table 2](#table-2).

| **Source**           | Route | **Dose [mg]/**  **Schedule \*** | **Pop.**     | **Sex** | **N** | **Form.** |
| -------------------- | ------------------------------- | ------------ | ------- | --------------------------------- | --------------------------------- | --------------------------------- |
| [Back 1981](#5-references)<sup>+</sup> | i.v.   | 0.03                            | HV       | F       | 5     | solution            |
| [Back 1981](#5-references)<sup>+</sup>          | p.o.  | 0.03                            | HV       | F       | 5     | tablet              |
| [Back 1979](#5-references)<sup>+</sup>          | i.v.  | 0.05                            | HV       | F       | 6     | solution            |
| [Back 1979](#5-references)<sup>+</sup>          | p.o.  | 0.05                            | HV       | F       | 6     | NA                  |
| [Back 1987](#5-references)                      | p.o.  | 0.05 q.d.                      | HV       | F       | 5     | tablet              |
| [Orme 1991](#5-references)<sup>+</sup>          | i.v.  | 0.03                            | HV       | F       | 10    | solution            |
| [Orme 1991](#5-references)<sup>+</sup>          | p.o.  | 0.03                            | HV       | F       | 10    | tablet              |
| [Kuhnz 1996](#5-references)                     | i.v.  | 0.06                            | HV       | F       | 19    | solution            |
| [Goebelsmann 1986](#5-references)<sup>+</sup>   | p.o.  | 0.03                            | HV       | F       | 24    | solution and tablet |
| [Stanczyk 1983](#5-references)<sup>+</sup>      | p.o.  | 0.12                            | HV       | F       | 24    | solution and tablet |
| [Zhang 2017](#5-references)<sup>+</sup>         | p.o.  | 0.03                            | HV       | F       | 12    | tablet              |
| [Martin 2016](#5-references)                    | p.o.  | 0.03 q.d.                       | HV       | F       | 27    | tablet              |
| [Stockis 2014](#5-references)                   | p.o.  | 0.03 q.d.                       | HV       | F       | 24    | tablet              |
| [Sidhu 2006](#5-references)                     | p.o.  | 0.03 q.d.                       | HV       | F       | 16    | tablet              |
| [Kothare 2012](#5-references)<sup>+</sup>       | p.o.  | 0.03/0.03 q.d.                  | HV       | F       | 20    | tablet              |
| [Timmer 2000](#5-references)<sup>+</sup>        | p.o.  | 0.03                            | HV       | F       | -     | tablet              |

**Table 2:**<a name="table-2"></a> Literature sources of clinical concentration data of ethinylestradiol used for model development and validation. *\*: single dose unless otherwise specified;<sup>+</sup>: Data used for final parameter identification*

## Literature in-vitro kinetics (given)

These are the measured in-vitro values from the literature - the fitted
catalytic rates are NOT given; you must determine them.

| molecule | process | parameter | value | unit |
| --- | --- | --- | --- | --- |
| UGT1A1 | MetabolizationLiverMicrosomes_MM | In vitro Vmax for liver microsomes | 408.5 | pmol/min/mg mic. protein |
| UGT1A1 | MetabolizationLiverMicrosomes_MM | Km | 19.22 | µmol/l |
