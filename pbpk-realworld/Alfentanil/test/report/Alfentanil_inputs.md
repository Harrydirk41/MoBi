# Building and Evaluation of a PBPK Model for alfentanil in Adults


> You are handed this report and the clinical data files. Extract the

> physicochemical and in-vitro inputs yourself, decide the model

> structure, and fit the parameters. The reference model's chosen

> methods and fitted values are withheld.


# Background

Alfentanil is a potent analgesic synthetic opioid. It is fast but short-acting and used for anesthesia during surgery. Alfentanil is metabolized solely by CYP3A4 ([Phimmasone 2001](#5-references)). Like midazolam, alfentanil is not a substrate for P-gp ([Wandel 2002](#5-references)) and less than 1% of an alfentanil dose is excreted unchanged in urine ([Meuldermans 1988](#5-references)).

Although in clinical use alfentanil is always administered intravenously (iv), some DDI studies published plasma concentration-time profiles of alfentanil following oral ingestion. The presented alfentanil model was established using clinical PK data of 8 publications, covering iv and oral (po) administration and a dosing range from 0.015 to 0.075 mg/kg as well as absolute doses of 1 mg iv and 4 mg po. The established model is based on the model developed by Hanke *et al.* ([Hanke 2018](#5-references)) and applies metabolism by CYP3A4 and glomerular filtration.

## 2.2 Data<a id="data"></a>

### 2.2.1 In vitro / physicochemical Data

A literature search was performed to collect available information on physicochemical properties of alfentanil. The obtained information from literature is summarized in the table below. 

| **Parameter**   | **Unit** | **Value**  | Source                            | **Description**                                           |
| :-------------- | -------- | ---------- | --------------------------------- | --------------------------------------------------------- |
| MW              | g/mol    | 416.52     | [DrugBank DB00802](#5-references) | Molecular weight                                          |
| pK<sub>a</sub>  |          | 6.5 (base) | [Jansson 2008](#5-references)     | Acid dissociation constant                                |
| Solubility (pH) | mg/L     | 992 (6.5)  | [Baneyx 2014](#5-references)      | Solubility                                                |
| fu              | %        | 8.6        | [Gertz 2010](#5-references)       | Fraction unbound in plasma                                |
|                 |          | 10.0       | [Edginton 2008](#5-references)    | Fraction unbound in plasma                                |
|                 |          | 12.0       | [Almond 2016](#5-references)      | Fraction unbound in plasma                                |

### 2.2.2 Clinical Data

A literature search was performed to collect available clinical data on alfentanil in healthy adults.

#### 2.2.2.1 Model Building

The following studies were used for model building:

| Publication                      | Arm / Treatment / Information used for model building        |
| :------------------------------- | :----------------------------------------------------------- |
| [Ferrier 1985](#5-references)    | Healthy subjects with a single iv dose of 0.05 mg/kg         |
| [Kharasch 1997](#5-references)   | Healthy subjects with a single iv dose of 0.02 mg/kg         |
| [Kharasch 2004](#5-references)   | Healthy subjects with a single iv dose of 0.015 mg/kg, healthy subjects with a single oral dose of 0.06 mg/kg |
| [Kharasch 2011](#5-references)   | Healthy subjects with a single iv dose of 0.015 mg/kg, healthy subjects with a single oral dose of 0.075 mg/kg |
| [Kharasch 2011b](#5-references)  | Healthy subjects with an iv dose of 1 mg, healthy subjects with an oral dose of 1 mg. Publication compares sequential and simultaneous dosing of oral deuterated and intravenous unlabeled alfentanil. Furthermore, IV and oral administration of alfentanil is combined with grapefruit juice. Grapefruit juice is considered to have no effect on hepatic clearance, and, hence, no effect on IV administered alfentanil |
| [Kharasch 2012](#5-references)   | Healthy subjects with a single iv dose of 0.02 mg/kg, healthy subjects with a single oral dose of 0.043 mg/kg |
| [Meistelman 1987](#5-references) | Healthy subjects with a single iv dose of 0.02 mg/kg         |
| [Phimmasone 2001](#5-references) | Healthy subjects with a single iv dose of 0.015 mg/kg        |
