# Building and evaluation of a PBPK model for sufentanil in adults


> You are handed this report and the clinical data files. Extract the

> physicochemical and in-vitro inputs yourself, decide the model

> structure, and fit the parameters. The reference model's chosen

> methods and fitted values are withheld.


# Background

The presented model building and evaluation report evaluates the performance of a PBPK model for sufentanil in adults.

Sufentanil is a potent synthetic opioid. It has a high lipid solubility, which accounts for the fast onset when given intravenously. The commercial solution comes as preservative‐free sufentanil citrate, injectable with a pH of 4.5–7.0 (Jansen‐Cilag AB, Sweden). Sufentanil is solely metabolised by CYP3A4. Due to the high hepatic extraction ratio, for the overall clearance of sufentanil both CYP3A4 activity as well as liver blood flow rate play dominant roles in the elimination in adult populations. The final sufentanil model applies metabolism by CYP3A4 and glomerular filtration and adequately described the pharmacokinetics of sufentanil in adults.

The sufentanil model is a whole-body PBPK model, allowing for dynamic translation between individuals. The sufentanil report demonstrates the level of confidence in the sufentanil PBPK model build with the OSP suite with regard to reliable predictions of sufentanil PK in adults during model-informed drug development.

## 2.2 Data used<a id="data"></a>

### 2.2.1 In vitro / physicochemical data

A literature search was performed to collect available information on physicochemical properties of sufentanil. The obtained information from literature is summarized in the table below, and is used for model building.

| **Parameter**   | **Unit**    | **Literature value (reference)** | **Description**                                  |
| :-------------- | ----------- | ----------------------------------- | ------------------------------------------------ |
| MW              | g/mol       | 386.6 ([Zhou 2017](#5-references))  | Molecular weight                                 |
| pKa             |             | 8 ([Zhou 2017](#5-references)) | Base dissociation constant                       |
| Solubility (pH) | mg/L        | 0.076 (7) ([Roy 1988](#5-references)) | Solubility                                       |
| fu              |             | 0.075 ([Zhou 2017](#5-references))  | Fraction unbound                                 |
| CLr*            | L/h         | 1 ([Zhou 2017](#5-references))      | Renal clearance                                  |
| CYP3A4 CLint* | µl/min/pmol | 20.74 ([Zhou 2017](#5-references))  | Cytochrome-P450 3A4 mediated intrinsic clearance |

*CLr and CYP3A4int parameters are built in PK-Sim as glomerular filtration (GF) and CYP3A4 - first order intrinsic clearance, respectively.

### 2.2.2 Clinical data

A literature search was performed to collect available clinical data on sufentanil in adults. 

The following publications were found in adults for model building and evaluation:

| Publication                       | Study description                                            |
| :-------------------------------- | :----------------------------------------------------------- |
| [Bovill 1984](#5-references) | The pharmacokinetics of sufentanil in elective surgical patients, without hepatic or renal dysfunction. |
| [Willsie 2015](#5-references) | Pharmacokinetic properties of single- and repeated-dose sufentanil sublingual tablets in healthy volunteers. |
| [Taverne 1992](#5-references) | Comparative absorption and distribution pharmacokinetics of intravenous and epidural sufentanil in elective surgical patients, without hepatic or renal dysfunction. |
