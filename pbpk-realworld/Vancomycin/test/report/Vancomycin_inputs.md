# Building and evaluation of a PBPK model for vancomycin in adults


> You are handed this report and the clinical data files. Extract the

> physicochemical and in-vitro inputs yourself, decide the model

> structure, and fit the parameters. The reference model's chosen

> methods and fitted values are withheld.


# Background

The presented model building and evaluation report evaluates the performance of a PBPK model for vancomycin in adults.

Vancomycin is a glycopeptide antibiotic related to ristocetin that inhibits bacterial cell wall assembly and is used to treat a number of bacterial infections. It can be administered intravenously, as well as orally in case of diarrhea therapy. Vancomycin is mainly eliminated via glomerular filtration (GF). A previous PBPK model for vancomycin using PK-Sim was reported by Radke et al. ([Radke 2017](#5-references)), with the dose fraction excreted unchanged into urine in adults being 90% with 10% hepatic elimination. Our final vancomycin model was rebuilt that applies only GFR mediated clearance that adequately described the pharmacokinetics in adults. No further improvement of vancomycin pharmacokinetics could be determined after introducing hepatic clearance.

The vancomycin model is a whole-body PBPK model, allowing for dynamic translation between individuals. The vancomycin report demonstrates the level of confidence in the vancomycin  PBPK model built with the OSP suite with regard to reliable predictions of vancomycin PK in adults during model-informed drug development.

## 2.2 Data used<a id="data"></a>

### 2.2.1 In vitro / physicochemical data

A literature search was performed to collect available information on physicochemical properties of vancomycin. The obtained information from literature is summarized in the table below, and is used for model building.

| **Parameter**      | **Unit**  | **Literature value (reference)**                             | **Description**                                 |
| :----------------- | --------- | ------------------------------------------------------------ | ----------------------------------------------- |
| MW                 | g/mol     | 1449.3 ([Radke 2017](#5-references))                         | Molecular weight                                |
| pKa                |           | Acid 2.18, Base 7.75, Base 8.89 ([Radke 2017](#5-references)) | Acid/base dissociation constant                 |
| Solubility (pH)    | mg/L      | 225 (7) ([Drugbank.ca](#5-references))                       | Solubility                                      |
| fu                 |           | 0.48 ([Zhou 2016](#5-references)), 0.67 ([Radke 2017](#5-references)) | Fraction unbound                                |
| GFR fraction       | µM        | 1 ([Zhou 2016](#5-references))                               | fraction of Glomerular filtration rate          |
| Hepatic clearance* | mL/min/kg | 0.11 ([Radke 2017](#5-references))                           | Hepatic clearance                               |
| Renal clearance*   | mL/min/kg | 0.95 ([Radke 2017](#5-references))                           | Renal clearance                                 |

*Both Hepatic and Renal clearance reported by others have not been used in the final model.

### 2.2.2 Clinical data

A literature search was performed to collect available clinical data on vancomycin in adults. 

The following publications were found in adults for model building and evaluation:

| Publication                       | Study description                                            |
| :-------------------------------- | :----------------------------------------------------------- |
| [Boeckh 1988](#5-references) | Pharmacokinetics and serum bactericidal activity of vancomycin alone and in combination with ceftazidime in healthy volunteers |
| [Healy 1987](#5-references) | Comparison of steady-state pharmacokinetics of two dosage regimens of vancomycin in normal volunteers |
