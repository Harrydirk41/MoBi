# Building and evaluation of a PBPK model for amikacin in adults


> You are handed this report and the clinical data files. Extract the

> physicochemical and in-vitro inputs yourself, decide the model

> structure, and fit the parameters. The reference model's chosen

> methods and fitted values are withheld.


# Background

The presented model building and evaluation report evaluates the performance of a PBPK model for amikacin in adults.

Amikacin is a semi-synthetic aminoglycoside antibiotic used for a number of bacterial infections. Amikacin is administered in several forms, including intravenous or intramuscular injection. The PBPK model for amikacin was previously developed in PK-sim for adults ([Wendl 2011](#5-references)) and preterm neonates ([Claassen 2015](#5-references)). As the latter model was built more recently, this PBPK model was used to evaluate the predictive performance of glomerular filtration rate (GFR) mediated clearance in adults without further changes. In this chapter we show that amikacin adequately described the pharmacokinetics of amikacin in adults, based on the PBPK model build and reported in preterm neonates.

The amikacin model is a whole-body PBPK model, allowing for dynamic translation between individuals with GFR based renal elimination. The amikacin report demonstrates the level of confidence in the amikacin PBPK model build with the OSP suite with regard to reliable predictions of amikacin PK in adults during model-informed drug development.

## 2.2 Data<a id="data"></a>

### 2.2.1 In vitro / physicochemical data

A literature search was performed to collect available information on physicochemical properties of amikacin. The obtained information from literature is summarized in the table below, and is used for model building.

| **Parameter**   | **Unit** | **Value (reference)**                            | **Description**                        |
| :-------------- | -------- | ------------------------------------------------ | -------------------------------------- |
| MW              | g/mol    | 588.6 ([Claassen 2015](#5-references))           | Molecular weight                       |
| pKa             |          | 9.7, 8.92, 8.13 ([Claassen 2015](#5-references)) | Acid dissociation constants            |
| Solubility (pH) | mg/L     | 50 (7)  ([Drugbank.ca](#5-references))           | Solubility                             |
| logMA           |          | -0.48 ([Claassen 2015](#5-references))           | The logarithm of the membrane affinity |
| fu              |          | 1 ([Claassen 2015](#5-references))               | Fraction unbound                       |
| GFR fraction    |          | 1 ([Claassen 2015](#5-references))               | Glomerular Filtration Rate fraction    |

### 2.2.2 Clinical data

A literature search was performed to collect available clinical data on amikacin in adults. 

The following publication was found in adults for the evaluation of the reported amikacin PBPK model:

| Publication                       | Study description                                            |
| :-------------------------------- | :----------------------------------------------------------- |
| [Walker 1979](#5-references) | The pharmacokinetics of amikacin and gentamicin in healthy volunteers |
