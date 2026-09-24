# Building and evaluation of a PBPK model for raltegravir in adults


> You are handed this report and the clinical data files. Extract the

> physicochemical and in-vitro inputs yourself, decide the model

> structure, and fit the parameters. The reference model's chosen

> methods and fitted values are withheld.


# Background

The presented model building and evaluation report evaluates the performance of a PBPK model for raltegravir in adults.

Raltegravir, sold under the brand name Isentress, is an antiretroviral medication used to treat HIV/AIDS by blocking the establishment of post-integration HIV latency. It is also used as part of post exposure prophylaxis to prevent HIV infection following potential exposure. Raltegravir is only taken orally and is mainly metabolized by UGT1A1 (~70%) [(Kassahun 2007](#5-references)). The final raltegravir model features metabolism by UGT1A1 and to a minor extent by UGT1A9. Additionally, there is excretion via glomerular filtration. The model adequately describes the pharmacokinetics of raltegravir in adults.

The raltegravir model is a whole-body PBPK model, allowing for dynamic translation between individuals with organs expressing UGT1A1. The raltegravir report demonstrates the level of confidence in the raltegravir PBPK model build with the OSP suite with regard to reliable predictions of raltegravir PK adults during model-informed drug development.

## 2.2 Data used<a id="data"></a>

### 2.2.1 In vitro / physico-chemical data

A literature search was performed to collect available information on physicochemical properties of raltegravir. The obtained information from literature is summarized in the table below, and is used for model building.

| **Parameter**   | **Unit**    | **Raltegravir literature**                                 | **Description**                                 |
| :-------------- | ----------- | ---------------------------------------------------------- | ----------------------------------------------- |
| MW              | g/mol       | 444.4163 ([drugbank.ca](#5-references))                    | Molecular weight                                |
| pKa             |             | 6.67 ([Moss 2012](#5-references))                          | Acid dissociation constant                      |
| Solubility (pH) | mg/L        | Reference pH-dependent table  ([Moss 2013](#5-references)) | Solubility                                      |
| fu              |             | 0.17 ([Laufer 2009](#5-references))                        | Fraction unbound                                |
| Km UGT1A1       | µM          | 99 ([Kassahun 2007](#5-references))                        | Michaelis-Menten constant                       |
| Vmax UGT1A1     | nmol/min/mg | 0.89 ([Kassahun 2007](#5-references))                      | Maximum rate of reaction                        |
| Km UGT1A9       | µM          | 296 ([Kassahun 2007](#5-references))                       | Michaelis-Menten constant                       |
| Vmax UGT1A9     | nmol/min/mg | 0.53 ([Kassahun 2007](#5-references))                      | Maximum rate of reaction                        |

### 2.2.2 Clinical data

A literature search was performed to collect available clinical data on Raltegravir in adults. 

The following publications were found in adults for model building and evaluation:

| Publication                       | Study description                                            |
| :-------------------------------- | :----------------------------------------------------------- |
| [Iwamoto 2008](#5-references)  | Single- and multiple-dose escalation study in healthy subjects |
| [Iwamoto 2009](#5-references)  | Effects of ritonavir and efavirenz on safety, tolerability and pharmacokinetics of raltegravir |
| [Markowitz 2006](#5-references) | Monotherapy, followed by a longer term combination therapy of raltegravir versus efavirenz |
| [Kassahun 2007](#5-references) | Pharmacokinetics study in healthy adults                     |
| [Rhee 2014](#5-references)     | Pediatric formulation study in healthy adults                |
| [Wenning 2009](#5-references)  | Effect of rifampin on the pharmacokinetics of raltegravir    |

## Literature in-vitro kinetics (given)

These are the measured in-vitro values from the literature - the fitted
catalytic rates are NOT given; you must determine them.

| molecule | process | parameter | value | unit |
| --- | --- | --- | --- | --- |
| UGT1A1 | MetabolizationLiverMicrosomes_MM | Km | 99.0 | µM |
| UGT1A9 | MetabolizationLiverMicrosomes_MM | Km | 296.0 | µM |
