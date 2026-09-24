# Building and evaluation of a PBPK model for propofol in adults


> You are handed this report and the clinical data files. Extract the

> physicochemical and in-vitro inputs yourself, decide the model

> structure, and fit the parameters. The reference model's chosen

> methods and fitted values are withheld.


# Background

The presented model building and evaluation report evaluates the performance of a PBPK model for propofol in adults.

Propofol is an anaesthetic agent used for induction and maintenance of general anaesthesia. Propofol is only given intravenously and is mainly metabolized by Uridine 5'-diphospho-glucuronosyltransferase 1A9 (UGT1A9) (53-70%) [(Al-Jahdari 2006, Restrepo 2009](#5-references)). The final propofol model features metabolism by UGT1A9 and to a minor extent by Cytochrome P450 2B6 (CYP2B6) [(Al-Jahdari 2006, Oda 2009](#5-references)). Additionally, there is excretion via glomerular filtration. The model adequately describes the pharmacokinetics of propofol in adults.

The propofol model is a whole-body PBPK model, allowing for dynamic translation between individuals with organs expressing UGT1A9. The propofol report demonstrates the level of confidence in the propofol PBPK model build with the OSP suite with regard to reliable predictions of propofol PK adults during model-informed drug development.

## 2.2 Data used<a id="data"></a>

### 2.2.1 In vitro / physicochemical data

A literature search was performed to collect available information on physicochemical properties of propofol. The obtained information from literature is summarized in the table below, and is used for model building.

| **Parameter**   | **Unit**    | **Literature value (reference)**        | **Description**                                 |
| :-------------- | ----------- | --------------------------------------- | ----------------------------------------------- |
| MW              | g/mol       | 178.2707 ([Drugbank.ca](#5-references)) | Molecular weight                                |
| pKa             |             | 10.1 ([Drugbank.ca](#5-references))     | Acid dissociation constant                      |
| Solubility (pH) | mg/L        | 124 (7) ([Drugbank.ca](#5-references))  | Solubility                                      |
| fu              |             | 0.024 ([Takizawa 2005](#5-references))  | Fraction unbound                                |
| Km,u UGT1A9     | mM          | 0.12 ([Al-Jahdari 2006](#5-references)) | Unbound Michaelis-Menten constant               |
| Vmax UGT1A9     | nmol/min/mg | 2.40 ([Al-Jahdari 2006](#5-references)) | Maximum rate of reaction                        |
| Km,u CYP2B7     | mM          | 0.03 ([Al-Jahdari 2006](#5-references)) | Unbound Michaelis-Menten constant               |
| Vmax CYP2B7     | nmol/min/mg | 1.08 ([Al-Jahdari 2006](#5-references)) | Maximum rate of reaction                        |

### 2.2.2 Clinical data

A literature search was performed to collect available clinical data on propofol in adults. 

The following publications were used in adults for model building and evaluation, of which individual patient data was available for download under [http://opentci.org/data/propofol](http://opentci.org/data/propofol):

| Publication                       | Study description                                            |
| :-------------------------------- | :----------------------------------------------------------- |
| [Gepts 1987](#5-references) | Disposition of propofol administered as constant rate intravenous infusion in humans |
| [Schnider 1998](#5-references) | Influence of administration rate on propofol plasma-effect site equilibrium |
| [Struys 2007](#5-references) | The influence of method of administration and covariates on the pharmacokinetics of propofol in adult volunteers |

## Literature in-vitro kinetics (given)

These are the measured in-vitro values from the literature - the fitted
catalytic rates are NOT given; you must determine them.

| molecule | process | parameter | value | unit |
| --- | --- | --- | --- | --- |
| UGT1A9 | MetabolizationLiverMicrosomes_MM | Km | 0.52 | mM |
| CYP2B6 | MetabolizationLiverMicrosomes_MM | Km | 0.02 | mM |
