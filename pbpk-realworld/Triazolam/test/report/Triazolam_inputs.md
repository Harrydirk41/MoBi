# Building and evaluation of a PBPK model for triazolam in healthy adults


> You are handed this report and the clinical data files. Extract the

> physicochemical and in-vitro inputs yourself, decide the model

> structure, and fit the parameters. The reference model's chosen

> methods and fitted values are withheld.


# Background

The presented model building and evaluation report evaluates the performance of a PBPK model for triazolam in healthy adults.

Triazolam, sold under the trade name Halcion, among others, belongs to the group of benzodiazepines and is used for short-term treatment of insomnia and circadian rhythm sleep disorders. It is generally administered orally as immediate release tablet, but other forms of administrations, e.g. intravenously or as sublingual tablet, exist as well.

Following oral administration, triazolam is rapidly absorbed with an absolute bioavailability of 44 ± 24% (mean ± standard deviation, [Kroboth 1995](#5-references)). Triazolam is widely distributed throughout the body. Its fraction unbound in human plasma averages around 17% and is, within the range of 20 to 1000 ng/mL, not influenced by total triazolam concentrations ([Eberts 1981](#5-references)). Triazolam is extensively metabolized via CYP3A4 to α-hydroxy-alprazolam and 4-hydroxy-alprazolam ([Eberts 1981](#5-references), [Kronbach 1989](#5-references)) and is therefore often used as victim compound in drug-drug interaction (DDI) studies.

The presented triazolam PBPK model was developed for intravenous (IV) administration and oral (PO) administration of the immediate release tablet given in fasted state in healthy, non-obese adults.

## 2.2 Data<a id="data"></a>

### 2.2.1 In vitro / physicochemical data

A literature search was carried out to collect available information on physicochemical properties of triazolam. The obtained information from the literature is summarized in the table below and is used for model building.

| **Parameter**          | **Unit** | **Literature**                                               | **Description**                                              |
| :--------------------- | -------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| Molecular weight       | g/mol    | 343.21 ([drugbank.ca](#5-references))                       | Molecular weight                                             |
| pK<sub>a</sub> (basic) |          | 1.52 ([Konishi 1982](#5-references))                        | Acid dissociation constant                                   |
| f<sub>u</sub>          |          | 0.099 ± 0.015<sup>a</sup> ([Jochemsen 1983](#5-references)); 0.11 ([Eberts 1981](#5-references)); 0.174 ± 0.020<sup>a</sup> ([Friedman 1988](#5-references)); 0.188 ± 0.139<sup>a</sup> ([Ochs 1987](#5-references)); 0.213 [0.193 - 0.264]<sup>b</sup> ([Greenblatt 1983b](#5-references)); 0.229 [0.204 - 0.259]<sup>c</sup> ([Greenblatt 1983b](#5-references)) | Fraction unbound in human plasma of healthy adults           |
| Water solubility       | mg/L     | 4.53 ([drugbank.ca](#5-references))                         | Estimated solubility in water                                |

<sup>a</sup> mean ± standard deviation

<sup>b</sup> mean [range] in young males

<sup>c</sup> mean [range] in young females

### 2.2.2 Clinical data

A literature search was carried out to collect triazolam PK data in healthy adults. 

The following publications were found and used for model building and evaluation:

| Publication                            | Study description                                            |
| :------------------------------------- | :----------------------------------------------------------- |
| [Friedman 1986](#5-references)        | PO single dose administration of 0.5 mg                      |
| [Friedman 1988](#5-references)        | PO single dose administration of 0.5 mg                      |
| [Greenblatt 1989](#5-references)      | PO single dose administration of 0.25 mg                     |
| [Greenblatt 1991](#5-references)      | PO single dose administration of 0.125 mg                    |
| [Greenblatt 2000](#5-references)      | PO single dose administration of 0.25 mg                     |
| [Greenblatt 2004](#5-references)      | PO single dose administration of 0.25 mg                     |
| [Hukkinen 1995](#5-references)        | PO single dose administration of 0.25 mg                     |
| [Lilja 2000](#5-references)           | PO single dose administration of 0.25 mg                     |
| [Kroboth 1985](#5-references)         | IV single dose administration of 0.25 mg and PO single dose administration of 0.25 mg |
| [O'Connor-Semmes 2001](#5-references) | PO single dose administration of 0.25 mg                     |
| [Ochs 1984](#5-references)            | PO single dose administration of 0.5 mg                      |
| [Phillips 1986](#5-references)        | PO single dose administration of 0.5 mg                      |
| [Smith 1987](#5-references)           | IV single dose administration of 0.125 mg, 0.25 mg, 0.5 mg, 0.75 mg, and 1 mg |
| [Varhe 1994](#5-references)           | PO single dose administration of 0.25 mg                     |
| [Varhe 1996a](#5-references)          | PO single dose administration of 0.25 mg                     |
| [Varhe 1996b](#5-references)          | PO single dose administration of 0.25 mg                     |
| [Varhe 1996c](#5-references)          | PO single dose administration of 0.25 mg                     |
| [Villikka 1997](#5-references)        | PO single dose administration of 0.5 mg                      |
| [Villikka 1998](#5-references)        | PO single dose administration of 0.5 mg                      |
| [von Moltke 1996](#5-references)      | PO single dose administration of 0.125 mg                    |

## Literature in-vitro kinetics (given)

These are the measured in-vitro values from the literature - the fitted
catalytic rates are NOT given; you must determine them.

| molecule | process | parameter | value | unit |
| --- | --- | --- | --- | --- |
| CYP3A4 | MetabolizationLiverMicrosomes_MM | In vitro Vmax for liver microsomes | 2.36 | nmol/min/mg mic. protein |
| CYP3A4 | MetabolizationLiverMicrosomes_MM | Km | 74.2 | µmol/l |
| CYP3A4 | MetabolizationLiverMicrosomes_MM | In vitro Vmax for liver microsomes | 10.27 | nmol/min/mg mic. protein |
| CYP3A4 | MetabolizationLiverMicrosomes_MM | Km | 305.0 | µmol/l |
