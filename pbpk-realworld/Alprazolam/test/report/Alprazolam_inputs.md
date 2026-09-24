# Building and evaluation of a PBPK model for alprazolam in healthy adults


> You are handed this report and the clinical data files. Extract the

> physicochemical and in-vitro inputs yourself, decide the model

> structure, and fit the parameters. The reference model's chosen

> methods and fitted values are withheld.


# Background

The presented model building and evaluation report evaluates the performance of a PBPK model for alprazolam in healthy adults.

Alprazolam, sold under the trade names Xanax and Solanax, among others, belongs to the group of benzodiazepines and is commonly used in short term management of anxiety disorders. It is generally administered orally as immediate release or extended release tablet, but other forms are also available, e.g. solution or sublingual tablet.

Following oral administration, alprazolam is rapidly absorbed with an absolute bioavailability ranging from 80% to 100% ([Greenblatt 1993](#5-references)). Absorption is independent of the dose and the relative bioavailability of solid and liquid dosage forms has been observed to be similar ([Dawson 1984](#5-references)). Alprazolam is widely distributed throughout the body and its free fraction in plasma, averaging around 30%, is not influenced by total alprazolam concentrations within the tested range of 0.01 to 10 mg/L ([Moschitto 1983](#5-references)). Alprazolam is extensively metabolized to various metabolites ([von Moltke 1993](#5-references)). The two major metabolites, α-hydroxy-alprazolam and 4-hydroxy-alprazolam, are formed through oxidation catalyzed by CYP3A ([Eberts 1980](#5-references), [von Moltke 1993](#5-references)). Within 72 h of a 2 mg oral dose of <sup>14</sup>C-alprazolam, 20% of the dose have been observed to be excreted unchanged in urine ([Eberts 1980](#5-references)). Alprazolam displays dose linear pharmacokinetics and does not accumulate during multiple dose treatment ([Dawson 1984](#5-references), [Greenblatt 1993](#5-references)). Because of the predominant role of CYP3A4 in alprazolam elimination, alprazolam is often used as victim compound in drug-drug interaction (DDI) studies.

The presented alprazolam PBPK model was developed for intravenous (IV) administration and oral (PO) administration of the immediate release tablet (Xanax) or extended-release formulation (Solanax) given in fasted state in healthy, non-obese adults; administration in fed state was not addressed here.

## 2.2 Data<a id="data"></a>

### 2.2.1 In vitro / physicochemical data

A literature search was carried out to collect available information on physicochemical properties of alprazolam. The obtained information from the literature is summarized in the table below and is used for model building.

| **Parameter**             | **Unit** | **Literature**                                               | **Description**                                              |
| :------------------------ | -------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| Molecular weight          | g/mol    | 308.765 ([drugbank.ca](#5-references))                       | Molecular weight                                             |
| pK<sub>a</sub> (basic)    |          | 2.40 ([Cho 1983](#5-references), [Raymond 1986](#5-references)); 2.48 ± 0.01 ([Manchester 2018](#5-references)) | Acid dissociation constant                                   |
| f<sub>u</sub>             |          | 0.20 ([Eberts 1980](#5-references)); 0.233 ± 0.028<sup>a</sup> ([Schmith 1991](#5-references)); 0.270 ± 0.017<sup>a</sup> ([Scavone 1988](#5-references)); 0.284 ± 0.017<sup>a</sup> ([Scavone 1988](#5-references)); 0.290 ± 0.025<sup>a</sup> ([Juhl 1984](#5-references)); 0.298 [0.259 - 0.316]<sup>b</sup> ([Abernethy 1983](#5-references)); 0.311 ± 0.026<sup>a</sup> ([Ochs 1986](#5-references)); 0.316<sup>c</sup> ([Moschitto 1983](#5-references)) | Fraction unbound in human plasma of healthy adults           |
| Water solubility (pH 1.2) | mg/L     | 12 ([drugbank.ca](#5-references))                            | Estimated solubility in water at pH 1.2                      |
| Water solubility (pH 7.0) | mg/L     | 40 ([drugbank.ca](#5-references))                            | Estimated solubility in water at pH 7.0                      |
| Water solubility          | mg/L     | 73 ([Loftsson 2006](#5-references))                          | Experimentally measured solubility in water at 22°C - 24°C   |

<sup>a</sup> mean ± SD

<sup>b</sup> mean [range]

<sup>c</sup> mean

### 2.2.2 Clinical data

A literature search was carried out to collect alprazolam PK data in healthy adults. 

The following publications were found and used for model building and evaluation:

| Publication                            | Study description                                            |
| :------------------------------------- | :----------------------------------------------------------- |
| [Adams 1984](#5-references)            | IV single dose administration of 0.25 mg and 4 mg            |
| [Bertz 1997](#5-references)           | IV single dose administration of 2 mg (young subjects group) |
| [Eberts 1980](#5-references)          | PO single dose administration of 2 mg <sup>14</sup>C-alprazolam (no plasma concentration-time profile was reported, but the dose fraction excreted unchanged in urine was quantified) |
| [Eller 1990](#5-references)           | PO single dose administration of 1 mg (Treatment C: IR tablet in fasted state) |
| [Fleishaker 1989](#5-references)      | IV single dose administration of 1 mg (Treatment A)          |
| [Fleishaker 1994](#5-references)      | PO multiple dose administration of 1 mg four times daily at irregular time intervals for 4 days (Control phase) |
| [Friedman 1991](#5-references)        | PO single dose administration of 1 mg                        |
| [Greenblatt 1988](#5-references)      | PO single dose administration of 1 mg                        |
| [Greenblatt 1992](#5-references)      | PO single dose administration of 1 mg (Control phase)        |
| [Greenblatt 1998](#5-references)      | PO single dose administration of 1 mg (Trial A)              |
| [Greenblatt 2000](#5-references)      | PO single dose administration of 1 mg (Control group)        |
| [Juhl 1984](#5-references)            | PO single dose administration of 1 mg (Healthy control group) |
| [Kaplan 1998](#5-references)          | PO single dose administration of 1 mg (young subjects group) |
| [Kirkwood 1991](#5-references)        | PO single dose administration of 1 mg                        |
| [Kroboth 1988](#5-references)         | IV single dose administration of 0.5 mg, 1 mg followed by 72 µg over 8 h, and 2 mg |
| [Lin 1988](#5-references)             | IV single dose administration of 0.5 mg and PO single dose administration of 0.5 mg |
| [Schmider 1999](#5-references)        | PO single dose administration of 1 mg (Control phase)        |
| [Schmith 1991](#5-references)         | PO single dose administration of 0.5 mg and 2 mg (normal subjects group) |
| [Smith 1984](#5-references)           | IV single dose administration of 1 mg and PO single dose administration of 1 mg |
| [Venkatakrishnan 2005](#5-references) | IV single dose administration of 1 mg                        |
| [Wennerholm 2005](#5-references)      | PO single dose administration of 1 mg                        |
| [Yasui 1996](#5-references)           | PO single dose administration of 0.8 mg (Control phase)      |
| [Yasui 1998](#5-references)           | PO single dose administration of 0.8 mg (Control phase)      |
| [Yasui 2000](#5-references)           | PO single dose administration of 0.8 mg (Control phase)      |

## Literature in-vitro kinetics (given)

These are the measured in-vitro values from the literature - the fitted
catalytic rates are NOT given; you must determine them.

| molecule | process | parameter | value | unit |
| --- | --- | --- | --- | --- |
| CYP3A4 | MetabolizationLiverMicrosomes_MM | In vitro Vmax for liver microsomes | 0.131 | nmol/min/mg mic. protein |
| CYP3A4 | MetabolizationLiverMicrosomes_MM | Km | 269.0 | µmol/l |
| CYP3A4 | MetabolizationLiverMicrosomes_MM | In vitro Vmax for liver microsomes | 2.23 | nmol/min/mg mic. protein |
| CYP3A4 | MetabolizationLiverMicrosomes_MM | Km | 704.0 | µmol/l |
