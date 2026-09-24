# Building and evaluation of a PBPK model for efavirenz in healthy adults


> You are handed this report and the clinical data files. Extract the

> physicochemical and in-vitro inputs yourself, decide the model

> structure, and fit the parameters. The reference model's chosen

> methods and fitted values are withheld.


# Background

Efavirenz is a non-nucleoside reverse transcriptase inhibitor (NNRTI) and is an antiretroviral drug to treat HIV.

Its major metabolizing enzyme is CYP2B6, but CYP3A4, CYP3A5, CYP1A2 and CYP2A6 also play a role ([Ward 2003](#5-references), [Ogburn 2010](#5-references)). CYP2B6 polymorphism is a major determinant of clinical efavirenz disposition and dose adjustment. Efavirenz activates the pregnane X receptor (PXR) and induces its target gene expression. As a consequence, some cytochrome P450 genes are upregulated, and, e.g. higher CYP3A4  ([Shou 2008](#5-references)) and CYP2B6 ([Ke 2016](#5-references)) activity levels can be measured.

It has a long half-life ranging from 52 to 76 hours following single oral doses and 40 to 55 hours following long term administration as a result of auto-induction of efavirenz metabolism. The long plasma half-life allows for once daily administration with long term administration of a single 600 mg daily dose ([Smith 2001](#5-references)). 

The presented efavirenz model was established using clinical PK data of 7 publications covering a dose range from 200 to 600 mg after single and multiple oral administration.  

The herein presented model building and evaluation report evaluates the performance of the PBPK model for efavirenz in (healthy) adults. 

The established efavirenz PBPK model is verified for the use as a perpetrator drug in drug-drug interaction simulations.

## 2.2 Data<a id="data"></a>

### 2.2.1 In vitro / physico-chemical Data

A literature search was performed to collect available information on physiochemical properties of efavirenz. The obtained information from literature is summarized in the table below. 

| **Parameter**   | **Unit** | **Value**       | Source                                                       | **Description**                                 |
| :-------------- | -------- | --------------- | ------------------------------------------------------------ | ----------------------------------------------- |
| MW              | g/mol    | 315.675         | https://www.drugbank.ca/                                     | Molecular weight                                |
| pK<sub>a</sub>  | 10.1     | (base)          | [Rabel 1996](#5-references)                                | Acid dissociation constant                      |
| Solubility (pH) | mg/L     | 11.5 (6.4) | [Cristofoletti 2013](#5-references)        | Water solubility                               |
| fu              |         | 0.006 [0.004 - 0.015] | [Almond 2005](#5-references) | Fraction unbound in plasma                      |
| Emax (CYP3A4) |          | 7.27, 3.15 (average 5.21) | [Shou 2008](#5-references) | Maximum induction effect |
| EC50 (CYP3A4) | µmol/l | 12.5, 2.18 (average 7.34) | [Shou 2008](#5-references) | Concentration at half maximum induction |
| Emax (CYP2B6) |          | 5.1       | [Ke 2016](#5-references) | Maximum induction effect |
| EC50 (CYP2B6) | µmol/l | 5.1       | [Ke 2016](#5-references) | Concentration at half maximum induction |

### 2.2.2 Clinical Data

A literature search was performed to collect available clinical data on efavirenz in healthy adults.

#### 2.2.2.1 Model Building

The following studies were used for model building:

| Publication                  | Arm / Treatment / Information used for model building        |
| :--------------------------- | :----------------------------------------------------------- |
| [Mouly 2002](#5-references)  | Healthy subjects receiving a single oral dose of 200 and 400 mg |
| [Ogburn 2013](#5-references) | Healthy subjects receiving a single oral dose of 600 mg      |
| [Xu 2013](#5-references)     | Healthy subjects with different CYP2B6 genotypes receiving a single oral dose of 600 mg |
| [Dooley 2012](#5-references) | Healthy subjects with different CYP2B6 genotypes receiving multiple doses of 600 mg |
| [Garg 2013](#5-references)   | Healthy subjects receiving multiple doses of 600 mg          |
| [Huang 2012](#5-references)  | Healthy subjects receiving multiple doses of 600 mg          |

#### 2.2.2.2 Midazolam interaction studies used to parameterize CYP3A4 interaction

The following studies were used for parameterization of CYP3A4 interaction:

| Publication                       | Arm / Treatment / Information used for model building        |
| :-------------------------------- | :----------------------------------------------------------- |
| [Mikus 2017](#5-references)       | Healthy subjects receiving a single oral dose of  400 mg Efavirenz at t=0h, 4 mg midazolam at t=12h and a single intravenous dose of 2 mg midazolam at t=18h. |
| [Katzenmaier 2010](#5-references) | Healthy subjects receiving multiple oral doses of 400 mg efavirenz QD. On day 14, subjects receive a single oral midazolam dose of 3 mg. |

## Literature in-vitro kinetics (given)

These are the measured in-vitro values from the literature - the fitted
catalytic rates are NOT given; you must determine them.

| molecule | process | parameter | value | unit |
| --- | --- | --- | --- | --- |
| CYP2B6 | rCYP450_MM | In vitro Vmax/recombinant enzyme | 3.5 | pmol/min/pmol rec. enzyme |
| CYP2B6 | rCYP450_MM | Km | 6.4 | µmol/l |
| CYP1A2 | rCYP450_MM | In vitro Vmax/recombinant enzyme | 0.6 | pmol/min/pmol rec. enzyme |
| CYP1A2 | rCYP450_MM | Km | 8.3 | µmol/l |
| CYP3A4 | rCYP450_MM | In vitro Vmax/recombinant enzyme | 0.16 | pmol/min/pmol rec. enzyme |
| CYP3A4 | rCYP450_MM | Km | 23.5 | µmol/l |
| CYP3A5 | rCYP450_MM | In vitro Vmax/recombinant enzyme | 0.6 | pmol/min/pmol rec. enzyme |
| CYP3A5 | rCYP450_MM | Km | 19.1 | µmol/l |
| CYP2A6 | rCYP450_MM | In vitro Vmax/recombinant enzyme | 1.0 | pmol/min/pmol rec. enzyme |
| CYP2A6 | rCYP450_MM | Km | 7.7 | µmol/l |
| CYP2B6 | rCYP450_MM | In vitro Vmax/recombinant enzyme | 2.268966 | pmol/min/pmol rec. enzyme |
| CYP2B6 | rCYP450_MM | Km | 6.4 | µmol/l |
| CYP2B6 | rCYP450_MM | In vitro Vmax/recombinant enzyme | 1.448276 | pmol/min/pmol rec. enzyme |
| CYP2B6 | rCYP450_MM | Km | 6.4 | µmol/l |
