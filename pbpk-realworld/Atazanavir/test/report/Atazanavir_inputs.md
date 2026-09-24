# Building and evaluation of a PBPK model for atazanavir in healthy adults


> You are handed this report and the clinical data files. Extract the

> physicochemical and in-vitro inputs yourself, decide the model

> structure, and fit the parameters. The reference model's chosen

> methods and fitted values are withheld.


# Background

The presented model building and evaluation report evaluates the performance of a PBPK model for atazanavir in healthy adults.

Atazanavir, sold under the trade name Reyataz among others, is an azapeptide protease inhibitor and used as antiretroviral medication to treat and prevent HIV/AIDS. It is taken orally once a day at a dose of 300 mg, if co-administered with ritonavir 100 mg orally once a day, and 400 mg, if administered without ritonavir. 

After oral administration, atazanavir is rapidly absorbed. A positive food effect has been observed, atazanavir is recommended to be taken with food. Protein binding is relatively high (86%) and independent of the concentration of serum proteins ([US Food and Drug Administration 2002](#5-references)). Atazanavir undergoes extensive metabolism by CYP3A isoenzymes with a dose fraction excreted unchanged in urine of approximately 7% ([US Food and Drug Administration 2002](#5-references), [Le Tiec 2005](#5-references)). Previous in vitro studies suggest that atazanavir is a mechanism-based inhibitor of CYP3A ([US Food and Drug Administration 2002](#5-references), [Perloff 2005](#5-references)) as well as a competitive inhibitor of CYP1A2, CYP2C9 and UGT1A1 ([US Food and Drug Administration 2002](#5-references), [Zhang 2005](#5-references)).

## 2.2 Data<a id="data"></a>

### 2.2.1 In vitro / physicochemical data

A literature search was carried out to collect available information on physicochemical properties of atazanavir. The obtained information from the literature is summarized in the table below and is used for model building.

| **Parameter**          | **Unit** | **Value**                                   | **Description**                                       |
| :--------------------- | -------- | ------------------------------------------------------------ | ----------------------------------------------------- |
| Molecular weight       | g/mol    | 704.9 ([drugbank.ca](#5-references))                         | Molecular weight                                      |
| pK<sub>a</sub> (basic) |          | 4.7 ([Berlin 2015](#5-references))                           | Acid dissociation constant                            |
| f<sub>u</sub>          |          | 0.14 ([US Food and Drug Administration 2002](#5-references)) | Fraction unbound in human plasma                      |
| Solubililty in FaSSIF  | µg/mL    | 2.74 ([Berlin 2015](#5-references))                          | Solubility in Fasted State Simulated Intestinal Fluid |
| Solubililty in FeSSIF  | µg/mL    | 4.13 ([Berlin 2015](#5-references))                          | Solubility in Fed State Simulated Intestinal Fluid    |

With regard to UGT1A1 inhibition, atazanavir inhibited 17β-Estradiol glucuronidation in recombinant UGT1A1 by a mixed-type mechanism (in-house data, [Jungmann 2019](#5-references)):

| **Parameter**    | **Unit** | **Value** | Source                         | **Description**                      |
| :--------------- | -------- | --------- | ------------------------------ | ------------------------------------ |
| K<sub>i</sub>    | µmol/L   | 0.22      | [Jungmann 2019](#5-references) | Inhibition constant                  |
| Alpha            |          | 4.5       | [Jungmann 2019](#5-references) | Alpha value in mixed-type inhibition |
| fu<sub>mic</sub> |          | 0.863     | [Fricke 2020](#5-references)   | determined *in vitro* at 0.22 µmol/L |

### 2.2.2 Clinical data

A literature search was carried out to collect available PK data on atazanavir in healthy adults. 

The following publications were found and used for model building and evaluation:

| Publication                                           | Study description                                            |
| :---------------------------------------------------- | :----------------------------------------------------------- |
| [Acosta 2007](#5-references)                          | 300 mg atazanavir BID, Period 1                              |
| [Agarwala 2003](#5-references)                        | 400 mg atazanavir QD, Day 6                                  |
| [Agarwala 2005a](#5-references)                       | 400 mg atazanavir QD, 400 mg AM                              |
| [Agarwala 2005b](#5-references)                       | 400 mg atazanavir QD, 400 mg (Treatment A)                   |
| [Martin 2008](#5-references)                          | 400 mg atazanavir QD, monotherapy                            |
| [Zhu 2010](#5-references)                             | 300 mg atazanavir QD                                         |
| [Zhu 2011](#5-references)                             | 400 mg atazanavir QD, 400 mg QPM and QAM                     |
| [US Food and Drug Administration 2002](#5-references) | Study AI424-004 (p. 94): 400 mg atazanavir single dose (treatment A);<br />Study AI424-014 (p. 77): 400 mg atazanavir single dose (young females & males);<br />Study AI424-015 (p. 81): 400 mg atazanavir single dose (normal subjects);<br />Study AI424-028 (p. 128): 200, 400, and 800 mg atazanavir QD (A-D Day6);<br />Study AI424-029 (p. 47): 400 mg [<sup>14</sup>C]atazanavir single dose;<br />Study AI424-040 (p. 64): 200, 400, and 800 mg atazanavir QD;<br />Study AI424-056 (p. 134): 300 mg atazanavir QD (without ritonavir, Day 10);<br />Study AI424-076 (p. 178): 400 and 800 mg atazanavir QD |
