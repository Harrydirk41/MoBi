# Building and evaluation of a PBPK model for fluvoxamine in healthy adults


> You are handed this report and the clinical data files. Extract the

> physicochemical and in-vitro inputs yourself, decide the model

> structure, and fit the parameters. The reference model's chosen

> methods and fitted values are withheld.


# Background

Fluvoxamine is a selective serotonin reuptake inhibitor used to treat major depression and obsessive compulsive disorder ([Perucca 1994](#5-references), [ANI Pharmaceuticals Inc. 2008](#5-references)) . Recommended doses are 50 to 300 mg once daily. The pharmacokinetics of orally administered single doses are linear. Following multiple oral administration, the pharmacokinetics at steady-state become non-linear, due to saturable Michaelis-Menten kinetics of the metabolic pathways ([Spigset 1998](#5-references)). Metabolism of fluvoxamine includes hydroxylation via CYP1A2 and O-demethylation via the very polymorphic CYP2D6 ([Miura 2007](#5-references), [Spigset 2001](#5-references)). Following oral administration fluvoxamine is excreted via the urine as metabolites ([DeBree 1983](#5-references)). The U.S. Food and Drug Administration (FDA) recommends fluvoxamine as strong clinical CYP1A2 and CYP2C19 index inhibitor to evaluate the impact of CYP1A2/CYP2C19 inhibition on CYP1A2/CYP2C19 substrates ([FDA 2017](#5-references)). Furthermore, the FDA lists fluvoxamine as moderate CYP3A4 inhibitor.

The aim of this project was to develop a PBPK model of fluvoxamine, mechanistically describing its metabolism by CYP1A2 and CYP2D6 and its inhibitory effect on CYP1A2 and CYP3A4, that can be used for drug-drug interaction (DDI) predictions.

The presented model was developed and evaluated by Britz et al. ([Britz 2019](#5-references))

## 2.2 Data<a id="22"></a>

### 2.2.1	In vitro / physico-chemical Data

A literature search was performed to collect available information on physicochemical properties of fluvoxamine. The obtained information from literature is summarized in the table below. 

| **Parameter**          | **Unit** | **Value**               | Source                                              | **Description**                                              |
| :--------------------- | -------- | ----------------------- | --------------------------------------------------- | ------------------------------------------------------------ |
| MW                     | g/mol    | 318.34                  | [Drugbank](#5-references)                           | Molecular weight                                             |
| pK<sub>a</sub>         |          | 9.40 (base)             | [Hallifax 2007](#5-references)                      | Acid dissociation constant                                   |
| Solubility (pH)        | mg/mL    | 14.66 (7.0)             | [MSDS](#5-references)                               | Solubility                                                   |
| f<sub>u</sub>          |          | 0.13 ± 0.01<sup>a</sup> | [Yao 2001](#5-references)                           | Fraction unbound in plasma                                   |
|                        |          | 0.14 ± 0.02<sup>a</sup> | [Yao 2001](#5-references)                           | Fraction unbound in plasma                                   |
|                        |          | 0.23                    | [Claassen 1983](#5-references)                      | Fraction unbound in plasma                                   |
| f<sub>u,mic</sub>      |          | 0.20 ± 0.05<sup>a</sup> | [Yao 2001](#5-references)                           | Fraction unbound in human liver microsomes at a protein concentration of 1 mg/mL |
|                        |          | 0.31 ± 0.03<sup>a</sup> | [Yao 2001](#5-references)                           | Fraction unbound in human liver microsomes at a protein concentration of 0.5 mg/mL |
|                        |          | 0.70 ± 0.03<sup>a</sup> | [Yao 2001](#5-references)                           | Fraction unbound in supersomes at a protein concentration of 0.3 mg/mL |
| CYP2D6 K<sub>m</sub>   | µmol/L   | 76.30                   | [Miura 2007](#5-references)                         | Michaelis-Menten constant                                    |
| CYP2D6 k<sub>cat</sub> | 1/min    | 0                       | [Crews 2014](#5-references)                         | The number of substrate molecule each enzyme site converts to product per unit time, and in which the enzyme is working at maximum efficiency |
| CYP1A2 K<sub>i</sub>   | µmol/L   | 0.011                   | [Karjalainen 2008](#5-references)                   | Competitive inhibition constant of the competitive inhibition model measured in human liver microsomes |
| CYP1A2 K<sub>i,u</sub> | nmol/L   | 35                      | [Yao 2001](#5-references)                           | Unbound competitive inhibition constant of the mixed inhibition model measured in human liver microsomes at a protein concentration of 1 mg/mL |
|                        | nmol/L   | 36                      | [Yao 2001](#5-references)                           | Competitive inhibition constant of the mixed inhibition model measured in human liver microsomes at a protein concentration of 0.5 mg/mL |
|                        | nmol/L   | 36                      | [Yao 2001](#5-references)                           | Competitive inhibition constant of the mixed inhibition model measured in supersomes at a protein concentration of 0.3 mg/mL |
| CYP3A4 K<sub>i</sub>   | µmol/L   | 1.60                    | [Olesen 2000](#5-references)                        | Competitive inhibition constant of the competitive inhibition model measured in human liver microsomes |

<sup>a</sup> denotes mean ± standard deviation

### 2.2.2 Clinical Data

A literature search was performed to collect available clinical data on fluvoxamine in healthy adults.

The fluvoxamine PBPK model was developed using 26 different clinical studies with pharmacokinetic (PK) blood sampling. These studies include 1 study of 30 mg fluvoxamine administered intravenously (iv) as a single-dose, and 25 studies of fluvoxamine administered orally (po) in single- or multiple-doses. In the single-dose po studies fluvoxamine was administered in doses of 25 - 200 mg. In the multiple-dose po studies fluvoxamine was administered once (q.d.) or twice daily (b.i.d.), in doses of 10 - 150 mg per administration.

#### 2.2.2.1	Model Building

The following studies were used for model building (training data):

| Publication                            | Arm / Treatment / Information used for model building        |
| :------------------------------------- | :----------------------------------------------------------- |
| [Japanese Society 2015](#5-references) | Healthy Japanese adults with 30 mg as 60 min infusion or oral administration of 200 mg |
| [de Vries 1993](#5-references)         | Healthy adults with oral administration of 25-100 mg         |
| [Orlando 2010](#5-references)          | Healthy adults with oral administration of 50 mg             |
| [Labellarte 2004](#5-references)       | Healthy CYP2D6 EM with oral administration of 50 mg twice a day |
| [Spigset 1998](#5-references)          | Healthy CYP2D6 EM (80%) and PM (20%) with oral administration of doses between 12.5-100 mg twice a day |
| [Fleishaker 1994](#5-references)       | Healthy adults with oral administration of 50 mg or 100 mg once daily |

#### 2.2.2.2	Model Verification

The following studies were used for model verification:

| Publication                            | Arm / Treatment / Information used for model building        |
| :------------------------------------- | :----------------------------------------------------------- |
| [Christensen 2002](#5-references)      | Healthy CYP2D6 EM with oral administration of 10 mg or 25 mg twice a day and healthy CYP2D6 PM with oral administration of 10 mg or 25 mg once daily |
| [Fukasawa 2006](#5-references)         | Healthy Japanese adults with single oral doses of 50 mg      |
| [Japanese Society 2015](#5-references) | Healthy Japanese adults with single oral doses of 25-100 mg  |
| [Kunii 2005](#5-references)            | Healthy CYP2D6 EM with single oral doses of 50 mg            |
| [Spigset 1995](#5-references)          | Healthy smokers or non-smokers with oral administration of 50 mg as single dose |
| [Spigset 1997](#5-references)          | Healthy CYP2D6 EM or PM with oral administration of 50 mg as single dose |
| [van Harten 1991](#5-references)       | Healthy adults  with oral administration of 50 mg as single dose |
| [de Vries 1992](#5-references)          | Healthy adults with oral administration of 50 mg twice a day |
| [Bahrami 2007](#5-references)          | Healthy adults with oral administration of 100 mg as single dose |
| [de Bree 1983](#5-references)          | Healthy adults with oral administration of 100 mg as single dose |

## Literature in-vitro kinetics (given)

These are the measured in-vitro values from the literature - the fitted
catalytic rates are NOT given; you must determine them.

| molecule | process | parameter | value | unit |
| --- | --- | --- | --- | --- |
| CYP1A2 | MetabolizationLiverMicrosomes_MM | In vitro Vmax for liver microsomes | 0.0 | pmol/min/mg mic. protein |
| CYP1A2 | MetabolizationLiverMicrosomes_MM | Km | 0.0073460807948 | µmol/l |
| CYP2D6 | rCYP450_MM | In vitro Vmax/recombinant enzyme | 0.69 | pmol/min/pmol rec. enzyme |
| CYP2D6 | rCYP450_MM | Km | 76.3 | µmol/l |
