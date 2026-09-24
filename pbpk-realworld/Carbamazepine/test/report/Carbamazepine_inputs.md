# Building and evaluation of a PBPK Model for carbamazepine in healthy adults


> You are handed this report and the clinical data files. Extract the

> physicochemical and in-vitro inputs yourself, decide the model

> structure, and fit the parameters. The reference model's chosen

> methods and fitted values are withheld.


# Background

Carbamazepine, sold under the trade name Tegretol<sup>®</sup> among others, is an anticonvulsant medication used primarily to treat epilepsy and neuropathic pain. Other indications include schizophrenia where it is used as an adjunctive treatment along with other medications, and bipolar disorder where it is used as a second-line agent. Carbamazepine is typically taken by mouth on empty stomach or together with meals, depending on the administered formulation. 

Carbamazepine is extensively metabolized by various enzymes including CYP2B6, 2C8, 3A4, and UGT2B7 ([Kerr 1994](#5-references), [Pelkonen 2001](#5-references), [Staines 2004](#5-references)). Following oral administration the major dose fraction is metabolized to carbamazepine-10,11-epoxide ([Eichelbaum 1985](#5-references), [Tomson 1983](#5-references)). This reaction is mainly catalyzed by CYP3A4, with some contribution from CYP2C8 ([Kerr 1994](#5-references)). After oral administration, a minor fraction of the dose (approximately 1 - 3%) is excreted unchanged in urine ([Bernus 1994](#5-references), [Morselli 1975](#5-references)), while approximately 1% of the dose can be recovered as unchanged drug in the bile ([Terhaag 1978](#5-references)).

Carbamazepine is classified by the U.S. Food and Drug Administration (FDA) as a strong CYP3A4 and CYP2B6 inducer and hence induces its own metabolism.

The herein presented model was developed independently of the model reported by Fuhr et al. ([Fuhr 2021](#5-references)). The main difference between the two models pertains to the metabolite carbamazepine-10,11-epoxide, which is included as separate compound in the model by Fuhr et al. ([Fuhr 2021](#5-references)), but not modeled in the herein presented model. Another structural model differences concerns the enzymatic elimination pathways of carbamazepine; the model by Fuhr et al. ([Fuhr 2021](#5-references)) includes five different metabolism pathways, whereas the herein presented model includes three different metabolism pathways. Additionally, the parameterization of CYP2B6 and 3A4 induction differs between the two models.

## 2.2 Data<a id="data"></a>

### 2.2.1 In vitro / physicochemical Data

A literature search was performed to collect available information on physicochemical properties of carbamazepine. The information is summarized in the table below. 

| **Parameter**                            | **Unit**                 | **Value**                                 | Source                            | **Description**                                         |
|:---------------------------------------- | ------------------------ | ----------------------------------------- | --------------------------------- | ------------------------------------------------------- |
| MW                                       | g/mol                    | 236.27                                    | [DrugBank DB00564](#5-references) | Molecular weight                                        |
| Solubility (pH)                          | µg/mL                    | 336 (6.2)                                 | [Annaert 2010](#5-references)     | Solubility in human intestinal fluid                    |
| Solubility (pH)                          | µg/mL                    | 283 (7.0)                                 | [Söderlind 2010](#5-references)   | Solubility in human intestinal fluid                    |
| Solubility (pH)                          | µg/mL                    | 306 (6.9)                                 | [Clarysse 2011](#5-references)    | Solubility in fasted human intestinal fluid             |
| f<sub>u</sub>                            |                          | 0.25                                      | [Pynnönen 1977](#5-references)    | Fraction unbound in plasma of healthy subjects          |
| f<sub>u</sub>                            |                          | 0.243 ± 0.013 [0.225 - 0.258]<sup>a</sup> | [Morselli 1975](#5-references)    | Fraction unbound in plasma of healthy male subjects     |
| f<sub>u</sub>                            |                          | 0.239                                     | [Di Salle 1974](#5-references)    | Fraction unbound in plasma of normal subjects           |
| f<sub>u</sub>                            |                          | 0.237 ± 0.031<sup>b</sup>                 | [Vinçon 1987](#5-references)      | Fraction unbound in plasma of epileptic patients        |
| f<sub>u</sub>                            |                          | 0.182 ± 0.05 [0.103 - 0.297]<sup>a</sup>  | [Hooper 1975](#5-references)      | Fraction unbound in plasma of normal subjects           |
| K<sub>m</sub> CYP2B6                     | µM                       | 420                                       | [Pearce 2002](#5-references)      | CYP2B6 Michaelis-Menten constant                        |
| V<sub>max</sub> CYP2B6                   | pmol/min/pmol rec enzyme | 0.429                                     | [Pearce 2002](#5-references)      | in vitro metabolic rate constant for recombinant CYP2B6 |
| K<sub>m</sub> CYP2C8                     | µM                       | 757                                       | [Cazali 2003](#5-references)      | CYP2C8 Michaelis-Menten constant                        |
| V<sub>max</sub> CYP2C8                   | pmol/min/pmol rec enzyme | 0.673                                     | [Cazali 2003](#5-references)      | in vitro metabolic rate constant for recombinant CYP2C8 |
| K<sub>m</sub> CYP3A4<sup>c</sup>         | µM                       | 282                                       | [Pearce 2002](#5-references)      | CYP3A4 Michaelis-Menten constant                        |
| K<sub>m</sub> CYP3A4 (→CBZE)<sup>d</sup> | µM                       | 248                                       | [Huang 2004](#5-references)       | CYP3A4 Michaelis-Menten constant                        |
| K<sub>m</sub> UGT2B7                     | µM                       | 214                                       | [Staines 2004](#5-references)     | UGT2B7 Michaelis-Menten constant                        |
| V<sub>max</sub> UGT2B7                   | pmol/min/mg mic enzyme   | 0.79                                      | [Staines 2004](#5-references)     | in vitro metabolic rate constant for microsomal enzymes |
| Microsomal UGT2B7                        | pmol/mg mic protein      | 82.9                                      | [Achour 2014](#5-references)      | Content of UGT2B7 proteins in liver microsomes          |
| Intestinal permeability                  | cm/min                   | 0.0258                                    | [Lennernäs 2007](#5-references)   | Transcellular intestinal permeability                   |

<sup>a</sup> denotes mean ± standard deviation [range]

<sup>b</sup> denotes mean ± standard deviation

<sup>c</sup> refers to CYP3A4-mediated reaction forming other metabolites than carbamazepine-10,11-epoxide

<sup>d</sup> refers to CYP3A4-mediated reaction forming carbamazepine-10,11-epoxide

### 2.2.2 Clinical Data

A literature search was conducted to collect available data on carbamazepine pharmacokinetics in healthy adult subjects after intravenous or oral administration in the fasted state.

The following studies were used for model building:

| Publication                    | Arm / Treatment / Information used for model building                                                                                                                    |
|:------------------------------ |:------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| [Bernus 1994](#5-references)   | Healthy subjects receiving two oral doses of 600 mg carbamazepine as IR tablet (only pharmacokinetic data following the first dose were used for model building)         |
| [Gérardin 1976](#5-references) | Healthy subjects receiving a single oral dose of 100 mg carbamazepine as IR tablet                                                                                       |
| [Gérardin 1990](#5-references) | Healthy subjects receiving a single oral dose of 100 mg [<sup>15</sup>N]-carbamazepine as suspension concomitantly with a single intravenous dose of 10 mg carbamazepine |
| [McLean 2001](#5-references)   | Healthy subjects receiving a single oral dose of 400 mg carbamazepine as XR formulation in fasted state                                                                  |
| [Møller 2001](#5-references)   | Healthy subjects receiving a multiple oral doses of carbamazepine, starting at 100 mg and escalating to 400 mg                                                           |
| [Wada 1978](#5-references)     | Healthy subjects receiving a single oral dose of 200 mg carbamazepine as syrup and IR tablet                                                                             |

The following studies were used for model evaluation:

| Publication                                                 | Arm / Treatment / Information used for model building                                                                                                  |
|:----------------------------------------------------------- |:------------------------------------------------------------------------------------------------------------------------------------------------------ |
| [Barzaghi 1987](#5-references)                              | Healthy subjects receiving a single oral dose of 400 mg carbamazepine                                                                                  |
| [Bedada 2015](#5-references)                                | Healthy subjects receiving a single oral dose of 200 mg carbamazepine                                                                                  |
| [Bedada 2016](#5-references)                                | Healthy subjects receiving a single oral dose of 200 mg carbamazepine                                                                                  |
| [Bernus 1994](#5-references)                                | Healthy subjects receiving two oral doses of 600 mg carbamazepine (only pharmacokinetic data following the second dose were used for model evaluation) |
| [Bianchetti 1987](#5-references)                            | Healthy subjects receiving a single oral dose of 400 mg carbamazepine                                                                                  |
| [Burstein 2000](#5-references)                              | Healthy subjects receiving a multiple oral doses of carbamazepine, starting at 100 mg and escalating to 400 mg                                         |
| [Caraco 1995](#5-references)                                | Healthy lean subjects receiving a single oral dose of 200 mg carbamazepine                                                                             |
| [Cawello 2000](#5-references)                               | Healthy subjects receiving a multiple oral doses of carbamazepine, starting at 100 mg and escalating to 200 mg                                         |
| [Cotter 1977](#5-references)                                | Healthy subject receiving a single oral dose of 800 mg carbamazepine                                                                                   |
| [Dalton 1985a](#5-references)                               | Healthy subjects receiving a single oral dose of 600 mg carbamazepine                                                                                  |
| [Dalton 1985b](#5-references)                               | Healthy subjects receiving a single oral dose of 600 mg carbamazepine                                                                                  |
| [Eichelbaum 1985](#5-references)                            | Healthy subjects receiving a single oral dose of 200 mg carbamazepine                                                                                  |
| [Elqidra 2004](#5-references)                               | Healthy subjects receiving a single oral dose of 200 mg carbamazepine                                                                                  |
| [European Patent Application EP 1044681 A2](#5-references)  | Healthy subjects receiving a single oral dose of 400 and 600 mg carbamazepine                                                                          |
| [Gérardin 1976](#5-references)                              | Healthy subjects receiving a single oral dose of 200, and 600 mg carbamazepine                                                                         |
| [Ji 2008](#5-references)                                    | Healthy subjects receiving a multiple oral doses of carbamazepine, starting at 200 mg and escalating to 400 mg                                         |
| [Kayali 1994](#5-references)                                | Healthy subjects receiving a single oral dose of 200 mg carbamazepine                                                                                  |
| [Kim 2005](#5-references)                                   | Healthy subjects receiving a single oral dose of 200 mg carbamazepine                                                                                  |
| [Kovacević 2009](#5-references)                             | Healthy subjects receiving a single oral dose of 400 mg carbamazepine                                                                                  |
| [Levy 1975](#5-references)                                  | Healthy subjects receiving a single oral carbamazepine dose of 6 mg/kg body weight                                                                     |
| [Meyer 1996](#5-references)                                 | Healthy subjects receiving a single oral dose of 200 mg carbamazepine                                                                                  |
| [Meyer 1998](#5-references)                                 | Healthy subjects receiving a single oral dose of 200 mg carbamazepine                                                                                  |
| [Miles 1989](#5-references)                                 | Healthy subjects receiving a multiple oral doses of 300 and 400 mg carbamazepine                                                                       |
| [Morselli 1975](#5-references)                              | Healthy subjects receiving a single oral dose of 400 mg carbamazepine                                                                                  |
| [Pynnönen 1977](#5-references)                              | Healthy subjects receiving a single oral dose of 400 mg carbamazepine                                                                                  |
| [Rawlins 1975](#5-references)                               | Healthy subject receiving a single oral dose of 50, 100, and 200 mg carbamazepine                                                                      |
| [Saint-Salvi 1987](#5-references)                           | Healthy subjects receiving a single oral dose of 200 mg carbamazepine                                                                                  |
| [Stevens 1998](#5-references)                               | Healthy subjects receiving multiple oral doses of 400 mg carbamazepine                                                                                 |
| [Strandjord 1975](#5-references)                            | Healthy subjects receiving a single oral dose of 400 mg carbamazepine                                                                                  |
| [Sumi 1987](#5-references)                                  | Healthy subjects receiving a single oral dose of 200 mg carbamazepine                                                                                  |
| [Tomson 1983](#5-references)                                | Healthy subject receiving a single oral doses of 200 mg carbamazepine                                                                                  |
| [US Patent Application - US 2009/0169619 A1](#5-references) | Healthy subjects receiving a single oral dose of 300 mg carbamazepine                                                                                  |
| [Wong 1983](#5-references)                                  | Healthy subjects receiving a single oral dose of 400 mg carbamazepine                                                                                  |

## Literature in-vitro kinetics (given)

These are the measured in-vitro values from the literature - the fitted
catalytic rates are NOT given; you must determine them.

| molecule | process | parameter | value | unit |
| --- | --- | --- | --- | --- |
| UGT2B7 | MetabolizationLiverMicrosomes_MM | In vitro Vmax for liver microsomes | 3.5 | pmol/min/mg mic. protein |
| UGT2B7 | MetabolizationLiverMicrosomes_MM | Km | 234.0 | µmol/l |
| CYP2B6 | MetabolizationLiverMicrosomes_MM | In vitro Vmax for liver microsomes | 49.0 | pmol/min/mg mic. protein |
| CYP2B6 | MetabolizationLiverMicrosomes_MM | Km | 235.0 | µmol/l |
| CYP3A4 | MetabolizationLiverMicrosomes_MM | In vitro Vmax for liver microsomes | 726.0 | pmol/min/mg mic. protein |
| CYP3A4 | MetabolizationLiverMicrosomes_MM | Km | 808.0 | µmol/l |
