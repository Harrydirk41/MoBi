# Building and Evaluation of a PBPK Model for Verapamil in Adults


> You are handed this report and the clinical data files. Extract the

> physicochemical and in-vitro inputs yourself, decide the model

> structure, and fit the parameters. The reference model's chosen

> methods and fitted values are withheld.


# Background

Verapamil is used for the treatment of high blood pressure, angina (chest pain from not enough blood flow to the heart), and supraventricular tachycardia.

Verapamil is administered as a 1:1 racemat of R- and S-verapamil which are metabolized mainly by CYP3A4 to S- and R-norverapamil. All four entities are mechanism-based inactivators of CYP3A4 and non-competitive inhibitors of P-gp. 

The presented verapamil model was established using observed concentration-time profiles of more than 45 clinical studies with doses from 0.1 mg to 250 mg in different verapamil dosing schedules including multiple doses and different routes of administration (intravenous, single and multiple oral administration). It includes enantioselective plasma protein binding, enantioselective metabolism by CYP3A4, non-stereospecific P-gp transport, and passive glomerular filtration 
The model building and application has been published by Hanke *et al.* 2020 ([Hanke 2020](#5-references)). 

The herein presented model building and evaluation report evaluates the performance of the PBPK model for verapamil in (healthy) adults.

## 2.2 Data<a id="methods-data"></a>

### In vitro / physico-chemical Data

A literature search was performed to collect available information on physiochemical properties of R- and S-verapamil and R- and S-norverapamil. The obtained information from literature is summarized in the tables below. 

#### R-verapamil

| **Parameter**   | **Unit** | **Value**       | Source                                                       | **Description**                                 |
| :-------------- | -------- | --------------- | ------------------------------------------------------------ | ----------------------------------------------- |
| MW              | g/mol    | 454.611     | [Wishart 2006](#5-references)                                    | Molecular weight                                |
| pK<sub>a</sub> (base)  | -  | 8.75          | [Hasegawa 1984](#5-references)              | Acid dissociation constant                      |
| Solubility (pH 6.54) | g/L     | 46.0 | [Vogelpoel 2004](#5-references)    | Water solubility                               |
| fu              |  %       | 5.1 | [Sanaee 2011](#5-references) | Fraction unbound in plasma                      |
| CYP3A4 Km -> Norvera            | µmol/L | 19.59 | [Wang 2013](#5-references) | CYP3A4 Michaelis-Menten constant for norverapamil formation    |
| CYP3A4 Km -> D617		  | µmol/L | 35.34 | [Wang 2013](#5-references) | CYP3A4 Michaelis-Menten constant for D617 formation    |
| P-gp Km		  | µmol/L | 1.01 | [Shirasaka 2008](#5-references) | Pgp Michaelis-Menten constant    |
| CYP3A4 MBI KI		  | µmol/L | 27.63 | [Wang 2013](#5-references) | Conc. for half-maximal inactivation    |
| CYP3A4 MBI kinact		  | 1/min | 0.038 | [Wang 2013](#5-references) | Maximum inactivation rate    |
| Pgp non-competitive Ki		  | µmol/L | 0.31 | [Döppenschmitt 1999](#5-references) | Conc. for half-maximal inactivation    |

#### S-verapamil

| **Parameter**   | **Unit** | **Value**       | Source                                                       | **Description**                                 |
| :-------------- | -------- | --------------- | ------------------------------------------------------------ | ----------------------------------------------- |
| MW              | g/mol    | 454.611     |[Wishart 2006](#5-references)                                    | Molecular weight                                |
| pK<sub>a</sub> (base)  | -  | 8.75          | [Hasegawa 1984](#5-references)            | Acid dissociation constant                      |
| Solubility (pH 6.54) | g/L     | 46.0 | [Vogelpoel 2004](#5-references) | Water solubility                               |
| fu              |  %       | 11 | [**Sanaee 2011**](#5-references) | Fraction unbound in plasma                      |
| CYP3A4 Km -> Norvera            | µmol/L | 9.72 | [Wang 2013](#5-references) | CYP3A4 Michaelis-Menten constant for norverapamil formation    |
| CYP3A4 Km -> D617		  | µmol/L | 23.64 | [Wang 2013](#5-references) | CYP3A4 Michaelis-Menten constant for D617 formation    |
| P-gp Km		  | µmol/L | 1.01 | [Shirasaka 2008](#5-references) | Pgp Michaelis-Menten constant    |
| CYP3A4 MBI KI		  | µmol/L | 3.85 | [Wang 2013](#5-references) | Conc. for half-maximal inactivation    |
| CYP3A4 MBI kinact		  | 1/min | 0.034 | [Wang 2013](#5-references) | Maximum inactivation rate    |
| Pgp non-competitive Ki		  | µmol/L | 0.31 | Döppenschmitt 1999 (#5-references) | Conc. for half-maximal inactivation    |

#### R-norverapamil

| **Parameter**   | **Unit** | **Value**       | Source                                                       | **Description**                                 |
| :-------------- | -------- | --------------- | ------------------------------------------------------------ | ----------------------------------------------- |
| MW              | g/mol    | 440.584     | [Wishart 2006](#5-references)                                    | Molecular weight                                |
| pK<sub>a</sub> (base)  | -  | 8.6 - 8.9          | [Sigma-Aldrich 2013](#5-references)             | Acid dissociation constant                      |
| fu              |  %       | 5.1 | assumed (from parent) | Fraction unbound in plasma                      |
| CYP3A4 Km -> D620            | µmol/L | 144 | [Tracy 1999](#5-references) | CYP3A4 Michaelis-Menten constant for norverapamil degradation |
| P-gp Km		  | µmol/L | 1.01 | assumed (from parent) | Pgp Michaelis-Menten constant    |
| CYP3A4 MBI KI		  | µmol/L | 6.1 | [Wang 2013](#5-references) | Conc. for half-maximal inactivation    |
| CYP3A4 MBI kinact		  | 1/min | 0.048 | [Wang 2013](#5-references) | Maximum inactivation rate    |
| Pgp non-competitive Ki		  | µmol/L | 0.30 | [Pauli-Magnus 2000](#5-references) | Conc. for half-maximal inactivation    |

#### S-norverapamil

| **Parameter**   | **Unit** | **Value**       | Source                                                       | **Description**                                 |
| :-------------- | -------- | --------------- | ------------------------------------------------------------ | ----------------------------------------------- |
| MW              | g/mol    | 440.584     | [Wishart 2006](#5-references)                                    | Molecular weight                                |
| pK<sub>a</sub> (base)  | -  | 8.6 - 8.9          | [Sigma-Aldrich 2013](#5-references)            | Acid dissociation constant                      |
| fu              |  %       | 11 | assumed (from parent) | Fraction unbound in plasma                      |
| CYP3A4 Km -> D620            | µmol/L | 36 | [Tracy 1999](#5-references) | CYP3A4 Michaelis-Menten constant for norverapamil degradation |
| P-gp Km		  | µmol/L | 1.01 | assumed (from parent) | Pgp Michaelis-Menten constant    |
| CYP3A4 MBI KI		  | µmol/L | 2.90 | [Wang 2013](#5-references) | Conc. for half-maximal inactivation    |
| CYP3A4 MBI kinact		  | 1/min | 0.080 | [Wang 2013](#5-references) | Maximum inactivation rate    |
| Pgp non-competitive Ki		  | µmol/L | 0.30 | [Pauli-Magnus 2000](#5-references) | Conc. for half-maximal inactivation    |

#### Individual

| **Parameter**   | **Unit** | **Value**       | Source                                                       | **Description**                                 |
| :-------------- | -------- | --------------- | ------------------------------------------------------------ | ----------------------------------------------- |
| EHC continuous fraction           |    | 1    | Assumption                                | EHC continuous fraction                             |

### Clinical Data

A literature search was performed to collect available clinical data on verapamil in healthy adults.

#### Model Building and parameterizing of CYP3A4 interaction

The following studies were used for model building and parameterization of CYP3A4 interaction:
If not stated otherwise, the drug was given as a 1:1 racemat of S- and R-verapamil.

| Publication                      | Arm / Treatment / Information used for model building        |
| :------------------------------- | :----------------------------------------------------------- |
| [Eichelbaum 1984](#5-references) | Healthy subjects receiving single intravenous doses of 5, 25 and 50 mg of R-verapamil and 5, 7.5 and 10 mg of S-verapamil |
| [Streit 2005](#5-references) | Healthy subjects receiving single intravenous doses of 5 mg |
| [Barbarash 1988](#5-references)  | Healthy subjects receiving single intravenous doses of 10 mg |
| [Abernethy 1993](#5-references) | Healthy subjects receiving single intravenous doses of 20 mg |
| [Maeda 2011](#5-references) | Healthy subjects receiving single oral doses of 0.1, 3 and 16 mg |
| [Blume 1989](#5-references) | Healthy subjects receiving single oral doses of 80 mg |
| [Ratiopharm 1988](#5-references) | Healthy subjects receiving single oral doses of 80 mg |
| [Johnson 2001](#5-references)   | Healthy subjects receiving multiple oral doses of 80 mg TID. | 
| [Härtter 2012](#5-references) | Healthy subjects receiving single oral doses of 120 mg and multiple oral doses of 120 mg BID |
| [Hla 1987](#5-references)       | Healthy subjects receiving multiple oral doses of 120 mg BID |

#### Model verification 

The following studies were used for model verification:

| Publication                      | Arm / Treatment / Information used for model building        |
| :------------------------------- | :----------------------------------------------------------- |
| [Mooy 1985](#5-references) |  Healthy subjects receiving single intravenous doses of 3 mg and single oral doses of 80 mg |
| [Johnston 1981](#5-references)   | Healthy subjects receiving single intravenous doses of 0.1 mg/kg and single oral doses of 120 mg |
| [Abernethy 1985](#5-references) | Healthy subjects receiving single intravenous doses of 10 mg and single oral doses of 120 mg |
| [Barbarash 1988](#5-references)  | Healthy subjects receiving single oral doses of 120 mg |
| [Wing 1985](#5-references) | Healthy subjects receiving single intravenous doses 10mg and single oral doses of 80 mg |
| [McAllister 1982](#5-references) | Healthy subjects receiving single intravenous doses of 10 mg |
| [Smith 1984](#5-references) | Healthy subjects receiving single intravenous doses of 10 mg and single oral doses of 120 mg |
| [Freedman 1981](#5-references) | Healthy subjects receiving single intravenous doses of 13.1 mg |
| [Vogelsang 1984](#5-references) | Healthy subjects receiving single oral doses of 250mg R-verapamil |
| [Blume 1983](#5-references) | Healthy subjects receiving single oral doses of 40 mg |
| [Blume 1990](#5-references) | Healthy subjects receiving single oral doses of 40 mg |
| [John 1992](#5-references) | Healthy subjects receiving single oral doses of 40 mg |
| [Sawicki 2002](#5-references) | Healthy subjects receiving single oral doses of 40 mg |
| [Choi 2008](#5-references) | Healthy subjects receiving single oral doses of 60 mg |
| [Wing 1985](#5-references) | Healthy subjects receiving single oral doses of 80 mg |
| [Maeda 2011](#5-references) | Healthy subjects receiving single oral doses of 80 mg |
| [Ratiopharm 1989](#5-references) | Healthy subjects receiving single oral doses of 80 mg |
| [Boehringer 2018](#5-references) | Healthy subjects receiving single oral doses of 120 mg |
| [Blume 1987](#5-references) | Healthy subjects receiving single oral doses of 120 mg |
| [Johnston 1981](#5-references) | Healthy subjects receiving single oral doses of 120 mg |
| [Mikus 1990](#5-references) | Healthy subjects receiving single oral doses of 160 mg |
| [van Haarst 2009](#5-references) | Healthy subjects receiving multiple oral doses of 180 mg BID |
| [Blume 1994](#5-references) | Healthy subjects receiving single oral doses of 240 mg QD |
| [Joergenson 1988](#5-references) | Healthy subjects receiving multiple oral doses of 120 mg BID |
| [Shand 1981](#5-references)     | Healthy subjects receiving multiple oral doses of 120 mg TID |
| [Karim 1995](#5-references)     | Healthy subjects receiving single oral doses of 240 mg |

## Literature in-vitro kinetics (given)

These are the measured in-vitro values from the literature - the fitted
catalytic rates are NOT given; you must determine them.

| molecule | process | parameter | value | unit |
| --- | --- | --- | --- | --- |
| P-gp | ActiveTransport_InVitro_VesicularAssay_MM | In vitro Vmax/transporter | 0.00057724 | pmol/min/pmol transporter |
| P-gp | ActiveTransport_InVitro_VesicularAssay_MM | Km | 1.01 | µmol/l |
| CYP3A4 | MetabolizationLiverMicrosomes_MM | In vitro Vmax for liver microsomes | 1.27 | nmol/min/mg mic. protein |
| CYP3A4 | MetabolizationLiverMicrosomes_MM | Km | 19.59 | µmol/l |
| CYP3A4 | MetabolizationLiverMicrosomes_MM | In vitro Vmax for liver microsomes | 1.17 | nmol/min/mg mic. protein |
| CYP3A4 | MetabolizationLiverMicrosomes_MM | Km | 35.34 | µmol/l |
