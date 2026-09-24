# Building and evaluation of a PBPK model for erythromycin in healthy adults


> You are handed this report and the clinical data files. Extract the

> physicochemical and in-vitro inputs yourself, decide the model

> structure, and fit the parameters. The reference model's chosen

> methods and fitted values are withheld.


# Background

The presented model building and evaluation report evaluates the performance of a PBPK model for erythromycin in healthy adults.

Erythromycin a macrolide antibiotic used for the treatment of a number of bacterial infections, including respiratory tract infections, skin infections, chlamydia infections, and others. Erythromycin is available in the form of various salts and formulations, for example as: 

- erythromycin lactobionate for intravenous injection
- erythromycin base in enteric-coated capsules or tablets for oral administration
- erythromycin stearate in filmcoated tablets for oral administration
- erythromycin ethylsuccinate in suspension or in filmcoated tablets for oral administration

In its free form as base, erythromycin is easily hydrolyzed in acidic aqueous solution ([Mordi 2000](#5-references)). Therefore, orally administered erythromycin is given in the form of enteric-coated formulations or as more acid-stable salts or esters of erythromycin (e.g. erythromycin ethylsuccinate). Once in the small intestine, erythromycin is rapidly absorbed displaying a highly variable bioavailability ([Chun 1977](#5-references), [Mather 1981](#5-references)). Erythromycin diffuses in most tissues and accumulates in leukocytes and phagocytes ([Miller 1984](#5-references), [Carlier 1987](#5-references)). About 70% of erythromycin is bound to plasma proteins ([Barre 1987](#5-references)). Erythromycin has been shown to be a substrate for various transporters including P-gp and OATP1B1. The latter has been shown to critically affect erythromycin disposition ([Lancaster 2012](#5-references)). Erythromycin is extensively metabolized through N-demethylation catalyzed by CYP3A. Metabolism via CYP4F11 has also been suggested ([Kalsotra 2004](#5-references)). Biliary excretion also appears to play an important role in erythromycin clearance ([Acocella 1968](#5-references), [Chelvan 1979](#5-references)), but its contribution to total elimination remains unknown. The dose fraction excreted unchanged in urine is minimal and highly variable; reported fractions after IV administration range from 0.018 ± 0.005 to 0.171 ± 0.11 (mean ± SD) ([Pasic 1987](#5-references), [Austin 1980](#5-references)). There is abundant evidence from *in vitro* studies that erythromycin irreversibly inhibits CYP3A (e.g. [Larrey 1983](#5-references)) and the FDA lists erythromycin as moderate index inhibitor for CYP3A. Findings from *in vivo* studies investigating the dose linearity of erythromycin pharmacokinetics are not fully conclusive, but some studies observed a slight dose dependency ([Austin 1980](#5-references), [Josefsson 1982](#5-references)).

## 2.2 Data<a id="data"></a>

### 2.2.1	In vitro / physicochemical data

A literature search was carried out to collect available information on physicochemical properties of erythromycin. The obtained information from the literature is summarized in the table below and is used for model building.

| **Parameter**                                | **Unit** | **Literature**                                               | **Description**                                  |
| :------------------------------------------- | -------- | ------------------------------------------------------------ | ------------------------------------------------ |
| Molecular weight                             | g/mol    | 733.9 ([drugbank.ca](#5-references))                        | Molecular weight                                 |
| pK<sub>a</sub> (basic)                       |          | 8.8 ([Lien 1974](#5-references)); 8.88 ([McFarland 1997](#5-references)) | Acid dissociation constant                       |
| f<sub>u</sub>                                |          | 0.27 ± 0.05<sup>a</sup> ([Sun 2010](#5-references)); 0.28 ± 0.04<sup>a</sup> ([Iliopoulou 1982](#5-references)); 0.305 ± 0.028<sup>a</sup> ([Barre 1987](#5-references)); 0.326 ([Xu 2009](#5-references)) | Fraction unbound in human plasma                 |
| Water solubility (erythromycin lactobionate) | mg/L     | 200000 ([U.S. Patent 2,761,859](#5-references))             | Solubility of erythromycin lactobionate in water |
| Water solubility (erythromycin stearate)     | mg/L     | 182 ([Jones 1969](#5-references))                           | Solubility of erythromycin stearate in water     |
| Water solubility (erythromycin base)         | mg/L     | 2100 ([Manna 1998](#5-references))                          | Solubility of free erythromycin as base in water |

 <sup>a</sup> mean ± SD

### 2.2.2	In vitro data on mechanism-based inhibition of CYP3A

A literature search was carried out to collect quantitative information the kinetic parameters describing the mechanism-based inhibition of CYP3A by erythromycin, K<sub>I</sub> and k<sub>inact</sub>. Ample data was collected from the literature and is summarized in the table below.

| Victim   compound | In   vitro system | Enzyme | K<sub>I</sub>   [µM] | k<sub>inact</sub>   [min<sup>-1</sup>] | Concentration   range of erythromycin [µM] | Concentration   of victim compound [µM] | Reference                           |
| ----------------- | ----------------- | ------ | -------------------- | -------------------------------------- | ------------------------------------------ | --------------------------------------- | ----------------------------------- |
| testosterone      | rec cDNA CYP3A4   | CYP3A4 | 1.04                 | 0.0293                                 | 0 - 20                                     | 100                                     | [Akiyoshi 2013](#5-references)     |
| testosterone      | rec cDNA CYP3A4   | CYP3A4 | 1.21                 | 0.0164                                 | 0 - 20                                     | 100                                     | [Akiyoshi 2013](#5-references)     |
| testosterone      | rec cDNA CYP3A4   | CYP3A4 | 0.415                | 0.0159                                 | 0 - 20                                     | 100                                     | [Akiyoshi 2013](#5-references)     |
| testosterone      | rec cDNA CYP3A4   | CYP3A4 | 2.24                 | 0.0293                                 | 0 - 20                                     | 100                                     | [Akiyoshi 2013](#5-references)     |
| testosterone      | rec cDNA CYP3A4   | CYP3A4 | 0.753                | 0.0248                                 | 0 - 20                                     | 100                                     | [Akiyoshi 2013](#5-references)     |
| midazolam         | rec cDNA CYP3A4   | CYP3A  | 9.5                  | 0.16                                   | 50-fold range                              | 10                                      | [Atkinson 2005](#5-references)     |
| midazolam         | rec cDNA CYP3A4   | CYP3A  | 8.82                 | 0.12                                   | 50-fold range                              | 10                                      | [Atkinson 2005](#5-references)     |
| midazolam         | rec cDNA CYP3A4   | CYP3A  | 8.3                  | 0.09                                   | 50-fold range                              | 10                                      | [Atkinson 2005](#5-references)     |
| testosterone      | rec cDNA CYP3A4   | CYP3A  | 9.5                  | 0.06                                   | 4 - 200                                    | 16.7                                    | [Atkinson 2005](#5-references)     |
| triazolam         | HLM               | CYP3A  | 5.4                  | 0.069                                  |                                            | 400                                     | [Aueviriyavit 2010](#5-references) |
| testosterone      | rec cDNA CYP3A4   | CYP3A4 | 5.3                  | 0.12                                   | 5 - 20                                     | 200                                     | [Chan 2000](#5-references)         |
| midazolam         | cryopreserved HEP | CYP3A  | 60                   | 0.081                                  | 3 - 300                                    | 30                                      | [Chen 2011](#5-references)         |
| midazolam         | cryopreserved HEP | CYP3A  | 67.9                 | 0.079                                  | 3 - 300                                    | 30                                      | [Chen 2011](#5-references)         |
| midazolam         | rec cDNA CYP3A4   | CYP3A4 | 0.762                | 0.0648                                 | 0.1 - 30                                   | 4                                       | [Ishikawa 2017](#5-references)     |
| testosterone      | rec cDNA CYP3A4   | CYP3A4 | 1.00                 | 0.0604                                 | 0.1 - 30                                   | 150                                     | [Ishikawa 2017](#5-references)     |
| nifedipine        | rec cDNA CYP3A4   | CYP3A4 | 0.794                | 0.0766                                 | 0.1 - 30                                   | 6                                       | [Ishikawa 2017](#5-references)     |
| triazolam         | HLM               | CYP3A4 | 15.9                 | 0.062                                  | 3 - 100                                    | 300                                     | [Kanamitsu 2000](#5-references)    |
| triazolam         | HLM               | CYP3A4 | 17.4                 | 0.055                                  | 3 - 100                                    | 300                                     | [Kanamitsu 2000](#5-references)    |
| triazolam         | rec cDNA CYP3A4   | CYP3A4 | 19.1                 | 0.173                                  | 3 - 100                                    | 300                                     | [Kanamitsu 2000](#5-references)    |
| triazolam         | rec cDNA CYP3A4   | CYP3A4 | 18.9                 | 0.097                                  | 3 - 100                                    | 300                                     | [Kanamitsu 2000](#5-references)    |
| testosterone      | HLM               | CYP3A4 | 29.4                 | 0.0271                                 |                                            | 250                                     | [Kosaka 2017](#5-references)       |
| testosterone      | HLM               | CYP3A4 | 30                   | 0.040                                  | 0.3 - 300                                  |                                         | [Kozakai 2013](#5-references)      |
| midazolam         | HLM               | CYP3A4 | 12                   | 0.035                                  | 0.3 - 300                                  |                                         | [Kozakai 2013](#5-references)      |
| midazolam         | HLM               | CYP3A4 | 20                   | 0.033                                  | 0.3 - 300                                  |                                         | [Kozakai 2013](#5-references)      |
| midazolam         | cryopreserved HEP | CYP3A  | 25.15                | 0.08                                   | 0.13 - 100                                 | 30                                      | [Mao 2011](#5-references)          |
| midazolam         | HLM               | CYP3A  | 10.8                 | 0.032                                  |                                            |                                         | [Mao 2016](#5-references)          |
| midazolam         | cryopreserved HEP | CYP3A  | 30.7                 | 0.05                                   | 0 - 300                                    | 20                                      | [Mao 2016](#5-references)          |
| midazolam         | cryopreserved HEP | CYP3A  | 59.2                 | 0.062                                  | 0 - 100                                    | 5                                       | [Mao 2016](#5-references)          |
| midazolam         | cryopreserved HEP | CYP3A  | 80.3                 | 0.052                                  | 0 - 100                                    | 20                                      | [Mao 2016](#5-references)          |
| midazolam         | rec cDNA CYP3A4   | CYP3A4 | 7.47                 | 0.042                                  | 2 - 50                                     | 8                                       | [McConn 2004](#5-references)       |
| midazolam         | HLM               | CYP3A4 | 10.9                 | 0.046                                  | 2 - 100                                    | 8                                       | [McConn 2004](#5-references)       |
| midazolam         | primary HEP       | CYP3A  | 11                   | 0.07                                   | 0.1 - 10                                   | 3                                       | [McGinnity 2006](#5-references)    |
| midazolam         | HLM               | CYP3A4 | 10                   | 0.036                                  | 0 - 25                                     |                                         | [Obach 2007](#5-references)        |
| testosterone      | HLM               | CYP3A4 | 9.8                  | 0.039                                  | 0 - 25                                     | 500                                     | [Obach 2007](#5-references)        |
| testosterone      | rec cDNA CYP3A4   | CYP3A4 | 0.92                 | 0.058                                  | 5 - 100                                    | 250                                     | [Polasek 2006](#5-references)      |
| testosterone      | HLM               | CYP3A4 | 12.8                 | 0.037                                  | 5 - 100                                    | 250                                     | [Polasek 2006](#5-references)      |
| midazolam         | rec cDNA CYP3A4   | CYP3A4 | 5.1                  | 0.30                                   | 0.5 - 50                                   | 100                                     | [Ring 2005](#5-references)         |
| testosterone      | rec cDNA CYP3A4   | CYP3A4 | 0.92                 | 0.058                                  |                                            |                                         | [Teng 2010](#5-references)         |
| testosterone      | HLM               | CYP3A4 | 4.579                | 0.0115                                 |                                            |                                         | [Teng 2010](#5-references)         |
| domperidone       | HLM               | CYP3A4 | 18.4                 | 0.022                                  | 2.5 - 200                                  | 500                                     | [Ung 2009](#5-references)          |
| domperidone       | rec cDNA CYP3A4   | CYP3A4 | 4.1                  | 0.026                                  | 2.5 - 200                                  | 500                                     | [Ung 2009](#5-references)          |
| midazolam         | HLM               | CYP3A4 | 12.1                 | 0.0215                                 | 0 - 100                                    | 25                                      | [Watanabe 2007](#5-references)     |
| nifedipine        | HLM               | CYP3A4 | 11.3                 | 0.0295                                 | 0 - 100                                    | 50                                      | [Watanabe 2007](#5-references)     |
| testosterone      | HLM               | CYP3A4 | 10.9                 | 0.0352                                 | 0 - 100                                    | 200                                     | [Watanabe 2007](#5-references)     |
| midazolam         | HLM               | CYP3A  | 1.48                 | 0.017                                  | 0.5 - 500                                  | 20                                      | [Xu 2009](#5-references)           |
| midazolam         | cryopreserved HEP | CYP3A  | 20.0                 | 0.016                                  | 0.5 - 500                                  | 20                                      | [Xu 2009](#5-references)           |
| midazolam         | cryopreserved HEP | CYP3A  | 109                  | 0.055                                  | 0.5 - 500                                  | 20                                      | [Xu 2009](#5-references)           |
| midazolam         | HLM               | CYP3A  | 81.8                 | 0.0665                                 | 20 - 400                                   | 10                                      | [Yamano 2001](#5-references)       |
| midazolam         | HLM               | CYP3A  | 15.7                 | 0.1                                    |                                            |                                         | [Zhang 2006](#5-references)        |
| testosterone      | rec cDNA CYP3A4   | CYP3A4 | 5                    | 0.34                                   | 1 - 50                                     | 200                                     | [Zhang 2009](#5-references)        |
| testosterone      | HLM               | CYP3A4 | 15.7                 | 0.09                                   | 1 - 50                                     | 200                                     | [Zhang 2009](#5-references)        |
| midazolam         | HLM               | CYP3A4 | 26.5                 | 0.041                                  | 2.5 - 50                                   | 10                                      | [Zimmerlin 2011](#5-references)    |

*Note:* Abbreviations: HEP: human hepatocytes; HLM: human liver microsomes; rec cDNA CYP3A4: human recombinant c-DNA CYP3A4 enzymes (e.g. supersomes, baculovirus-insect cell system, E. coli transfected cells)

The data listed in the Table above can be statistically summarized as follows:

| Parameter [unit]                     | Min    | Q1     | Geometric mean | Median | Arithmetic mean | Q3     | Max   |
| ------------------------------------ | ------ | ------ | -------------- | ------ | --------------- | ------ | ----- |
| K<sub>I</sub> [µM]                   | 0.420  | 4.89   | 8.71           | 10.9   | 18.4            | 19.3   | 109   |
| k<sub>inact</sub> [min<sup>-1</sup>] | 0.0115 | 0.0314 | 0.0504         | 0.0535 | 0.0664          | 0.0772 | 0.340 |

### 2.2.3 Clinical data

A literature search was carried out to collect available PK data on erythromycin in healthy adults. The following data from the publications listed below were used for model building and evaluation:

| Publication                       | Study description                                            |
| :-------------------------------- | :----------------------------------------------------------- |
| [Austin 1980](#5-references)     | IV administration of 125 mg, 250 mg, 500 mg, and 900 mg as erythromycin lactobionate; single dose |
| [Barre 1987](#5-references)      | IV administration of 500 mg as erythromycin lactobionate; single dose |
| [Berend 1979](#5-references)     | PO administration of 500 mg erythromycin stearate; multiple dose |
| [Birkett 1990](#5-references)    | PO administration of 250 mg as enteric-coated capsules containing pellets of erythromycin base; single and multiple dose |
| [Brannan 1995](#5-references)    | PO administration of 500 mg erythromycin stearate; multiple dose |
| [DiSanto 1981](#5-references)    | PO administration of 500 mg as unprotected tablets containing erythromycin base, as film-coated tablets containing erythromycin base, as enteric-coated tablets containing erythromycin base and as film-coated tablets containing erythromycin stearate; single and multiple dose |
| [Henry 1980](#5-references)      | PO administration of 500 mg erythromycin stearate; single dose |
| [Huppertz 2011](#5-references)   | IV administration of 1000 mg as erythromycin lactobionate; single dose |
| [Iliopoulou 1982](#5-references) | PO administration of 500 mg as film-coated tablet containing erythromycin stearate; multiple dose |
| [Josefsson 1982](#5-references)  | PO administration of 500 mg as film-coated tablet containing erythromycin stearate and 250 mg, 500 mg, 1000 mg as enteric coated capsules containing pellets of erythromycin base; single dose |
| [Kavi 1988](#5-references)       | PO administration of 500 mg as film-coated tablet containing erythromycin stearate; single dose |
| [Kivistö 1997](#5-references)    | PO administration of 500 mg as enteric-coated capsules containing erythromycin base; multiple dose |
| [Kroboth 1982](#5-references)    | PO administration of 500 mg as enteric-coated tablets containing erythromycin base; single dose |
| [Malmborg 1978](#5-references)   | PO administration of 500 mg as film-coated tablet containing erythromycin stearate; multiple dose |
| [Miglioli 1990](#5-references)   | PO administration of 1000 mg erythromycin stearate; multiple dose |
| [Olkkola 1993](#5-references)    | Midazolam-erythromycin interaction study; PO administration of 500 mg as enteric-coated tablet containing erythromycin base; multiple dose |
| [Parsons 1977](#5-references)    | PO administration of 500 mg erythromycin stearate; single dose |
| [Pasic 1987](#5-references)      | IV administration of 300 mg as erythromycin lactobionate; single dose |
| [Posti 1983](#5-references)      | PO administration of 500 mg as film-coated tablets containing erythromycin stearate, as enteric-coated tablets containing erythromycin base, and as enteric-coated tablets containing erythromycin stearate; single dose |
| [Schreiner 1984](#5-references)  | PO administration of 500 mg as film-coated tablets containing erythromycin stearate and as enteric-coated capsule containing pellets of erythromycin base; single dose |
| [Shanson 1984](#5-references)    | PO administration of 1500 mg erythromycin stearate; single dose |
| [Simon 1980](#5-references)      | IV administration of 500 mg erythromycin lactobionate and PO administration of 500 mg erythromycin stearate; single dose |
| [Sun 2010](#5-references)        | IV administration of 125 mg and PO administration of 250 mg (salt and formulation type not specified); single dose |
| [Yakatan 1979](#5-references)    | PO administration of 250 mg as film-coated tablet containing erythromycin stearate; single dose |
| [Yakatan 1980](#5-references)    | PO administration of 250 mg as film-coated tablet containing erythromycin stearate and as enteric-coated tablet containing erythromycin base; single and multiple dose |

## Literature in-vitro kinetics (given)

These are the measured in-vitro values from the literature - the fitted
catalytic rates are NOT given; you must determine them.

| molecule | process | parameter | value | unit |
| --- | --- | --- | --- | --- |
| CYP3A4 | MetabolizationLiverMicrosomes_MM | In vitro Vmax for liver microsomes | 918.33333 | pmol/min/mg mic. protein |
| CYP3A4 | MetabolizationLiverMicrosomes_MM | Km | 70.0 | µM |
| OATP1B1 | ActiveTransportSpecific_MM | Vmax | 11.66 | pmol/ml/min |
