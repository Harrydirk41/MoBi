# Building and evaluation of a PBPK model for S-Mephenytoin in adults


> You are handed this report and the clinical data files. Extract the

> physicochemical and in-vitro inputs yourself, decide the model

> structure, and fit the parameters. The reference model's chosen

> methods and fitted values are withheld.


# Background

The presented PBPK model of S-mephenytoin has been developed to be used in a PBPK Drug-Drug-Interactions (DDI) network with S-mephenytoin as a substrate of CYP2C19.

Mephenytoin is a hydantoin-derivative anticonvulsant used to control various partial seizures and was first used in the 1940s ([Troupin 1979](#5-references)).

Only limited clinical PK and ADME data are available. Mephenytoin is soluble and rapidly absorbed with a Tmax of 1 hour. The mean half-life in human is 6.8 hours. No hints for dose non-linearity could be found in literature.

Mephenytoin is the mixture of the two enantiomers S- and R-Mephenytoin. S-Mephenytoin is mainly metabolized via CYP2C19. Only a very minor part is metabolized by CYP2C9. The R-enantiomer is not metabolized by CYP2C19. The clearance of S-Mephenytoin in CYP2C19 EM is 40 to 100-fold higher than in PM.

## 2.2 Data<a id="data"></a>

### 2.2.1 In vitro and physico-chemical data

A literature search was performed to collect available information on physico-chemical properties of S-mephenytoin ([Table 1](#table-1)).

| **Parameter**                   | **Unit** | **Value**     | Source                            | **Description**              |
| :------------------------------ | -------- | ------------- | --------------------------------- | ---------------------------- |
| MW<sup>+</sup>                  | g/mol    | 218.52        | [DrugBank DB00532](#5-references) | Molecular weight.            |
| pK<sub>a,acid</sub><sup>+</sup> |          | 8.51          | [DrugBank DB00532](#5-references) | Acidic dissociation constant |
| Solubility (pH)<sup>+</sup>     | mg/mL    | 1.27<br />(7) | [DrugBank DB00532](#5-references) | Aqueous Solubility           |
| fu<sup>+</sup>                  | %        | 70.2          | [Steere 2015](#5-references)      | Fraction unbound in plasma   |

**Table 1:**<a name="table-1"></a> Physico-chemical and *in-vitro* metabolization properties of S-mephenytoin extracted from literature. *<sup>+</sup>: Value used in final model*

### 2.2.2 Clinical data

A literature search was performed to collect available clinical data on S-mephenytoin. Data used for model development and validation are listed in [Table 2](#table-2).

| **Source**           | **Route** | **Dose [mg]/**  **Schedule** | **Pop.**     | **Age [yrs] (mean) /range** | **Weight [kg] (mean) /range** | **Sex** | **N** | **Form.** | **Comment**                       |
| -------------------- | --------- | ------------------------------- | ------------ | ------- | ----- | --------- | --------------------------------- | --------------------------------- | --------------------------------- |
| [Adedoyin 1998](#5-references) | p.o.      | 100 mg s.d.                   | HV, all EM               | 54.7 / 32-73                | -                             | m/f     | 8     | tablet    | 50 g S-mephenytoin simulated                           |
| [Jacqz 1986](#5-references)    | p.o.      | 100 mg s.d.                   | HV, 6 EM, 1  IM and 1 PM | 25-76                       | -                             | m/f     | 8     | tablet    | 50 g S-mephenytoin simulated                            |
| [Yao 2003](#5-references)      | p.o.      | 100 mg s.d.                   | HV                       | 23-49                       | -                             | m/f     | 12    | -         | S-mephenytoin, with  and without Fluvoxamine MD of 37.5, 62.5 and 87.5 mg/day |
| [Iga 2016](#5-references)      | p.o.      | 100 mg s.d.                   | -                        | -                           | -                             | -       | -     | -         | 50 g S-mephenytoin simulated |
|                                |           |                               |                          |                             |                               |         |       |           |                                                              |
| [Wedlund 1985](#5-references)  | p.o.      | 100 mg s.d.                  | HV                       | 21-7                        | 54-108                        | male    | 8     | tablet    | 50 g S-mephenytoin simulated                    |

**Table 2:**<a name="table-2"></a> Literature sources of clinical concentration data of S-mephenytoin used for model development and validation. *s.d.: single dose*
