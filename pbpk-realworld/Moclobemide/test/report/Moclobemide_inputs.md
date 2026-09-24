# Building and evaluation of a PBPK model for Moclobemide in adults


> You are handed this report and the clinical data files. Extract the

> physicochemical and in-vitro inputs yourself, decide the model

> structure, and fit the parameters. The reference model's chosen

> methods and fitted values are withheld.


# Background

The presented PBPK model of moclobemide has been developed to be used in a PBPK Drug-Drug-Interactions (DDI) network with moclobemide as a substrate and a moderate inhibitor of CYP2C19.

Moclobemide is a reversible inhibitor of monoamine oxidase A (RIMA), a drug primarily used to treat depression and social anxiety ([Mayersohn 1995](#5-references)). Moclobemide pharmacokinetics is characterized by non-linearity in dose and time. Cmax concentrations decrease with dose for doses above 100 mg. Furthermore, a saturation in clearance at higher doses (200 mg and up) could be seen, as indicated by a longer terminal phase. In addition, multiple doses administration resulted in higher moclobemide concentrations compared to single dose. This could be indicative of the previously reported auto-inhibitory effect ([Nair 1993](#5-references)).

**Absorption**: Moclobemide is highly soluble, and consequently fast and completely absorbed. Absolute bioavailability has been reported to be dependent on the dose, likely due to saturable (first-pass) metabolism ([Mayersohn 1995](#5-references)).

**Distribution**: Moclobemide is moderately bound to plasma proteins and due to its lipophilicity distributes widely in the body (Vss ~ 1.2 L/kg) ([MHRA Label Moclobemide](#5-references)).

**Metabolism**: About 99% of a dose is metabolized mainly via CYP2C19 (C-oxidation, producing metabolite RO12-8095) and FMO3 (N-oxidation, producing metabolite RO12-5637) ([Gram 1995](#5-references), [Hoskins 2001](#5-references), [Mayersohn 1995](#5-references)).

**Excretion**: Less than 1% of a dose is excreted unchanged via the kidneys.

## 2.2 Data<a id="data"></a>

### 2.2.1 In vitro and physico-chemical data

A literature search was performed to collect available information on physico-chemical properties of moclobemide, see [Table 1](#table-1).

| **Parameter**                   | **Unit**          | **Value**    | Source                                  | **Description**                                              |
| :------------------------------ | ----------------- | ------------ | --------------------------------------- | ------------------------------------------------------------ |
| MW<sup>+</sup>                  | g/mol             | 268.74       | [DrugBank DB01171](#5-references)       | Molecular weight                                             |
| pK<sub>a,base</sub><sup>+</sup> |                   | 6.2          | [IPCS INCHEM](#5-references)             | Acidic dissociation constant                                 |
| Solubility (pH)<sup>+</sup>     | mg/mL             | 3<br />(6.8) | [Moclobemide, INCHEM](#5-references)    | Aqueous Solubility                                           |
| logD                            |                   | 1.79         | [Pons 1990](#5-references)              | Distribution coefficient                                     |
| fu<sup>+</sup>                  | %                 | 50           | [MHRA Label Moclobemide](#5-references) | Fraction unbound in plasma                                   |
| Km_FMO (microsomes)             | mmol/L            | 0.77         | [Hoskins 2001](#5-references)           |                                                              |
| Vmax_FMO (microsomes)           | nmol/min/mg prot. | 1.39         | [Hoskins 2001](#5-references)           |                                                              |
| Renal Elimination               | ml/min/kg         | 0.03         | [Schoerlin 1987](#5-references)         | Schoerlin reports 2.6 ml/min/76kg                            |
| Ki_CYP2C19 (free)               | µmol/L            | 203.8        | [Kramer-Nielsen 1996](#5-references)    | The total ki value reported by Kramer was 210 umol/L and corrected with an fu_mic of 0.97 |

**Table 1:**<a name="table-1"></a> Physico-chemical and *in-vitro* metabolization properties of moclobemide extracted from literature. *<sup>+</sup>: Value used in final model*

### 2.2.2 Clinical data

A literature search was performed to collect available clinical data on moclobemide in adults, see [Table 2](#table-2).

| **Source**           | **Route** | **Dose [mg]/**  **Schedule \*** | **Pop.**     | Age [yrs] (mean or range) | Weight [kg] (mean or range) | **Sex** | **N** | **Form.** | **Comment**                       |
| -------------------- | --------- | ------------------------------- | ------------ | ------- | ----- | --------- | --------------------------------- | --------------------------------- | --------------------------------- |
| [Gram 1995](#5-references)<sup>+</sup> | p.o.                  | 300 s.d. / b.i.d.               | HV                | 26                        | -                           | m/f     | 8     | tablet           | EM + PM             |
| [Yu 2001](#5-references)<sup>+</sup> | p.o.                  | 300 s.d.                        | HV-Asian          | -                         | 60.3                        | m       | 8     | tablet           | EM, PM and EM+OMP40 |
| [Wiesel 1985](#5-references)<sup>+</sup> | p.o.                  | 50, 100, 200 s.d.               | HV or patient etc | 26.3                      | 75.8                        | m       | 6     | tablet           |                        |
| [Ignjatovic 2009](#5-references) | p.o.                  | 150 t.i.d.                      | Pat               | -                         | -                           | m/f     | 6     | tablet           |                        |
| [Dingemanse 1998](#5-references) | p.o.                  | 300 s.d.                        | HV                | -                         | -                           | f/m     | 12    | tablet           |                        |
| [Schoerlin 1987](#5-references)<sup>+</sup> | p.o. &  i.v. infusion | 150 t.i.d. /s.d.                | HV                | 27                        | 76                          | m       | 12    | tablet/ solution |                        |
| [Guentert 1990](#5-references)<sup>+</sup> | p.o.                  | 150 t.i.d.                      | HV                | 19-29                     | 59-86                       | m/f     | 14    | tablet           |                        |
| [Raaflaub 1984](#5-references)<sup>+</sup> | p.o. & i.v. infusion  | 50 s.d.                         | HV                | 42                        | 4                           | m       | 6     | tablet/ solution |                        |

**Table 2:**<a name="table-2"></a> Literature sources of clinical concentration data of moclobemide used for model development and validation. *-: respective information was not provided in the literature source; \*:single dose unless otherwise specified; EM: extensive metabolizers; PM: poor metabolizers; <sup>+</sup>: Data used for final parameter identification*
