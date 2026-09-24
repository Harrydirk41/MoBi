# Building and evaluation of a PBPK model for digoxin in adults


> You are handed this report and the clinical data files. Extract the

> physicochemical and in-vitro inputs yourself, decide the model

> structure, and fit the parameters. The reference model's chosen

> methods and fitted values are withheld.


# Background

The presented model building and evaluation report evaluates the performance of a PBPK model for digoxin in adults.

Digoxin is a cardiac glycoside used to treat atrial fibrillation, atrial flutter and heart failure.
Digoxin is transported by P-glycoprotein 1 (P-gp), also known as multidrug resistance protein 1 (MDR1) or ATP-binding cassette sub-family B member 1 (ABCB1) or cluster of differentiation 243 (CD243) poly-glycoprotein. P-gp and mainly excreted unchanged via the kidneys with a small fraction eliminated via biliary excretion and only a very low degree of hepatic metabolism  ([Greiner 1999, Ochs 1978 ](#References)).
Many other substrates of P-gp are metabolized by CYP3A4, setting digoxin apart as an
exception and thereby turning it into a model victim drug of P-gp-mediated DDIs.

Digoxin is reported to have a large volume of distribution due to extensive tissue binding and to be mainly excreted unchanged to urine (50 - 70%) while the remainder of a dose is eliminated by hepatic metabolism and biliary excretion ([Ochs 1978, Bauer 2008](#References))]. The final digoxin model applies target-binding, transport by P-gp in various organs including gut, liver and kidney, an unspecific hepatic metabolic clearance and glomerular filtration, and adequately described the pharmacokinetics of digoxin in adults.

The digoxin model is a whole-body PBPK model, allowing for dynamic translation between individuals. The digoxin report demonstrates the level of confidence in the digoxin PBPK model with the OSP suite with regard to reliable predictions of digoxin PK in adults during model-informed drug development.

## 2.2 Data<a id="data"></a>

### 2.2.1 In vitro / physico-chemical data

A literature search was performed to collect available information on physicochemical properties of digoxin. The obtained information from the literature is summarized in the table below, and is used for model building.

| **Parameter**   | **Unit**    | **Digoxin literature** | **Description**                                  |
| :-------------- | ----------- | ----------------------------------- | ------------------------------------------------ |
| MW              | g/mol       | 780.93 ([Drugbank](#References)) | Molecular weight                                 |
| pKa             |             | none | Acid dissociation constant                   |
| Solubility (pH) | mg/L       | 64.8 (7) ([Drugbank](#References)) | Solubility                                       |
| fu              |             | 70.0, 71.0, 77.7 ([Hinderling 1984, Obach 2008, Neuhoff 2013](#References)) | Fraction unbound                                 |
| ATP1A2 KD | µmol/L | 0.0256 ([Katz 2010](#References)) | Dissociation constant |
| ATP1A2 koff | 1/min | n.a. | Dissociation rate constant |
| P-gp KM | µmol/L | 73.0, 177.0 ([Collett 2004, Troutman 2003](#References)) | Michaelis-Menten constant |
| P-gp kcat       | 1/min    | n.a. | P-gp catalytic rate constant |
| CLhep       | mL/min | n.a. | Hepatic plasma clearance |
| GFR fraction       | 1/min    | 1 | Fraction of filtered drug reaching the urine |
| Formulation       |     | solution | Formulation used in predictions |
| Cell permeabilities       |     | PK-Sim | Permeation across cell membranes |
| Specific intest. perm.       | dm/min   | n.a. | Normalized to surface area |
| Specific organ perm.       | dm/min   | n.a. | Normalized to surface area |


  

### 2.2.2 Clinical data

A literature search was performed to collect available clinical data on digoxin in adults. 

The following publications were found in adults for model building and evaluation:

| Publication                         | Study  description                                           |
| ----------------------------------- | ------------------------------------------------------------ |
| [Becquemont   2001](#References)    | Becquemont, L. et al. Effect  of grapefruit juice on digoxin pharmacokinetics in humans. Clin. Pharmacol.  Ther. 70, 311–6 (2001). |
| [Ding   2004](#References)          | Ding, R. et al. Substantial  pharmacokinetic interaction between digoxin and ritonavir in healthy  volunteers. Clin. Pharmacol. Ther. 76, 73–84 (2004). |
| [Eckermann   2012](#References)     | Eckermann, G., Lahu, G.,  Nassr, N. & Bethke, T.D. Absence of pharmacokinetic interaction between  roflumilast and digoxin in healthy adults. J. Clin. Pharmacol. 52, 251–7  (2012). |
| [Friedrich   2011](#References)     | Friedrich, C. et al.  Evaluation of the pharmacokinetic interaction after multiple oral doses of  linagliptin and digoxin in healthy volunteers. Eur. J. Drug Metab.  Pharmacokinet. 36, 17–24 (2011). |
| [Greiner   1999](#References)       | Greiner, B. et al. The role of  intestinal P-glycoprotein in the interaction of digoxin and rifampin. J.  Clin. Invest. 104, 147–53 (1999). |
| [Gurley   2008b](#References)       | Gurley, B.J., Swain, A.,  Williams, D.K., Barone, G. & Battu, S.K. Gauging the clinical  significance of P-glycoprotein-mediated herb-drug interactions: comparative  effects of St. John’s wort, Echinacea, clarithromycin, and rifampin on  digoxin pharmacokinetics. Mol. Nutr. food Res. 52, 772–9 (2008). |
| [Hayward   1978](#References)       | Hayward, R.P., Greenwood, H.  & Hamer, J. Comparison of digoxin and medigoxin in normal subjects. Br.  J. Clin. Pharmacol. 6, 81–6 (1978). |
| [Jalava   1997 ](#References)       | Jalava, K.M., Partanen, J.  & Neuvonen, P.J. Itraconazole decreases renal clearance of digoxin. Ther.  Drug Monit. 19, 609–13 (1997). |
| [Johne   1999](#References)         | Johne, A. et al.  Pharmacokinetic interaction of digoxin with an herbal extract from St John’s  wort (Hypericum perforatum). Clin. Pharmacol. Ther. 66, 338–45 (1999). |
| [Kirby   2012](#References)         | Kirby B.J., Collier A.C., Kharasch E.D., Whittington D., Thummel K.E., Unadkat J.D. Complex drug interactions of the HIV protease inhibitors 3: effect of simultaneous or staggered dosing of digoxin and ritonavir, nelfinavir, rifampin, or bupropion. Drug Metab Dispos. 2012 Mar;40(3):610-6. |
| [Kirch   1986](#References)         | Kirch, W., Hutt, H.J.,  Dylewicz, P., Gräf, K.J. & Ohnhaus, E.E. Dose-dependence of the  nifedipine-digoxin interaction? Clin. Pharmacol. Ther. 39, 35–9 (1986). |
| [Koup   1975](#References)          | Koup, J.R., Greenblatt, D.J.,  Jusko, W.J., Smith, T.W. & Koch-Weser, J. Pharmacokinetics of digoxin in  normal subjects after intravenous bolus and infusion doses. J. Pharmacokinet.  Biopharm. 3, 181–92 (1975). |
| [Kramer   1979](#References)        | Kramer, W.G. et al.  Pharmacokinetics of digoxin: relationship between response intensity and  predicted compartmental drug levels in man. J. Pharmacokinet. Biopharm. 7,  47–61 (1979). |
| [Lalonde   1985](#References)       | Lalonde, R.L., Deshpande, R.,  Hamilton, P.P., McLean, W.M. & Greenway, D.C. Acceleration of digoxin  clearance by activated charcoal. Clin. Pharmacol. Ther. 37, 367–71 (1985). |
| [Larsen   2007](#References)        | Larsen, U.L. et al. Human  intestinal P-glycoprotein activity estimated by the model substrate digoxin.  Scand. J. Clin. Lab. Invest. 67, 123–34 (2007). |
| [Martin   1997](#References)        | Martin, D.E. et al. Lack of  effect of eprosartan on the single dose pharmacokinetics of orally  administered digoxin in healthy male volunteers. Br. J. Clin. Pharmacol. 43,  661–4 (1997). |
| [Ochs   1975](#References)          | Ochs, H., Bodem, G., Schäfer, P.K., Kodrat, G., Dengler, H.J. Absorption of digoxin from the distal parts of the intestine in man. Eur J Clin Pharmacol. 1975 Dec 19;9(2-3):95-7. |
| [Ochs   1978](#References)          | Ochs, H., Greenblatt, D.J.,  Bodem, G. & Harmatz, J.S. Dose-independent pharmacokinetics of digoxin in  humans. Am. Heart J. 96, 507–11 (1978). |
| [Oosterhuis   1991](#References)    | Oosterhuis, B., Jonkman, J.H.,  Andersson, T., Zuiderwijk, P.B. & Jedema, J.N. Minor effect of multiple  dose omeprazole on the pharmacokinetics of digoxin after a single oral dose.  Br. J. Clin. Pharmacol. 32, 569–72 (1991). |
| [Qiu   2010](#References)           | Qiu, R. et al. Lack of a  pharmacokinetic interaction between dimebon (latrepirdine) and digoxin in  healthy subjects. Am. Soc. Clin. Pharmacol. Ther. Meet. Atlanta, GA, USA  (2010). |
| [Ragueneau   1999](#References)     | Ragueneau, I. et al.  Pharmacokinetic and pharmacodynamic drug interactions between digoxin and  macrogol 4000, a laxative polymer, in healthy volunteers. Br. J. Clin.  Pharmacol. 48, 453–6 (1999). |
| [Rengelshausen   2003](#References) | Rengelshausen, J. et al.  Contribution of increased oral bioavailability and reduced nonglomerular  renal clearance of digoxin to the digoxin-clarithromycin interaction. Br. J.  Clin. Pharmacol. 56, 32–8 (2003). |
| [Rodin   1988](#References)         | Rodin, S.M., Johnson, B.F.,  Wilson, J., Ritchie, P. & Johnson, J. Comparative effects of verapamil  and isradipine on steady-state digoxin kinetics. Clin. Pharmacol. Ther. 43,  668–72 (1988). |
| [Steiness   1982](#References)      | Steiness, E., Waldorff, S.  & Hansen, P.B. Renal digoxin clearance: dependence on plasma digoxin and  diuresis. Eur. J. Clin. Pharmacol. 23, 151–4 (1982). |
| [Tayrouz   2003](#References)       | Tayrouz, Y. et al.  Pharmacokinetic and pharmaceutic interaction between digoxin and Cremophor  RH40. Clin. Pharmacol. Ther. 73, 397–405 (2003). |
| [Tsutsumi   2002](#References)      | Tsutsumi, K. et al. The effect  of erythromycin and clarithromycin on the pharmacokinetics of intravenous  digoxin in healthy volunteers. J. Clin. Pharmacol. 42, 1159–64 (2002). |
| [Vaidyanathan   2008](#References)  | Vaidyanathan, S. et al.  Pharmacokinetics of the oral direct renin inhibitor aliskiren in combination  with digoxin, atorvastatin, and ketoconazole in healthy subjects: the role of  P-glycoprotein in the disposition of aliskiren. J. Clin. Pharmacol. 48, 1323–38  (2008). |
| [Verstuyft   2003](#References)     | Verstuyft, C. et al.  Dipyridamole enhances digoxin bioavailability via P-glycoprotein inhibition.  Clin. Pharmacol. Ther. 73, 51–60 (2003). |
| [Wagner   1981](#References)        | Wagner, J.G., Popat, K.D.,  Das, S.K., Sakmar, E. & Movahhed, H. Evidence of nonlinearity in digoxin  pharmacokinetics. J. Pharmacokinet. Biopharm. 9, 147–66 (1981). |
| [Westphal   2000](#References)      | Westphal, K. et al. Oral  bioavailability of digoxin is enhanced by talinolol: evidence for involvement  of intestinal P-glycoprotein. Clin. Pharmacol. Ther. 68, 6–12 (2000). |

## Literature in-vitro kinetics (given)

These are the measured in-vitro values from the literature - the fitted
catalytic rates are NOT given; you must determine them.

| molecule | process | parameter | value | unit |
| --- | --- | --- | --- | --- |
| P-gp | ActiveTransportSpecific_MM | Vmax | 8.67 | µmol/l/min |
| P-gp | ActiveTransportSpecific_MM | Km | 177.0 | µmol/l |
