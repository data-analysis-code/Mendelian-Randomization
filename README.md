# Mendelian-Randomization

We performed univariable, multivariable, and two-step two-sample MR mediation analyses to investigate whether genetically predicted educational attainment was causally associated with the risk of EC, BE, and GERD, and to assess the proportion mediated by modifiable risk factors in the above associations. The raw data used to support the findings of this study is included within Supplementary Raw Data.

Example: two-step two-sample MR mediation analysis with education as exposure, BMI and smoking as mediators, and EC as outcome.Statistical analysis codes are as follows.

Total Effect βc: 01the univariable MR analysis of the causal effect of EduYears on EC

Direct Effect βa: 02the univariable MR analysis of the causal effect of EduYears on BMI

Direct Effect βb1: 03the multivariable MR analysis of the causal effect of BMI on EC after adjusting for EduYears

Direct Effect βb2: 04the multivariable MR analysis of the causal effect of Smoking on EC after adjusting for EduYears

Direct Effect βb1+βb2: 05the multivariable MR analysis of the causal effect of BMI and Smoking on EC after adjusting for EduYears

Citation: When using these codes in publications, please remember to cite this paper: Zhang, et al. Association between Education Attainment and Esophageal Cancer, Barrett´s Esophagus, and Gastroesophageal Reflux Disease in European Ancestry Individuals and the Mediating Effect of Modifiable Risk Factors: A Mendelian Randomization Study
