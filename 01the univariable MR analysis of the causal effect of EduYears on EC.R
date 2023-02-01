
library(devtools)
library(tidyr)
library(MRInstruments)
library(MendelianRandomization) 
library(TwoSampleMR) 
#install.packages("simex")
library(simex) 
#install_github("rondolab/MR-PRESSO")
library(MRPRESSO)

#1   Read the final mydata file

mydata <-read.table('01mydata_edu_cancer.txt',header=T,sep="\t",check.names = F,stringsAsFactors = F)


#2   Main MR analysis
heterogeneity <- mr_heterogeneity(mydata)
heterogeneity
#2.1 If I2<25% or Q_P>0.05,fixed IVW Inverse variance weighted (fixed effects) was used
mr_results<-mr(mydata,method_list=c('mr_ivw_fe')) 
mr_results
generate_odds_ratios(mr_results)
#2.2 若I2>25% and Q_P<0.05,random IVW Inverse variance weighted (multiplicative random effects) was used
mr_results<-mr(mydata,method_list=c('mr_ivw_mre')) 
mr_results
generate_odds_ratios(mr_results)

#3   Supplementary MR analysis

mr_results1 <-mr(mydata,method_list=c('mr_two_sample_ml','mr_egger_regression',
                                      "mr_simple_median","mr_weighted_median","mr_penalised_weighted_median",
                                      "mr_simple_mode","mr_weighted_mode")) 

mr_results1
#4   Pleiotropic analysis

#4.1  MR-Egger intercept test
mr_pleiotropy_test(mydata)

mr_egger(mr_input(bx = mydata$beta.exposure, bxse = mydata$se.exposure, by = mydata$beta.outcome, byse = mydata$se.outcome)) 

#4.2  MR pleiotropy residual sum and outlier (MR-PRESSO) test
# Run MR-PRESSO global method
mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", 
          OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = mydata, NbDistribution = 10000,  SignifThreshold = 0.05,seed = 123456) 
# Run MR-PRESSO on a multi-variable MR (MMR) model specifying several exposures
#mr_presso(BetaOutcome = "beta.outcome", BetaExposure = c("E1_effect", "E2_effect"), SdOutcome = "Y_se", SdExposure = c("E1_se", "E2_se"), OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = SummaryStats, NbDistribution = 1000,  SignifThreshold = 0.05)


#5  Weak tool variable validation

#5.1    IVW radial: mr_ivw_radial
mr(mydata,method_list=c("mr_ivw_radial")) 

#5.2    Egger-SIMEXS
mr_egger(mr_input(bx = mydata$beta.exposure, bxse = mydata$se.exposure, by = mydata$beta.outcome, byse = mydata$se.outcome))
BetaXG <- mydata$beta.exposure
BetaYG <- mydata$beta.outcome
seBetaXG <- mydata$se.exposure
seBetaYG <- mydata$se.outcome
BYG <- BetaYG*sign(BetaXG)# Pre-processing steps to ensure all gene--exposure estimates are positive
BXG <- abs(BetaXG)         
# MR-Egger regression (unweighted)
Fit1 <- lm(BYG~BXG,x=TRUE,y=TRUE)
mod.sim1 <- simex(Fit1,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
summary(mod.sim1)
# MR-Egger regression (weighted)
Fit2 <- lm(BYG~BXG,weights=1/seBetaYG^2,x=TRUE,y=TRUE) #Fit2其实就是mr-egger
# Simulation extrapolation 
mod.sim2 <- simex(Fit2,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
summary(mod.sim2)

Fit3 <- lm(BYG~-1+BXG,x=TRUE,y=TRUE) #Fit3其实就是IVW不加权
# Simulation extrapolation 
mod.sim3 <- simex(Fit3,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
summary(mod.sim3)

Fit4 <- lm(BYG~-1+BXG,weights=1/seBetaYG^2,x=TRUE,y=TRUE) #Fit3其实就是IVW加权
# Simulation extrapolation 
mod.sim4 <- simex(Fit4,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
summary(mod.sim4)


#7.Visualization of MR results
#（1）	scatter plot；

plot1 <- mr_scatter_plot(mr_results1, mydata) 
plot1

#（2）	single snp forest plot；
res_single <- mr_singlesnp(
  mydata,
  parameters = default_parameters(),
  single_method = "mr_wald_ratio",
  all_method = c("mr_ivw",'mr_two_sample_ml',"mr_weighted_median")
) 
plot2 <- mr_forest_plot(res_single)
plot2

#（3） Leave-one-out sensitivity test
single <- mr_leaveoneout(mydata)
plot3<-mr_leaveoneout_plot(single)
plot3

#（4）funnel plot
plot4 <-mr_funnel_plot(res_single)
plot4

