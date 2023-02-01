library(devtools)
library(MRInstruments)
library(TwoSampleMR) #加载R包
#install_github("WSpiller/MVMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
library(MVMR)
# install_github("WSpiller/MRPracticals",build_opts = c("--no-resave-data", "--no-manual"),build_vignettes = TRUE)
library(MRPracticals)
vignette("MVMR Tutorial")


library(MVMR)
rawdat_mvmr <- read.table('03edu_BMI_cancer.txt',header=T,sep="\t",check.names = F,stringsAsFactors = F)

F.data<-format_mvmr(BXGs=rawdat_mvmr[,c(4,6)],
                    BYG=rawdat_mvmr[,8],
                    seBXGs=rawdat_mvmr[,c(5,7)],
                    seBYG=rawdat_mvmr[,9],
                    RSID=rawdat_mvmr[,1])
head(F.data)

#  Step 3: Test for weak instruments
# 方法1 conventional F-statistic
sres<-strength_mvmr(r_input=F.data,gencov=0)
#  方法2 conditional F-statistics
mvmrcovmatrix<-matrix(c(1,-0.05,-0.05,1), nrow = 2, ncol = 2)
Xcovmat<-phenocov_mvmr(mvmrcovmatrix,F.data[,c(6,7)])
sres2 <- strength_mvmr(r_input = F.data, gencov = Xcovmat)

#Step 4: Test for horizontal pleiotropy using conventional Q-statistic estimation
# 方法1 conventional 
pres <- pleiotropy_mvmr(r_input = F.data, gencov = 0)
# 方法1 conventional 
pres <- pleiotropy_mvmr(r_input = F.data, gencov = Xcovmat)

#Step 5: Estimate causal effects
res <- ivw_mvmr(r_input = F.data)

#Step 6: Robust causal effect estimation.
res1 <- qhet_mvmr(F.data, mvmrcovmatrix, CI = F, iterations = 1000)
res1

library(MendelianRandomization)

mr_mvegger(mr_mvinput(bx = cbind(F.data$betaX1, F.data$betaX2), bxse = cbind(F.data$sebetaX1, F.data$sebetaX2),
                      by = F.data$betaYG, byse = F.data$sebetaYG), orientate = 1)


