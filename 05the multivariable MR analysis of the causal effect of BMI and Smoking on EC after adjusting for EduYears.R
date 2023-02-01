library(devtools)
library(MRInstruments)
library(TwoSampleMR) #加载R包
#install_github("WSpiller/MVMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
library(MVMR)
# install_github("WSpiller/MRPracticals",build_opts = c("--no-resave-data", "--no-manual"),build_vignettes = TRUE)
library(MRPracticals)
vignette("MVMR Tutorial")


rawdat_mvmr <- read.table('04edu_BMI_smoke_cancer.txt',header=T,sep="\t",check.names = F,stringsAsFactors = F)


F.data<-format_mvmr(BXGs=rawdat_mvmr[,c(4,6,8)],
                    BYG=rawdat_mvmr[,10],
                    seBXGs=rawdat_mvmr[,c(5,7,9)],
                    seBYG=rawdat_mvmr[,11],
                    RSID=rawdat_mvmr[,1])
res <- ivw_mvmr(r_input = F.data)
head(F.data)
#  Step 3: Test for weak instruments用
# Method 1 conventional F-statistic
sres<-strength_mvmr(r_input=F.data,gencov=0)
#  Method 2 conditional F-statistics
mvmrcovmatrix<-matrix(c(1,-0.1,-0.05,-0.1,1,0.2,-0.05,0.2,1), nrow = 3, ncol = 3)
#Xcovmat<-phenocov_mvmr(mvmrcovmatrix,F.data[,c(6,7)])
#sres2 <- strength_mvmr(r_input = F.data, gencov = Xcovmat)

#Step 4: Test for horizontal pleiotropy using conventional Q-statistic estimation
# 方法1 conventional 
pres <- pleiotropy_mvmr(r_input = F.data, gencov = 0)
# 方法1 conventional 
pres <- pleiotropy_mvmr(r_input = F.data, gencov = Xcovmat)

#Step 5: Estimate causal effects用
res <- ivw_mvmr(r_input = F.data)

#Step 6: Robust causal effect estimation.
res1 <- qhet_mvmr(F.data, mvmrcovmatrix, CI = F, iterations = 1000)

res1

#mr_mvegger用
library(MendelianRandomization)

mr_mvegger(mr_mvinput(bx = cbind(F.data$betaX1, F.data$betaX2,  F.data$betaX3), bxse = cbind(F.data$sebetaX1, F.data$sebetaX2, F.data$sebetaX3),
                      by = F.data$betaYG, byse = F.data$sebetaYG), orientate = 1)

