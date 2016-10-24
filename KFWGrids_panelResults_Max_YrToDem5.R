

#-------------------------------------------------
#-------------------------------------------------
#Panel Model - KFW Grid
#Testing in Panel the impact of demarcation treatment, grouping by years pre- or post-demarcation
#Outcome: Max Level of NDVI, measured as the yearly max NDVI value (LTDR)
#-------------------------------------------------
#-------------------------------------------------

library(devtools)
devtools::install_github("itpir/SCI@master")
library(SCI)
library(stargazer)
library(lmtest)
library(multiwayvcov)
loadLibs()


#---------------
#Dataset
#---------------

#Load grid panel dataset 
#"trimmed" - dropped cells from communities outside of common support for full propensity score matching

psm_Long <- read.csv("/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/GridDataProcessed/psm_Long_TRTDEM.csv")

#drop incorrect treatment vars (something went wrong when panel dataset created from shapefile)
#correct treatment variable is "trtdem" for demarcation or "trtenf" for existence of enforcement support
psm_Long<-psm_Long[-c(32:34)]

#Create years to demarcation variable that groups everything 5 years or longer (pre- or post- demarcation) into one category
#yrtodem is positive post-demarcation (e.g. 3 means 3 years after demarcation)
psm_Long$yrtodem5<-psm_Long$yrtodem
psm_Long$yrtodem5[psm_Long$yrtodem>=5]<-5
psm_Long$yrtodem5[psm_Long$yrtodem<=-5]<--5


#-----------------
# Analytic Models
#-----------------
#standard errors clustered by community and year

#replacing TrtMnt_demend_y with factor(yrtodem)
pModelMax_A <- "MaxL_ ~ factor(yrtodem) + factor(reu_id)"
pModelMax_B <- "MaxL_ ~ factor(yrtodem) + Pop_ + MeanT_ + MeanP_ +MaxT_ + MaxP_ + MinT_ + MinP_ + factor(reu_id) "

pModelMax_C <- "MaxL_ ~ factor(yrtodem) + Pop_ + MeanT_ + MeanP_ +MaxT_ + MaxP_ + MinT_ + MinP_ + Year + factor(reu_id)"
pModelMax_D <- "MaxL_ ~ factor(yrtodem) + Pop_ +MeanT_ + MeanP_ +MaxT_ + MaxP_ + MinT_ + MinP_ + factor(Year) + factor(reu_id)"


pModelMax_A_fit <- Stage2PSM(pModelMax_A ,psm_Long,type="lm", table_out=TRUE, opts=c("reu_id","Year"))
pModelMax_B_fit <- Stage2PSM(pModelMax_B ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMax_C_fit <- Stage2PSM(pModelMax_C ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMax_D_fit <- Stage2PSM(pModelMax_D ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
