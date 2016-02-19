#-------------------------------------------------
#-------------------------------------------------
#Panel Models - KFW Grid
#Testing in Panel the impact of being treated with demarcation
#On the Max Level of NDVI, measured as the yearly max NDVI value (LTDR)
#-------------------------------------------------
#-------------------------------------------------
library(devtools)
devtools::install_github("itpir/SCI@master")
library(SCI)
library(stargazer)
library(lmtest)
library(multiwayvcov)
loadLibs()

#load in shpfile for grids

shpfile = "/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/GridDataProcessed/OhFive_gridanalysis_inputs_wpretrends.shp"
dta_Shp = readShapePoly(shpfile)

names(dta_Shp@data)[names(dta_Shp@data)=="pre_trend_"] <- "pre_trend_NDVI_max"
names(dta_Shp@data)[names(dta_Shp@data)=="pre_trend_.1"] <- "pre_trend_temp_mean"
names(dta_Shp@data)[names(dta_Shp@data)=="pre_trend_.2"] <- "pre_trend_temp_max"
names(dta_Shp@data)[names(dta_Shp@data)=="pre_trend_.3"] <- "pre_trend_temp_min"
names(dta_Shp@data)[names(dta_Shp@data)=="pre_trend_.4"] <- "pre_trend_precip_mean"
names(dta_Shp@data)[names(dta_Shp@data)=="pre_trend_.5"] <- "pre_trend_precip_max"
names(dta_Shp@data)[names(dta_Shp@data)=="pre_trend_.6"] <- "pre_trend_precip_min"

names(dta_Shp@data)[names(dta_Shp@data)=="pre_trend_.7"] <- "pre_trend_ntl"
names(dta_Shp@data)[names(dta_Shp@data)=="pre_trend_.8"] <- "pre_trend_cv"
names(dta_Shp@data)[names(dta_Shp@data)=="pre_trend_.9"] <- "pre_trend_cy"
names(dta_Shp@data)[names(dta_Shp@data)=="pre_trend_.10"] <- "pre_trend_rv"
names(dta_Shp@data)[names(dta_Shp@data)=="pre_trend_.11"] <- "pre_trend_ry"
names(dta_Shp@data)[names(dta_Shp@data)=="pre_trend_.12"] <- "pre_trend_sov"
names(dta_Shp@data)[names(dta_Shp@data)=="pre_trend_.13"] <- "pre_trend_soy"
names(dta_Shp@data)[names(dta_Shp@data)=="pre_trend_.14"] <- "pre_trend_suv"
names(dta_Shp@data)[names(dta_Shp@data)=="pre_trend_.15"] <- "pre_trend_suy"
names(dta_Shp@data)[names(dta_Shp@data)=="pre_trend_.16"] <- "pre_trend_wv"


#-----------------------------------------------
#----------------------------------------------
#Load in psm_Long untrimmed version created from KFWGrides_panelResults_Max_Pre2001_UNTRIMMED.R
#-----------------------------------------------
#----------------------------------------------

psm_Long_Untrimmed <- read.csv("/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/GridDataProcessed/psm_Long_Untrimmed.csv")

psmtest <- psm_Long_Untrimmed
dtatest <- subset(dta_Shp@data, select=c(GridID, demend_y, enforce_st))
psmtest2=merge(psmtest, dtatest, by.x="GridID", by.y="GridID")
psm_Long <- psmtest2

#Create years to demarcation
psm_Long$yrtodem <- NA
psm_Long$yrtodem=psm_Long$Year - psm_Long$demend_y

#Create correct demarcation treatment variable
psmtest3 <- psm_Long
psmtest3$trtdem <- NA
psmtest3$trtdem[which(psmtest3$Year<psmtest3$demend_y)]<-0
psmtest3$trtdem[which(psmtest3$Year>=psmtest3$demend_y)]<-1

psm_Long <- psmtest3

#create categorical variable for distance to boundary
psm_Long$HubDistCat<-0
psm_Long$HubDistCat[psm_Long$HubDist>5]<-1

#create arc of deforestation variable
psm_Long$arc<-0
psm_Long$arc[which(psm_Long$UF=="PA" | psm_Long$UF=="TO")] <- 1

#SUBSET for communities demarcated 2002 or later (to match comms in KFW_Points)

psm_Long<- psm_Long[psm_Long$demend_y>=2002,]

#-------------------------------------------------
#-------------------------------------------------
#Run Panel Models
#-------------------------------------------------
#-------------------------------------------------

pModelMax_A = lm(MaxL_ ~ trtdem + factor(reu_id),data=psm_Long, weights=terrai_are)
summary(pModelMax_A)
clusterA <- cluster.vcov(pModelMax_A,cbind(psm_Long$reu_id,psm_Long$Year),force_posdef=TRUE)
CMREG_A <- coeftest(pModelMax_A, clusterA)
print(CMREG_A)

pModelMax_B = lm(MaxL_ ~ trtdem+ MeanT_ + MeanP_ + Pop_ + MaxT_ + MaxP_ + MinT_ + MinP_ + factor(reu_id),
                 data=psm_Long, weights=terrai_are)
summary(pModelMax_B)
clusterB <- cluster.vcov(pModelMax_B,cbind(psm_Long$reu_id,psm_Long$Year),force_posdef=TRUE)
CMREG_B <- coeftest(pModelMax_B, clusterB)
print(CMREG_B)

pModelMax_C = lm(MaxL_ ~ trtdem+ MeanT_ + MeanP_ + Pop_ + MaxT_ + MaxP_ + MinT_ + MinP_ + Year + factor(reu_id),
                 data=psm_Long, weights=terrai_are)
summary(pModelMax_C)
clusterC <- cluster.vcov(pModelMax_C,cbind(psm_Long$reu_id,psm_Long$Year),force_posdef=TRUE)
CMREG_C <- coeftest(pModelMax_C, clusterC)
print(CMREG_C)

pModelMax_D = lm(MaxL_ ~ trtdem+ MeanT_ + MeanP_ + Pop_ + MaxT_ + MaxP_ + MinT_ + MinP_ + Year + 
                   trtdem*arc + factor(reu_id),
                 data=psm_Long, weights=terrai_are)
summary(pModelMax_D)
clusterD <- cluster.vcov(pModelMax_D,cbind(psm_Long$reu_id,psm_Long$Year),force_posdef=TRUE)
CMREG_D <- coeftest(pModelMax_D, clusterD)
print(CMREG_D)


## Stargazer Output

stargazer(CMREG_A,CMREG_B,CMREG_C,CMREG_D,
          type="html",align=TRUE,keep=c("trtdem","MeanT_","MeanP_","Pop_","MaxT_","MaxP_","MinT_","MinP_","Year"),
          #covariate.labels=c("TrtMnt_demend_y","MeanT","MeanP","Pop","MaxT","MaxP","MinT","MinP","Year"),
          omit.stat=c("f","ser"),
          title="Regression Results",
          dep.var.labels=c("Max NDVI")
)
