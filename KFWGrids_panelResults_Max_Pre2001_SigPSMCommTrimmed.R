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
#-------------------------------------------------
#-------------------------------------------------
#Load in Processed Data - produced from script KFW_dataMerge.r
#-------------------------------------------------
#-------------------------------------------------

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


#-------------------------------------------------
#-------------------------------------------------
#Define the Treatment Variable and Population
#-------------------------------------------------
#-------------------------------------------------
#Make a binary to test treatment..
dta_Shp@data["TrtBin"] <- 0
dta_Shp@data$TrtBin[dta_Shp@data$demend_y <= 2001] <- 1
dta_Shp@data$TrtBin[(dta_Shp@data$demend_m > 4) & (dta_Shp@data$demend_y==2001)] <- 0
trttable <- table (dta_Shp@data$TrtBin)
View(trttable)


#-----------------------------------------------
#-----------------------------------------------
#Create variables for high pressure and predicted high pressure communities based on NDVI
#-----------------------------------------------
#-----------------------------------------------

#Creating climate averages for 1982-1995

dta_Shp$MeanT_82_95 <- timeRangeAvg(dta_Shp@data, "MeanT_",1982,1995)
dta_Shp$MinT_82_95 <- timeRangeAvg(dta_Shp@data, "MinT_",1982,1995)
dta_Shp$MaxT_82_95 <- timeRangeAvg(dta_Shp@data, "MaxT_",1982,1995)
dta_Shp$MeanP_82_95 <- timeRangeAvg(dta_Shp@data, "MeanP_",1982,1995)
dta_Shp$MinP_82_95 <- timeRangeAvg(dta_Shp@data, "MinP_",1982,1995)
dta_Shp$MaxP_82_95 <- timeRangeAvg(dta_Shp@data, "MaxP_",1982,1995)

#Creating ntl and ag averages for available years

dta_Shp$ntl_92_95 <- timeRangeAvg(dta_Shp@data, "ntl_",1992,1995)
dta_Shp$cv_94_95 <- timeRangeAvg(dta_Shp@data, "cv",1994,1995)
dta_Shp$cy_91_95 <- timeRangeAvg(dta_Shp@data, "cy",1991,1995)
dta_Shp$rv_94_95 <- timeRangeAvg(dta_Shp@data, "rv", 1994,1995)
dta_Shp$ry_90_95 <- timeRangeAvg(dta_Shp@data, "ry",1990,1995)
dta_Shp$sov_94_95 <- timeRangeAvg(dta_Shp@data, "sov",1994,1995)
dta_Shp$soy_91_95 <- timeRangeAvg(dta_Shp@data, "soy",1991,1995)
dta_Shp$suv_94_95 <- timeRangeAvg(dta_Shp@data, "suv",1994,1995)
dta_Shp$suy_91_95 <- timeRangeAvg(dta_Shp@data, "suy",1991,1995)
dta_Shp$wv_94_95 <- timeRangeAvg(dta_Shp@data, "wv",1994,1995)

#Creating outcome variable options for NDVI 1982-95
#level change, 1995 minus 1982
dta_Shp$MaxL_levchange_95_82 <- dta_Shp$MaxL_1995 - dta_Shp$MaxL_1982
#negative pre_trend_NDVI_max
dta_Shp$BinNDVI=0
dta_Shp$BinNDVI[dta_Shp$pre_trend_NDVI_max<0]=1


HPModel = glm(BinNDVI ~ Pop_1995 + pre_trend_temp_mean + pre_trend_temp_min + 
                pre_trend_temp_max + pre_trend_precip_min + pre_trend_precip_mean + pre_trend_precip_max + 
                pre_trend_ntl + Slope + Elevation + Riv_Dist + urbtravtim + Road_dist + 
                log_dist + mine_dist + fedcon_dis + stcon_dist + rail_dist + 
                pre_trend_cv + pre_trend_cy + pre_trend_rv + pre_trend_ry + pre_trend_sov + pre_trend_soy +
                pre_trend_suv + pre_trend_suy + pre_trend_wv,
              family=binomial(logit), data=dta_Shp@data)
#pre-trends on right hand side, for pre_trend_NDVI and level change outcomes 
HPModel = lm(pre_trend_NDVI_max ~ Pop_1995 + pre_trend_temp_mean + pre_trend_temp_min + 
               pre_trend_temp_max + pre_trend_precip_min + pre_trend_precip_mean + pre_trend_precip_max + 
               pre_trend_ntl + Slope + Elevation + Riv_Dist + urbtravtim + Road_dist + 
               log_dist + mine_dist + fedcon_dis + stcon_dist + rail_dist + 
               pre_trend_cv + pre_trend_cy + pre_trend_rv + pre_trend_ry + pre_trend_sov + pre_trend_soy +
               pre_trend_suv + pre_trend_suy + pre_trend_wv, data=dta_Shp@data)
HPModel = lm(MaxL_levchange_95_82 ~ Pop_1995 + pre_trend_temp_mean + pre_trend_temp_min + 
               pre_trend_temp_max + pre_trend_precip_min + pre_trend_precip_mean + pre_trend_precip_max + 
               pre_trend_ntl + Slope + Elevation + Riv_Dist + urbtravtim + Road_dist + 
               log_dist + mine_dist + fedcon_dis + stcon_dist + rail_dist + 
               pre_trend_cv + pre_trend_cy + pre_trend_rv + pre_trend_ry + pre_trend_sov + pre_trend_soy +
               pre_trend_suv + pre_trend_suy + pre_trend_wv + ,data=dta_Shp@data)
#averages over pre time period on right hand side, for pre_trend and level change outcomes
HPModel = lm(pre_trend_NDVI_max ~ Pop_1995 + Slope + Elevation + Riv_Dist + urbtravtim + Road_dist + 
               log_dist + mine_dist + fedcon_dis + stcon_dist + rail_dist + 
               MeanT_82_95 + MinT_82_95 + MaxT_82_95 + MinP_82_95 + MeanP_82_95 + MaxP_82_95 +
               ntl_92_95 + cv_94_95 + cy_91_95 + rv_94_95 + ry_90_95 + suv_94_95 + 
               suy_91_95 + wv_94_95,data=dta_Shp@data)
sov_94_95 + soy_91_95 +
  HPModel = lm(MaxL_levchange_95_82 ~ Pop_1995 + Slope + Elevation + Riv_Dist + urbtravtim + Road_dist + 
                 log_dist + mine_dist + fedcon_dis + stcon_dist + rail_dist + 
                 MeanT_82_95 + MinT_82_95 + MaxT_82_95 + MinP_82_95 + MeanP_82_95 + MaxP_82_95 +
                 ntl_92_95 + cv_94_95 + cy_91_95 + rv_94_95 + ry_90_95 + sov_94_95 + soy_91_95 + suv_94_95 + 
                 suy_91_95 + wv_94_95,data=dta_Shp@data)
#pre-trends and avgs for ag, but only avgs for climate
HPModel = lm(pre_trend_NDVI_max ~ Pop_1995 + 
               pre_trend_ntl + Slope + Elevation + Riv_Dist + urbtravtim + Road_dist + 
               log_dist + mine_dist + fedcon_dis + stcon_dist + rail_dist + 
               pre_trend_cv + pre_trend_cy + pre_trend_ry + pre_trend_sov  +
               pre_trend_suv + pre_trend_suy + 
               MeanT_82_95 + MinT_82_95 + MaxT_82_95 + MinP_82_95 + MeanP_82_95 + MaxP_82_95 +
               ntl_92_95 + cv_94_95 + cy_91_95 + + ry_90_95 + suv_94_95 + 
               suy_91_95,data=dta_Shp@data)
#everything, pre_trend and level change
HPModel = lm(pre_trend_NDVI_max ~ Pop_1995 + pre_trend_temp_mean + pre_trend_temp_min + 
               pre_trend_temp_max + pre_trend_precip_min + pre_trend_precip_mean + pre_trend_precip_max + 
               pre_trend_ntl + Slope + Elevation + Riv_Dist + urbtravtim + Road_dist + 
               log_dist + mine_dist + fedcon_dis + stcon_dist + rail_dist + 
               pre_trend_cv + pre_trend_cy + pre_trend_ry + pre_trend_sov  +
               pre_trend_suv + pre_trend_suy + 
               MeanT_82_95 + MinT_82_95 + MaxT_82_95 + MinP_82_95 + MeanP_82_95 + MaxP_82_95 + ntl_92_95 +
               cv_94_95 + cy_91_95 + ry_90_95 + suv_94_95 + suy_91_95,
             data=dta_Shp@data)
#pre_trend_wv sov_94_95 + soy_91_95 +rv_94_95 wv_94_95pre_trend_rv
#HPModel = lm(MaxL_levchange_95_82 ~ Pop_1995 + pre_trend_temp_mean + pre_trend_temp_min + 
pre_trend_temp_max + pre_trend_precip_min + pre_trend_precip_mean + pre_trend_precip_max + 
  pre_trend_ntl + Slope + Elevation + Riv_Dist + urbtravtim + Road_dist + 
  log_dist + mine_dist + fedcon_dis + stcon_dist + rail_dist + 
  pre_trend_cv + pre_trend_cy + pre_trend_rv + pre_trend_ry + pre_trend_sov + pre_trend_soy +
  pre_trend_suv + pre_trend_suy + pre_trend_wv + 
  MeanT_82_95 + MinT_82_95 + MaxT_82_95 + MinP_82_95 + MeanP_82_95 + MaxP_82_95 +
  ntl_92_95 + cv_94_95 + cy_91_95 + rv_94_95 + ry_90_95 + sov_94_95 + soy_91_95 + suv_94_95 + 
  suy_91_95 + wv_94_95,data=dta_Shp@data)

#Running the model with cmreg
HPModel$Id <- cluster.vcov(HPModel,c(dta_Shp@data$Id))
CMREG <- coeftest(HPModel, HPModel$Id)
print(CMREG)
summary(HPModel)

#------------------------------------------------
#------------------------------------------------
#Predicted pre_trend_NDVI_max (to interact with treatment binary)
#-----------------------------------------------
#-----------------------------------------------

dta_Shp@data$model_int_early_1 <- CMREG[1]
dta_Shp@data$model_int_early_2 <- CMREG[2] * dta_Shp@data$Pop_1995
dta_Shp@data$model_int_early_3 <- CMREG[3] * dta_Shp@data$pre_trend_temp_mean
dta_Shp@data$model_int_early_4 <- CMREG[4] * dta_Shp@data$pre_trend_temp_min
dta_Shp@data$model_int_early_5 <- CMREG[5] * dta_Shp@data$pre_trend_temp_max
dta_Shp@data$model_int_early_6 <- CMREG[6] * dta_Shp@data$pre_trend_precip_min
dta_Shp@data$model_int_early_7 <- CMREG[7] * dta_Shp@data$pre_trend_precip_mean
dta_Shp@data$model_int_early_8 <- CMREG[8] * dta_Shp@data$pre_trend_precip_max
dta_Shp@data$model_int_early_9 <- CMREG[9] * dta_Shp@data$pre_trend_ntl
dta_Shp@data$model_int_early_10 <- CMREG[10] * dta_Shp@data$Slope
dta_Shp@data$model_int_early_11 <- CMREG[11] * dta_Shp@data$Elevation
dta_Shp@data$model_int_early_12 <- CMREG[12] * dta_Shp@data$Riv_Dist
dta_Shp@data$model_int_early_13 <- CMREG[13] * dta_Shp@data$urbtravtim
dta_Shp@data$model_int_early_14 <- CMREG[14] * dta_Shp@data$Road_dist
dta_Shp@data$model_int_early_15 <- CMREG[15] * dta_Shp@data$log_dist
dta_Shp@data$model_int_early_16 <- CMREG[16] * dta_Shp@data$mine_dist
dta_Shp@data$model_int_early_17 <- CMREG[17] * dta_Shp@data$fedcon_dis
dta_Shp@data$model_int_early_18 <- CMREG[18] * dta_Shp@data$stcon_dist
dta_Shp@data$model_int_early_19 <- CMREG[19] * dta_Shp@data$rail_dist
dta_Shp@data$model_int_early_20 <- CMREG[20] * dta_Shp@data$pre_trend_cv
dta_Shp@data$model_int_early_21 <- CMREG[21] * dta_Shp@data$pre_trend_cy
dta_Shp@data$model_int_early_22 <- CMREG[22] * dta_Shp@data$pre_trend_ry
dta_Shp@data$model_int_early_23 <- CMREG[23] * dta_Shp@data$pre_trend_sov
dta_Shp@data$model_int_early_24 <- CMREG[24] * dta_Shp@data$pre_trend_suv
dta_Shp@data$model_int_early_25 <- CMREG[25] * dta_Shp@data$pre_trend_suy
dta_Shp@data$model_int_early_26 <- CMREG[26] * dta_Shp@data$MeanT_82_95
dta_Shp@data$model_int_early_27 <- CMREG[27] * dta_Shp@data$MinT_82_95
dta_Shp@data$model_int_early_28 <- CMREG[28] * dta_Shp@data$MaxT_82_95
dta_Shp@data$model_int_early_29 <- CMREG[29] * dta_Shp@data$MinP_82_95
dta_Shp@data$model_int_early_30 <- CMREG[30] * dta_Shp@data$MeanP_82_95
dta_Shp@data$model_int_early_31 <- CMREG[31] * dta_Shp@data$MaxP_82_95
dta_Shp@data$model_int_early_32 <- CMREG[32] * dta_Shp@data$ntl_92_95
dta_Shp@data$model_int_early_33 <- CMREG[33] * dta_Shp@data$cv_94_95
dta_Shp@data$model_int_early_34 <- CMREG[34] * dta_Shp@data$cy_91_95
dta_Shp@data$model_int_early_35 <- CMREG[35] * dta_Shp@data$ry_90_95
dta_Shp@data$model_int_early_36 <- CMREG[36] * dta_Shp@data$suv_94_95
dta_Shp@data$model_int_early_37 <- CMREG[37] * dta_Shp@data$suy_91_95


dta_Shp@data$predict_NDVI_max_pre <- dta_Shp@data$model_int_early_1+dta_Shp@data$model_int_early_2+dta_Shp@data$model_int_early_3+dta_Shp@data$model_int_early_4+dta_Shp@data$model_int_early_5+dta_Shp@data$model_int_early_6+
  dta_Shp@data$model_int_early_7+dta_Shp@data$model_int_early_8+dta_Shp@data$model_int_early_9+dta_Shp@data$model_int_early_10+dta_Shp@data$model_int_early_11+dta_Shp@data$model_int_early_12+
  dta_Shp@data$model_int_early_13+dta_Shp@data$model_int_early_14+dta_Shp@data$model_int_early_15+dta_Shp@data$model_int_early_16+dta_Shp@data$model_int_early_17+dta_Shp@data$model_int_early_18+dta_Shp@data$model_int_early_19+
  dta_Shp@data$model_int_early_20+dta_Shp@data$model_int_early_21+dta_Shp@data$model_int_early_22+dta_Shp@data$model_int_early_23+dta_Shp@data$model_int_early_24+dta_Shp@data$model_int_early_25+dta_Shp@data$model_int_early_26+
  dta_Shp@data$model_int_early_27+dta_Shp@data$model_int_early_28+dta_Shp@data$model_int_early_29+dta_Shp@data$model_int_early_30+dta_Shp@data$model_int_early_31+dta_Shp@data$model_int_early_32+
  dta_Shp@data$model_int_early_33+dta_Shp@data$model_int_early_34+dta_Shp@data$model_int_early_35+dta_Shp@data$model_int_early_36+dta_Shp@data$model_int_early_37


#-------------------------------------------------
#-------------------------------------------------
#Define and run the first-stage of the PSM, calculating propensity scores
#-------------------------------------------------
#-------------------------------------------------
psmModel <-  "TrtBin ~ terrai_are + Pop_1995 + MeanT_1995 + pre_trend_temp_mean + pre_trend_temp_min + 
pre_trend_temp_max + MeanP_1995 + pre_trend_precip_min + pre_trend_NDVI_max + ntl_1995 +Slope + Elevation + 
MaxL_1995 + Riv_Dist + Road_dist + pre_trend_precip_mean + pre_trend_precip_max"

psmRes <- SCI::SpatialCausalPSM(dta_Shp,mtd="logit",psmModel,
                                drop="none",
                                visual=TRUE)

dta_Shp_psm = psmRes$data

#Creating categorical measure of pre_trend_NDVI_max to interact with treatment binary
#(where lower than median value of pre_trend_NDVI_max is considered high pressure)

pretrend_NDVI_median<-fivenum(dta_Shp_psm$pre_trend_NDVI_max)[3]
dta_Shp_psm$pre_trend_NDVI_max_cat <- NA
dta_Shp_psm$pre_trend_NDVI_max_cat <-ifelse(dta_Shp_psm$pre_trend_NDVI_max<pretrend_NDVI_median,1,0)

predict_NDVI_median<-fivenum(dta_Shp$predict_NDVI_max_pre)[3]
dta_Shp_psm$predict_NDVI_max_pre_cat <- NA
dta_Shp_psm$predict_NDVI_max_pre_cat <-ifelse(dta_Shp_psm$predict_NDVI_max_pre<predict_NDVI_median,1,0)


#-------------------------------------------------
#-------------------------------------------------
#Based on the Propensity Score Matches, pair comprable treatment and control units.
#-------------------------------------------------
#-------------------------------------------------
#drop_set<- c(drop_unmatched=TRUE,drop_method="None",drop_thresh=0.25)
#psm_Pairs <- SAT(dta = psmRes$data, mtd = "fastNN",constraints=c(groups="UF"),psm_eq = psmModel, ids = "GridID", drop_opts = drop_set, visual="TRUE", TrtBinColName="TrtBin")

#trttable <- table (psm_Pairs@data$TrtBin)
#View(trttable)

#-------------------------------------------------
#-------------------------------------------------
#Convert from a wide-form dataset for the Cross-sectional 
#to a long-form dataset for the panel model.
#-------------------------------------------------
#-------------------------------------------------

# varList=c("MaxL_")
# psm_Long <- BuildTimeSeries(dta=dta_Shp_psm,idField="GridID",varList_pre=varList,1982,2010,colYears=c("demend_y","apprend_y","regend_y"),
#                             interpYears=c("Slope","Road_dist","Riv_Dist","UF","Elevation","terrai_are","Pop_","MeanT_","MeanP_","MaxT_",
#                                           "MaxP_","MinP_","MinT_","ntl_", "fedcon_dis", "stcon_dist", "log_dist", "mine_dist", "rail_dist",
#                                           "urbtravtim", "pre_trend_NDVI_max","predict_NDVI_max_pre", "pre_trend_NDVI_max_cat",
#                                           "predict_NDVI_max_pre_cat", "reu_id", "Id" ))
# psm_Long$Year <- as.numeric(psm_Long$Year)
# 
# write.csv(psm_Long,file="/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/GridDataProcessed/psm_Long.csv")

varList=c("MaxL_")
psm_Long_Untrimmed <- BuildTimeSeries(dta=dta_Shp,idField="GridID",varList_pre=varList,1982,2010,colYears=c("demend_y","enforce_st"),
                                      interpYears=c("Slope","Road_dist","Riv_Dist","UF","Elevation","terrai_are","Pop_","MeanT_","MeanP_","MaxT_",
                                                    "MaxP_","MinP_","MinT_", "reu_id", "Id" ))
psm_Long_Untrimmed$Year <- as.numeric(psm_Long_Untrimmed$Year)

write.csv(psm_Long_Untrimmed,file="/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/GridDataProcessed/psm_Long_Untrimmed.csv")
write.csv(psm_Long_5yr,file="/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/GridDataProcessed/psm_Long_Untrimmed_5Yr.csv")

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

#Check new dem treatment variable
psm_demyear <- psm_Long
psm_demyear <- psm_Long[psm_Long$Year==psm_Long$demend_y,]
#this should be equal to 0:
summary(psm_demyear$yrtodem)

#Create correct enforcement treatment variable

#change 1 community (reu_id=84) where enforcement starts 1 year before demarcation to start in year of demarcation
psmtest4 <- psm_Long
psmtest4$enfdiff= psmtest4$enforce_st - psmtest4$demend_y
table(psmtest4$enfdiff)
psmenf <- subset(psmtest4, psmtest4$enfdiff<0)
table(psmenf$reu_id)

psmtest4$enforce_st[which(psmtest4$reu_id==84)]<-2007
psm_Long <- psmtest4

#create enforcement treatment var
psmtest5 <- psm_Long
psmtest5$trtenf <- 0
psmtest5$trtenf[which(psmtest5$Year>=psmtest5$enforce_st)]<-1

psm_Long <- psmtest5

#Create subset that only includes years within -5 and +5 years of demarcation
# psm_Long_5yr <- psm_Long
# test <- psm_Long_5yr[psm_Long_5yr$yrtodem>=-5,]
# test1 <- test[test$yrtodem<=5,]
# psm_Long <- test

#Create subset that only includes reu_ids for the pairs made from 1st stage PSM (at community level) with 
# only the predictors that were significant (slope, elevation, road distance, pre-trend in min annual temp) at the COMMUNITY level
psm_Long_sigpsm <- psm_Long
psm_Long_sigpsm1 <- psm_Long_sigpsm[psm_Long_sigpsm$reu_id %in% c(131,114,118,142,117, 121, 105, 148,  
                                                                  93, 107, 152, 150, 154, 112, 158, 159, 
                                                                  160, 161, 162, 163, 164, 146, 110, 180, 
                                                                  168, 151, 157,173, 176, 115,  80,  92,74,
                                                                  119, 132,  88, 128, 155, 129, 156, 100, 123,
                                                                  106, 172,  87,  78,  73, 122, 169, 144, 133, 111,  85,  89,
                                                                  79, 86,  91, 175,  81,  82, 125, 126, 141,  96, 109, 
                                                                  103, 143, 137, 135, 136, 134, 179, 178,  95),]

psm_Long <- psm_Long_sigpsm1

## Run Models

pModelMax_A <- "MaxL_ ~ trtdem + trtenf + factor(reu_id)"
pModelMax_B <- "MaxL_ ~ ttrtdem + trtenf  + Pop_ + MeanT_ + MeanP_ + MaxT_ + MaxP_ + MinT_ + MinP_  + factor(reu_id) "
pModelMax_C <- "MaxL_ ~ trtdem + trtenf + Pop_ + MeanT_ + MeanP_ + MaxT_ + MaxP_ + MinT_ + MinP_  + Year + factor(reu_id)"
pModelMax_C1 <- "MaxL_ ~ trtdem + Pop_ + MeanT_ + MeanP_+ MaxT_ + MaxP_ + MinT_ + MinP_  + factor(Year) + factor(reu_id)"
pModelMax_C2 <- "MaxL_ ~ trtdem + trtenf + Pop_ + MeanT_ + MeanP_+ MaxT_ + MaxP_ + MinT_ + MinP_  + factor(Year) + factor(reu_id)"



pModelMax_A_fit <- Stage2PSM(pModelMax_A ,psm_Long,type="cmreg", table_out=TRUE,opts=c("reu_id","Year"))
pModelMax_B_fit <- Stage2PSM(pModelMax_B ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMax_C_fit <- Stage2PSM(pModelMax_C ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMax_C1_fit <- Stage2PSM(pModelMax_C1 ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMax_C2_fit <- Stage2PSM(pModelMax_C2 ,psm_Long,type="cmreg", table_out=TRUE,opts=c("reu_id","Year"))


##Stargazer

stargazer(pModelMax_A_fit$cmreg,pModelMax_B_fit$cmreg,pModelMax_C_fit$cmreg,
          pModelMax_C1_fit$cmreg,pModelMax_C2_fit$cmreg,
          type="html", align=TRUE,
          keep=c("TrtMnt","Pop","Mean","Max","Min","Year"),
          covariate.labels=c("Treatment (Demarcation)","Treatment (Demarcation + Enforcement Support)","Population","Mean Temp",
                             "Mean Precip","Max Temp","Max Precip","Min Temp","Min Precip","Year"),
          omit.stat=c("f","ser"),
          add.lines=list(c("Observations","2146","2146","2146","2146","2146"),
                         c("Community Fixed Effects?","Yes","Yes","Yes","Yes","Yes","Yes"),
                         c("Year Fixed Effects?","No","No","No","Yes","Yes")),
          title="Regression Results",
          dep.var.labels=c("Max NDVI"))
