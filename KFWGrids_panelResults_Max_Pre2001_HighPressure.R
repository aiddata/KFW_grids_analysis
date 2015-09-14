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
               ntl_92_95 + cy_94_95 + cv_91_95 + rv_94_95 + ry_90_95 + sov_94_95 + soy_91_95 + suv_94_95 + 
               suy_91_95 + wv_94_95,data=dta_Shp@data)
HPModel = lm(MaxL_levchange_95_82 ~ Pop_1995 + Slope + Elevation + Riv_Dist + urbtravtim + Road_dist + 
               log_dist + mine_dist + fedcon_dis + stcon_dist + rail_dist + 
               MeanT_82_95 + MinT_82_95 + MaxT_82_95 + MinP_82_95 + MeanP_82_95 + MaxP_82_95 +
               ntl_92_95 + cy_94_95 + cv_91_95 + rv_94_95 + ry_90_95 + sov_94_95 + soy_91_95 + suv_94_95 + 
               suy_91_95 + wv_94_95,data=dta_Shp@data)
#everything, pre_trend and level change
HPModel = lm(pre_trend_NDVI_max ~ Pop_1995 + pre_trend_temp_mean + pre_trend_temp_min + 
               pre_trend_temp_max + pre_trend_precip_min + pre_trend_precip_mean + pre_trend_precip_max + 
               pre_trend_ntl + Slope + Elevation + Riv_Dist + urbtravtim + Road_dist + 
               log_dist + mine_dist + fedcon_dis + stcon_dist + rail_dist + 
               pre_trend_cv + pre_trend_cy + pre_trend_rv + pre_trend_ry + pre_trend_sov + pre_trend_soy +
               pre_trend_suv + pre_trend_suy + pre_trend_wv + 
               MeanT_82_95 + MinT_82_95 + MaxT_82_95 + MinP_82_95 + MeanP_82_95 + MaxP_82_95 +
               ntl_92_95 + cy_94_95 + cv_91_95 + rv_94_95 + ry_90_95 + sov_94_95 + soy_91_95 + suv_94_95 + 
               suy_91_95 + wv_94_95,data=dta_Shp@data)
HPModel = lm(MaxL_levchange_95_82 ~ Pop_1995 + pre_trend_temp_mean + pre_trend_temp_min + 
               pre_trend_temp_max + pre_trend_precip_min + pre_trend_precip_mean + pre_trend_precip_max + 
               pre_trend_ntl + Slope + Elevation + Riv_Dist + urbtravtim + Road_dist + 
               log_dist + mine_dist + fedcon_dis + stcon_dist + rail_dist + 
               pre_trend_cv + pre_trend_cy + pre_trend_rv + pre_trend_ry + pre_trend_sov + pre_trend_soy +
               pre_trend_suv + pre_trend_suy + pre_trend_wv + 
               MeanT_82_95 + MinT_82_95 + MaxT_82_95 + MinP_82_95 + MeanP_82_95 + MaxP_82_95 +
               ntl_92_95 + cy_94_95 + cv_91_95 + rv_94_95 + ry_90_95 + sov_94_95 + soy_91_95 + suv_94_95 + 
               suy_91_95 + wv_94_95,data=dta_Shp@data)


#Running the model with cmreg
HPModel$Id <- cluster.vcov(HPModel,c(dta$Id))
CMREG <- coeftest(HPModel, HPModel$Id)
print(CMREG)
summary(HPModel)
summary(HPModel)$r.squared


#-------------------------------------------------
#-------------------------------------------------
#Define and run the first-stage of the PSM, calculating propensity scores
#-------------------------------------------------
#-------------------------------------------------
psmModel <-  "TrtBin ~ terrai_are + Pop_1995 + MeanT_1995 + pre_trend_temp_mean + pre_trend_temp_min + 
pre_trend_temp_max + MeanP_1995 + pre_trend_precip_min + pre_trend_NDVI_max + ntl_1995 +Slope + Elevation + 
MaxL_1995 + Riv_Dist + Road_dist + pre_trend_precip_mean + pre_trend_precip_max"

psmRes <- SCI::SpatialCausalPSM(dta_Shp,mtd="logit",psmModel,drop="support",visual=TRUE)


#-------------------------------------------------
#-------------------------------------------------
#Based on the Propensity Score Matches, pair comprable treatment and control units.
#-------------------------------------------------
#-------------------------------------------------
drop_set<- c(drop_unmatched=TRUE,drop_method="None",drop_thresh=0.25)
psm_Pairs <- SAT(dta = psmRes$data, mtd = "fastNN",constraints=c(groups="UF"),psm_eq = psmModel, ids = "GridID", drop_opts = drop_set, visual="TRUE", TrtBinColName="TrtBin")

trttable <- table (psm_Pairs@data$TrtBin)
View(trttable)

#Creating categorical measure of pre_trend_NDVI_max to interact with treatment binary
#(where <median value of pre_trend_NDVI_max is high pressure)

pretrend_NDVI_median<-fivenum(dta_Shp$pre_trend_NDVI_max)[3]
psm_Long$pre_trend_NDVI_max_cat <- NA
psm_Long$pre_trend_NDVI_max_cat <-ifelse(psm_Long$pre_trend_NDVI_max<pretrend_NDVI_median,1,0)

predict_NDVI_median<-fivenum(psm_Long$predict_NDVI_max_pre.y)[3]
psm_Long$predict_NDVI_max_pre_cat <- NA
psm_Long$predict_NDVI_max_pre_cat <-ifelse(psm_Long$predict_NDVI_max_pre.y<predict_NDVI_median,1,0)

#-------------------------------------------------
#-------------------------------------------------
#Convert from a wide-form dataset for the Cross-sectional 
#to a long-form dataset for the panel model.
#-------------------------------------------------
#-------------------------------------------------
#Clean up data entry
#psm_Pairs$enforce_st[psm_Pairs$enforce_st == "1998-1999"] <- NA
#psm_Pairs$enforce_st <- as.numeric(paste(psm_Pairs$enforce_st))

varList = c("MaxL_")
psm_Long <- BuildTimeSeries(dta=psm_Pairs,idField="GridID",varList_pre=varList,1982,2010,colYears=c("demend_y","apprend_y","regend_y"),interpYears=c("Slope","Road_dist","Riv_Dist","UF","Elevation","terrai_are","Pop_","MeanT_","MeanP_","MaxT_","MaxP_","MinP_","MinT_"))
psm_Long$Year <- as.numeric(psm_Long$Year)

write.csv(psm_Long,file="/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/GridDataProcessed/psm_Long.csv")


psm_Long_HP <- BuildTimeSeries(dta=psm_Pairs,idField="GridID",varList_pre=varList,1982,2010,colYears=c("demend_y","apprend_y","regend_y"),
                               interpYears=c("Slope","Road_dist","Riv_Dist","UF","Elevation","terrai_are","Pop_","MeanT_","MeanP_","MaxT_",
                                             "MaxP_","MinP_","MinT_","ntl_", "fedcon_dis", "stcon_dist", "log_dist", "mine_dist", "rail_dist",
                                             "urbtravtim", "pre_trend_NDVI_max", "reu_id", "Id" ))
psm_Long_HP$Year <- as.numeric(psm_Long_HP$Year)

write.csv(psm_Long_HP,file="/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/GridDataProcessed/psm_Long_HP.csv")


dta_Shp_pre <-subset(dta_Shp, select=c(GridID, pre_trend_NDVI_max, pre_trend_temp_mean, pre_trend_precip_mean))
psm_Long_pre=merge(dta_Shp_pre@data, psm_Long_HP, by.x="GridID", by.y="GridID")
psm_Long <- psm_Long_pre

psm_Long$pre_trend_NDVI_cat=0
psm_Long$pre_trend_NDVI_cat[psm_Long$pre_trend_NDVI_max<=30.39]=1

pModelMax_A <- "MaxL_ ~ TrtMnt_demend_y + factor(reu_id)"
pModelMax_B <- "MaxL_ ~ TrtMnt_demend_y + MeanT_ + MeanP_ + Pop_ + MaxT_ + MaxP_ + MinT_ + MinP_  + factor(reu_id) "
pModelMax_C <- "MaxL_ ~ TrtMnt_demend_y + MeanT_ + MeanP_ + Pop_ + MaxT_ + MaxP_ + MinT_ + MinP_  + factor(reu_id) + Year"
pModelMax_D <- "MaxL_ ~ TrtMnt_demend_y + MeanT_ + MeanP_ + Pop_ + MaxT_ + MaxP_ + MinT_ + MinP_ + pre_trend_NDVI_cat*TrtMnt_demend_y + factor(reu_id) + Year"

pModelMax_A_fit <- Stage2PSM(pModelMax_A ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMax_B_fit <- Stage2PSM(pModelMax_B ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMax_C_fit <- Stage2PSM(pModelMax_C ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMax_D_fit <- Stage2PSM(pModelMax_D ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))


#------------------------------------------------------------------------
#------------------------------------------------------------------------
## Creating predicted high pressure regions






# Creating predicted NDVI values to interact with TrtBin

dta_Shp2@data$model_int_early_1 <- CMREG[1]
dta_Shp2@data$model_int_early_2 <- CMREG[2] * dta_Shp2@data$Pop_1990
dta_Shp2@data$model_int_early_3 <- CMREG[3] * dta_Shp2@data$Pop_1995
dta_Shp2@data$model_int_early_4 <- CMREG[4] * dta_Shp2@data$MeanT_82_95
dta_Shp2@data$model_int_early_5 <- CMREG[5] * dta_Shp2@data$MinT_82_95
dta_Shp2@data$model_int_early_6 <- CMREG[6] * dta_Shp2@data$MaxT_82_95
dta_Shp2@data$model_int_early_7 <- CMREG[7] * dta_Shp2@data$MinP_82_95
dta_Shp2@data$model_int_early_8 <- CMREG[8] * dta_Shp2@data$MeanP_82_95
dta_Shp2@data$model_int_early_9 <- CMREG[9] * dta_Shp2@data$MaxP_82_95
dta_Shp2@data$model_int_early_10 <- CMREG[10] * dta_Shp2@data$Slope
dta_Shp2@data$model_int_early_11 <- CMREG[11] * dta_Shp2@data$Elevation
dta_Shp2@data$model_int_early_12 <- CMREG[12] * dta_Shp2@data$Riv_Dist
dta_Shp2@data$model_int_early_13 <- CMREG[13] * dta_Shp2@data$ntl_1992
dta_Shp2@data$model_int_early_14 <- CMREG[14] * dta_Shp2@data$ntl_1993
dta_Shp2@data$model_int_early_15 <- CMREG[15] * dta_Shp2@data$ntl_1994
dta_Shp2@data$model_int_early_16 <- CMREG[16] * dta_Shp2@data$ntl_1995
dta_Shp2@data$model_int_early_17 <- CMREG[17] * dta_Shp2@data$urbtravtim
dta_Shp2@data$model_int_early_18 <- CMREG[18] * dta_Shp2@data$Road_dist
dta_Shp2@data$model_int_early_19 <- CMREG[19] * dta_Shp2@data$log_dist
dta_Shp2@data$model_int_early_20 <- CMREG[20] * dta_Shp2@data$mine_dist
dta_Shp2@data$model_int_early_21 <- CMREG[21] * dta_Shp2@data$fedcon_dis
dta_Shp2@data$model_int_early_22 <- CMREG[22] * dta_Shp2@data$stcon_dist
dta_Shp2@data$model_int_early_23 <- CMREG[23] * dta_Shp2@data$rail_dist
dta_Shp2@data$model_int_early_24 <- CMREG[24] * dta_Shp2@data$cv1995
dta_Shp2@data$model_int_early_25 <- CMREG[25] * dta_Shp2@data$cv1994
dta_Shp2@data$model_int_early_26 <- CMREG[26] * dta_Shp2@data$cy1995
dta_Shp2@data$model_int_early_27 <- CMREG[27] * dta_Shp2@data$cy1994
dta_Shp2@data$model_int_early_28 <- CMREG[28] * dta_Shp2@data$cy1993
dta_Shp2@data$model_int_early_29 <- CMREG[29] * dta_Shp2@data$cy1992
dta_Shp2@data$model_int_early_30 <- CMREG[30] * dta_Shp2@data$cy1991
dta_Shp2@data$model_int_early_31 <- CMREG[31] * dta_Shp2@data$rv1995
dta_Shp2@data$model_int_early_32 <- CMREG[32] * dta_Shp2@data$rv1994
dta_Shp2@data$model_int_early_33 <- CMREG[33] * dta_Shp2@data$ry1995
dta_Shp2@data$model_int_early_34 <- CMREG[34] * dta_Shp2@data$ry1994
dta_Shp2@data$model_int_early_35 <- CMREG[35] * dta_Shp2@data$ry1993
dta_Shp2@data$model_int_early_36 <- CMREG[36] * dta_Shp2@data$ry1992
dta_Shp2@data$model_int_early_37 <- CMREG[37] * dta_Shp2@data$ry1991
dta_Shp2@data$model_int_early_38 <- CMREG[38] * dta_Shp2@data$ry1990
dta_Shp2@data$model_int_early_39 <- CMREG[39] * dta_Shp2@data$sov1995
dta_Shp2@data$model_int_early_40 <- CMREG[40] * dta_Shp2@data$suv1995
dta_Shp2@data$model_int_early_41 <- CMREG[41] * dta_Shp2@data$suv1994
dta_Shp2@data$model_int_early_42 <- CMREG[42] * dta_Shp2@data$suy1995
dta_Shp2@data$model_int_early_43 <- CMREG[43] * dta_Shp2@data$suy1994
dta_Shp2@data$model_int_early_44 <- CMREG[44] * dta_Shp2@data$suy1993
dta_Shp2@data$model_int_early_45 <- CMREG[45] * dta_Shp2@data$suy1992
dta_Shp2@data$model_int_early_46 <- CMREG[46] * dta_Shp2@data$suy1991


dta_Shp2@data$predict_NDVI_max_pre <- dta_Shp2@data$model_int_early_1+dta_Shp2@data$model_int_early_2+dta_Shp2@data$model_int_early_3+dta_Shp2@data$model_int_early_4+dta_Shp2@data$model_int_early_5+dta_Shp2@data$model_int_early_6+
  dta_Shp2@data$model_int_early_7+dta_Shp2@data$model_int_early_8+dta_Shp2@data$model_int_early_9+dta_Shp2@data$model_int_early_10+dta_Shp2@data$model_int_early_11+dta_Shp2@data$model_int_early_12+
  dta_Shp2@data$model_int_early_13+dta_Shp2@data$model_int_early_14+dta_Shp2@data$model_int_early_15+dta_Shp2@data$model_int_early_16+dta_Shp2@data$model_int_early_17+dta_Shp2@data$model_int_early_18+dta_Shp2@data$model_int_early_19+
  dta_Shp2@data$model_int_early_20+dta_Shp2@data$model_int_early_21+dta_Shp2@data$model_int_early_22+dta_Shp2@data$model_int_early_23+dta_Shp2@data$model_int_early_24+dta_Shp2@data$model_int_early_25+dta_Shp2@data$model_int_early_26+
  dta_Shp2@data$model_int_early_27+dta_Shp2@data$model_int_early_28+dta_Shp2@data$model_int_early_29+dta_Shp2@data$model_int_early_30+dta_Shp2@data$model_int_early_31+dta_Shp2@data$model_int_early_32+
  dta_Shp2@data$model_int_early_33+dta_Shp2@data$model_int_early_34+dta_Shp2@data$model_int_early_35+dta_Shp2@data$model_int_early_36+dta_Shp2@data$model_int_early_37+dta_Shp2@data$model_int_early_38+dta_Shp2@data$model_int_early_39+
  dta_Shp2@data$model_int_early_40+dta_Shp2@data$model_int_early_41+dta_Shp2@data$model_int_early_42+dta_Shp2@data$model_int_early_43+dta_Shp2@data$model_int_early_44+dta_Shp2@data$model_int_early_45+dta_Shp2@data$model_int_early_46

dta_Shp2_predict <-subset(dta_Shp2@data, select=c(GridID, predict_NDVI_max_pre))
psm_Long_predict=merge(dta_Shp2_predict, psm_Long, by.x="GridID", by.y="GridID")
psm_Long <- psm_Long_predict




pModelMax_E <- "MaxL_ ~ TrtMnt_demend_y + MeanT_ + MeanP_ + Pop_ + MaxT_ + MaxP_ + MinT_ + MinP_ + 
                predict_NDVI_max_pre_cat*TrtMnt_demend_y + factor(reu_id) + Year"
pModelMax_F <- "MaxL_ ~ TrtMnt_demend_y + MeanT_ + MeanP_ + Pop_ + MaxT_ + MaxP_ + MinT_ + MinP_ + 
                predict_NDVI_max_pre.y*TrtMnt_demend_y + factor(reu_id) + Year"

pModelMax_G <- "MaxL_ ~ TrtMnt_demend_y + MeanT_ + MeanP_ + Pop_ + MaxT_ + MaxP_ + MinT_ + MinP_ + 
                pre_trend_NDVI_max_cat*TrtMnt_demend_y + factor(reu_id) + Year"
pModelMax_H <- "MaxL_ ~ TrtMnt_demend_y + MeanT_ + MeanP_ + Pop_ + MaxT_ + MaxP_ + MinT_ + MinP_ + 
                pre_trend_NDVI_max*TrtMnt_demend_y + factor(reu_id) + Year"

pModelMax_G_fit <- Stage2PSM(pModelMax_G,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMax_H_fit <- Stage2PSM(pModelMax_H ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))


### -------------------

View(psm_Long$MaxL)
temp_TS_median <- fivenum(psm_Long$MaxL[1041:1120])[3]
high_pressure_regions_1995 <- ifelse(psm_Long$MaxL[1041:1120] > temp_TS_median, 1, 0)


high_pressure_regions <- ifelse(psm_Long$reu_id == 118 | psm_Long$reu_id == 142 | 
                                  psm_Long$reu_id == 105 | psm_Long$reu_id == 148 | 
                                  psm_Long$reu_id == 154 | psm_Long$reu_id == 159 | 
                                  psm_Long$reu_id == 160 | psm_Long$reu_id == 161 | 
                                  psm_Long$reu_id == 162 | psm_Long$reu_id == 163 | 
                                  psm_Long$reu_id == 146 | psm_Long$reu_id == 168 | 
                                  psm_Long$reu_id == 151 | psm_Long$reu_id == 157 | 
                                  psm_Long$reu_id == 170 | psm_Long$reu_id == 174 | 
                                  psm_Long$reu_id == 115 | psm_Long$reu_id == 80 | 
                                  psm_Long$reu_id == 147 | psm_Long$reu_id == 74 | 
                                  psm_Long$reu_id == 88 | psm_Long$reu_id == 155 | 
                                  psm_Long$reu_id == 100 | psm_Long$reu_id == 123 | 
                                  psm_Long$reu_id == 172 | psm_Long$reu_id == 133 | 
                                  psm_Long$reu_id == 85 | psm_Long$reu_id == 89 | 
                                  psm_Long$reu_id == 171 | psm_Long$reu_id == 86 | 
                                  psm_Long$reu_id == 91 | psm_Long$reu_id == 175 | 
                                  psm_Long$reu_id == 130 | psm_Long$reu_id == 113 | 
                                  psm_Long$reu_id == 109 | psm_Long$reu_id == 103 | 
                                  psm_Long$reu_id == 134 | psm_Long$reu_id == 179 | 
                                  psm_Long$reu_id == 94 | psm_Long$reu_id == 95, 1, 0)

high_pressure_regions_int <- (high_pressure_regions * psm_Long$TrtMnt_demend_y)

pModelMax_HP <- "MaxL_ ~ TrtMnt_demend_y + MeanT_ + MeanP_ + Pop_ + MaxT_ + MaxP_ + MinT_ + MinP_  + factor(reu_id) + Year + high_pressure_regions + high_pressure_regions_int"
pModelMax_HP_fit <- Stage2PSM(pModelMax_C ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))

#temp_HPR <- ifelse(psm_Long$Year <= 1995 & high_pressure_regions == 1, 1, 0)

## Stargazer Output

stargazer(pModelMax_A_fit $cmreg,pModelMax_B_fit $cmreg,pModelMax_C_fit $cmreg,type="html",align=TRUE,keep=c("TrtMnt","MeanT_","MeanP_","Pop_","MaxT_","MaxP_","MinT_","MinP_","Year"),
          covariate.labels=c("TrtMnt_regend_y","MeanT","MeanP","Pop","MaxT","MaxP","MinT","MinP","Year"),
          omit.stat=c("f","ser"),
          title="Regression Results",
          dep.var.labels=c("Max NDVI")
)

stargazer(pModelMax_E_fit $cmreg,pModelMax_F_fit $cmreg,pModelMax_G_fit$cmreg, pModelMax_H_fit$cmreg,type="html",align=TRUE,
          keep=c("TrtMnt","MeanT_","MeanP_","Pop_","MaxT_","MaxP_","MinT_","MinP_","Year","predict_NDVI_max_pre_cat","TrtMnt_demend_y:predict_NDVI_max_pre_cat", "predict_NDVI_max_pre.y","TrtMnt_demend_y:predict_NDVI_max_pre.y","pre_trend_NDVI_max_cat","TrtMnt_demend_y:pre_trend_NDVI_max_cat","pre_trend_NDVI_max","TrtMnt_demend_y:pre_trend_NDVI_max"),
          omit.stat=c("f","ser"),
          title="Regression Results",
          dep.var.labels=c("Max NDVI")
)

stargazer(pModelMax_G_fit $cmreg,pModelMax_H_fit $cmreg,type="html",align=TRUE,keep=c("TrtMnt","MeanT_","MeanP_","Pop_","MaxT_","MaxP_","MinT_","MinP_","Year","pre_trend_NDVI_max_cat","TrtMnt_demend_y:pre_trend_NDVI_max_cat","pre_trend_NDVI_max","TrtMnt_demend_y:pre_trend_NDVI_max"),
          omit.stat=c("f","ser"),
          title="Regression Results",
          dep.var.labels=c("Max NDVI")
)

