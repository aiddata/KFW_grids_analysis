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

#-------------------------------------------------
#-------------------------------------------------
#Pre-processing to create cross-sectional variable summaries
#-------------------------------------------------
#-------------------------------------------------
#Calculate NDVI Trends
#dta_Shp$pre_trend_NDVI_mean <- timeRangeTrend(dta_Shp,"MeanL_[0-9][0-9][0-9][0-9]",1982,1995,"SP_ID")
#dta_Shp$pre_trend_NDVI_max <- timeRangeTrend(dta_Shp,"MaxL_[0-9][0-9][0-9][0-9]",1982,1995,"GridID")
#dta_Shp$NDVIslope_95_10 <- timeRangeTrend(dta_Shp,"MaxL_[0-9][0-9][0-9][0-9]",1995,2010,"SP_ID")
#dta_Shp@data["NDVIslopeChange_95_10"] <- dta_Shp$MeanL_2010 - dta_Shp$MeanL_1995

#NDVI Trends for 1995-2001
#dta_Shp$post_trend_NDVI_95_01 <- timeRangeTrend(dta_Shp,"MeanL_[0-9][0-9][0-9][0-9]",1995,2001,"SP_ID")
#dta_Shp@data["NDVIslopeChange_95_01"] <- dta_Shp$MeanL_2001 - dta_Shp$MeanL_1995
#NDVI Trends for 2001-2010
#dta_Shp$post_trend_NDVI_01_10 <- timeRangeTrend(dta_Shp,"MeanL_[0-9][0-9][0-9][0-9]",2001,2010,"SP_ID")
#dta_Shp@data["NDVIslopeChange_01_10"] <- dta_Shp$MeanL_2010 - dta_Shp$MeanL_2001
#dta_Shp@data["NDVIslopeChange_01_10"] <- dta_Shp@data["post_trend_NDVI_01_10"] - dta_Shp@data["pre_trend_NDVI_max"]

#Calculate Temp and Precip Pre and Post Trends
# dta_Shp$pre_trend_temp_mean <- timeRangeTrend(dta_Shp,"MeanT_[0-9][0-9][0-9][0-9]",1982,1995,"SP_ID")
# dta_Shp$pre_trend_temp_max <- timeRangeTrend(dta_Shp,"MaxT_[0-9][0-9][0-9][0-9]",1982,1995,"SP_ID")
# dta_Shp$pre_trend_temp_min <- timeRangeTrend(dta_Shp,"MinT_[0-9][0-9][0-9][0-9]",1982,1995,"SP_ID")

#dta_Shp$post_trend_temp_mean <- timeRangeTrend(dta_Shp,"MeanT_[0-9][0-9][0-9][0-9]",1995,2010,"SP_ID")
#dta_Shp$post_trend_temp_max <- timeRangeTrend(dta_Shp,"MaxT_[0-9][0-9][0-9][0-9]",1995,2010,"SP_ID")
#dta_Shp$post_trend_temp_min <- timeRangeTrend(dta_Shp,"MinT_[0-9][0-9][0-9][0-9]",1995,2010,"SP_ID")
#dta_Shp$post_trend_temp_95_01 <- timeRangeTrend(dta_Shp,"MeanT_[0-9][0-9][0-9][0-9]",1995,2001,"SP_ID")
#dta_Shp$post_trend_temp_01_10 <- timeRangeTrend(dta_Shp,"MeanT_[0-9][0-9][0-9][0-9]",2001,2010,"SP_ID")

# dta_Shp$pre_trend_precip_mean <- timeRangeTrend(dta_Shp,"MeanP_[0-9][0-9][0-9][0-9]",1982,1995,"SP_ID")
# dta_Shp$pre_trend_precip_max <- timeRangeTrend(dta_Shp,"MaxP_[0-9][0-9][0-9][0-9]",1982,1995,"SP_ID")
# dta_Shp$pre_trend_precip_min <- timeRangeTrend(dta_Shp,"MinP_[0-9][0-9][0-9][0-9]",1982,1995,"SP_ID")

# dta_Shp$post_trend_precip_mean <- timeRangeTrend(dta_Shp,"MeanP_[0-9][0-9][0-9][0-9]",1995,2010,"SP_ID")
# dta_Shp$post_trend_precip_max <- timeRangeTrend(dta_Shp,"MaxP_[0-9][0-9][0-9][0-9]",1995,2010,"SP_ID")
# dta_Shp$post_trend_precip_min <- timeRangeTrend(dta_Shp,"MinP_[0-9][0-9][0-9][0-9]",1995,2010,"SP_ID")
# dta_Shp$post_trend_precip_95_01 <- timeRangeTrend(dta_Shp,"MeanP_[0-9][0-9][0-9][0-9]",1995,2001,"SP_ID")
# dta_Shp$post_trend_precip_01_10 <- timeRangeTrend(dta_Shp,"MeanP_[0-9][0-9][0-9][0-9]",2001,2010,"SP_ID")


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

#Remove units that did not ever receive any treatment (within-sample test)
dta_Shp@data$NA_check <- 0
dta_Shp@data$NA_check[is.na(dta_Shp@data$demend_y)] <- 1
int_Shp <- dta_Shp[dta_Shp@data$NA_check != 1,]
dta_Shp <- int_Shp



#-------------------------------------------------
#-------------------------------------------------
#Define and run the first-stage of the PSM, calculating propensity scores
#-------------------------------------------------
#-------------------------------------------------
psmModel <-  "TrtBin ~ terrai_are + Pop_1990 + MeanT_1995 + pre_trend_temp_mean + pre_trend_temp_min + 
pre_trend_temp_max + MeanP_1995 + pre_trend_precip_min + pre_trend_NDVI_max + Slope + Elevation + MaxL_1995 + Riv_Dist + Road_dist +
pre_trend_precip_mean + pre_trend_precip_max"

psmRes <- SCI::SpatialCausalPSM(dta_Shp,mtd="logit",psmModel,drop="support",visual=FALSE)


#-------------------------------------------------
#-------------------------------------------------
#Based on the Propensity Score Matches, pair comprable treatment and control units.
#-------------------------------------------------
#-------------------------------------------------
drop_set<- c(drop_unmatched=TRUE,drop_method="None",drop_thresh=0.25)
psm_Pairs <- SAT(dta = psmRes$data, mtd = "fastNN",constraints=c(groups="UF"),psm_eq = psmModel, ids = "GridID", drop_opts = drop_set, visual="TRUE", TrtBinColName="TrtBin")
#c(groups=c("UF"),distance=NULL)
trttable <- table (psm_Pairs@data$TrtBin)
View(trttable)


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
                               interpYears=c("Slope","Road_dist","Riv_Dist","UF","Elevation","terrai_are","Pop_","MeanT_","MeanP_","MaxT_","MaxP_","MinP_","MinT_","ntl_", "fedcon_dis", "stcon_dist", "log_dist", "mine_dist", "rail_dist", "reu_id", "Id" ))
psm_Long_HP$Year <- as.numeric(psm_Long_HP$Year)

write.csv(psm_Long_HP,file="/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/GridDataProcessed/psm_Long_HP.csv")



dta_Shp_id <-subset(dta_Shp, select=c(GridID, reu_id, Id))
psm_Long_reuid=merge(dta_Shp_id@data, psm_Long, by.x="GridID", by.y="GridID")
psm_Long <- psm_Long_reuid

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


#**Note pre_trend_NDVI_max works from NDVI values that are 10,000 times the actual values, so slope is inflated as well

#Creating climate averages for 1982-1995

dta_Shp2$MeanT_82_95 <- timeRangeAvg(dta_Shp2@data, "MeanT_",1982,1995)
dta_Shp2$MinT_82_95 <- timeRangeAvg(dta_Shp2@data, "MinT_",1982,1995)
dta_Shp2$MaxT_82_95 <- timeRangeAvg(dta_Shp2@data, "MaxT_",1982,1995)
dta_Shp2$MeanP_82_95 <- timeRangeAvg(dta_Shp2@data, "MeanP_",1982,1995)
dta_Shp2$MinP_82_95 <- timeRangeAvg(dta_Shp2@data, "MinP_",1982,1995)
dta_Shp2$MaxP_82_95 <- timeRangeAvg(dta_Shp2@data, "MaxP_",1982,1995)

#Creating var for level change in NDVI, 1995-1982
dta_Shp2$MaxL_levchange_95_82 <- dta_Shp2$MaxL_1995 - dta_Shp2$MaxL_1982

#Subsetting data to create data frame with observations that have negative pre_trend_NDVI_max
dta_Shp2 <- dta_Shp
dta_Shp3=dta_Shp2[dta_Shp2$pre_trend_NDVI_max<0,]
#Creating binary measure where negative pre_trend_NDVI_max is negative (to reflect deforestation)
dta_Shp2$BinNDVI=0
dta_Shp2$BinNDVI[dta_Shp2$pre_trend_NDVI_max<0]=1

dta = dta_Shp2@data

#logit with binary where projects with neg pre-trends=1, all else=0
HPModel = logit(BinNDVI ~ terrai_are + Pop_1995 + pre_trend_temp_mean + pre_trend_temp_min + 
                  pre_trend_temp_max + pre_trend_precip_min + pre_trend_precip_mean + pre_trend_precip_max + Slope + 
                  Elevation + Riv_Dist + ntl_1992 + ntl_1993 + ntl_1994 + ntl_1995 + urbtravtim, data=dta_Shp2@data)
#lm with only grid cells that have neg pre-trends, only 375 obs
HPModel = lm(pre_trend_NDVI_max ~ terrai_are + Pop_1995 + pre_trend_temp_mean + pre_trend_temp_min + 
               pre_trend_temp_max + pre_trend_precip_min + pre_trend_precip_mean + pre_trend_precip_max + Slope + 
               Elevation + Riv_Dist + ntl_1992 + ntl_1993 + ntl_1994 + ntl_1995 + urbtravtim, data=dta_Shp2@data)
#lm with all grid cells, i.e. pos and neg pre-trends
HPModel = lm(pre_trend_NDVI_max ~ terrai_are + Pop_1995 + pre_trend_temp_mean + pre_trend_temp_min + 
               pre_trend_temp_max + pre_trend_precip_min + pre_trend_precip_mean + pre_trend_precip_max + Slope + 
               Elevation + Riv_Dist + ntl_1992 + ntl_1993 + ntl_1994 + ntl_1995 + urbtravtim, data=dta_Shp2@data)
#lm with all grid cells, i.e. pos and neg pre-trends, and all covars from any time period
HPModel = lm(pre_trend_NDVI_max ~ terrai_are + Pop_1995 + pre_trend_temp_mean + pre_trend_temp_min + 
               pre_trend_temp_max + pre_trend_precip_min + pre_trend_precip_mean + pre_trend_precip_max + Slope + 
               Elevation + Riv_Dist + ntl_1992 + ntl_1993 + ntl_1994 + ntl_1995 + urbtravtim + Road_dist + log_dist +
               mine_dist + fedcon_dis + stcon_dist + rail_dist, data=dta_Shp3@data)
HPModel = glm(BinNDVI ~ terrai_are + Pop_1995 + pre_trend_temp_mean + pre_trend_temp_min + 
                pre_trend_temp_max + pre_trend_precip_min + pre_trend_precip_mean + pre_trend_precip_max + Slope + 
                Elevation + Riv_Dist + ntl_1992 + ntl_1993 + ntl_1994 + ntl_1995 + urbtravtim + Road_dist + log_dist +
                mine_dist + fedcon_dis + stcon_dist + rail_dist, family=binomial(logit), data=dta_Shp2@data)
HPModel = glm(BinNDVI ~ terrai_are + Pop_1995 + pre_trend_temp_mean + pre_trend_temp_min + 
                pre_trend_temp_max + pre_trend_precip_min + pre_trend_precip_mean + pre_trend_precip_max + Slope + 
                Elevation + Riv_Dist + ntl_1992 + ntl_1993 + ntl_1994 + ntl_1995 + urbtravtim + Road_dist + log_dist +
                mine_dist + fedcon_dis + stcon_dist + rail_dist + cv1995 + cv1994 + cy1995 + cy1994 + cy1993 + cy1992 +
                cy1991 + rv1995 + rv1994 + ry1995 + ry1994 + ry1993 + ry1992 + ry1991 + ry1990 + sov1995 + sov1994 +
                soy1995 + soy1994 + soy1993 + soy1992 + soy1991 + suv1995 + suv1994 + suy1995 + suy1994 + suy1993 +
                suy1992 + suy1991 + wv1995 + wv1994, family=binomial(logit), data=dta_Shp2@data)
HPModel = lm(pre_trend_NDVI_max ~ terrai_are + Pop_1995 + pre_trend_temp_mean + pre_trend_temp_min + 
               pre_trend_temp_max + pre_trend_precip_min + pre_trend_precip_mean + pre_trend_precip_max + Slope + 
               Elevation + Riv_Dist + ntl_1992 + ntl_1993 + ntl_1994 + ntl_1995 + urbtravtim + Road_dist + log_dist +
               mine_dist + fedcon_dis + stcon_dist + rail_dist + cv1995 + cv1994 + cy1995 + cy1994 + cy1993 + cy1992 +
               cy1991 + rv1995 + rv1994 + ry1995 + ry1994 + ry1993 + ry1992 + ry1991 + ry1990 + sov1995 + sov1994 +
               soy1995 + soy1994 + soy1993 + soy1992 + soy1991 + suv1995 + suv1994 + suy1995 + suy1994 + suy1993 +
               suy1992 + suy1991 + wv1995 + wv1994 , data=dta_Shp2@data)
HPModel = lm(pre_trend_NDVI_max ~ Pop_1990 + Pop_1995 + MeanT_82_95 + MinT_82_95 + 
               MaxT_82_95 + MinP_82_95 + MeanP_82_95 + MaxP_82_95 + Slope + 
               Elevation + Riv_Dist + ntl_1992 + ntl_1993 + ntl_1994 + ntl_1995 + urbtravtim + Road_dist + log_dist +
               mine_dist + fedcon_dis + stcon_dist + rail_dist + cv1995 + cv1994 + cy1995 + cy1994 + cy1993 + cy1992 +
               cy1991 + rv1995 + rv1994 + ry1995 + ry1994 + ry1993 + ry1992 + ry1991 + ry1990 + sov1995 + suv1995 + suv1994 + suy1995 + suy1994 + suy1993 +
               suy1992 + suy1991, data=dta_Shp2@data)
HPModel = lm(MaxL_levchange_95_82 ~ Pop_1990 + Pop_1995 + MeanT_82_95 + MinT_82_95 + 
               MaxT_82_95 + MinP_82_95 + MeanP_82_95 + MaxP_82_95 + Slope + 
               Elevation + Riv_Dist + ntl_1992 + ntl_1993 + ntl_1994 + ntl_1995 + urbtravtim + Road_dist + log_dist +
               mine_dist + fedcon_dis + stcon_dist + rail_dist + cv1995 + cv1994 + cy1995 + cy1994 + cy1993 + cy1992 +
               cy1991 + rv1995 + rv1994 + ry1995 + ry1994 + ry1993 + ry1992 + ry1991 + ry1990 + sov1995 + suv1995 + 
               suv1994 + suy1995 + suy1994 + suy1993 +
               suy1992 + suy1991, data=dta_Shp2@data)


#Running the model with cmreg
HPModel$Id <- cluster.vcov(HPModel,c(dta$Id))
CMREG <- coeftest(HPModel, HPModel$Id)
print(CMREG)
summary(HPModel)
summary(HPModel)$r.squared

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

psm_Long$predict_NDVI_max_pre_cat <- NA
psm_Long$predict_NDVI_max_pre_cat[which(psm_Long$predict_NDVI_max_pre)]

pModelMax_E <- "MaxL_ ~ TrtMnt_demend_y + MeanT_ + MeanP_ + Pop_ + MaxT_ + MaxP_ + MinT_ + MinP_ + 
                predict_NDVI_max_pre*TrtMnt_demend_y + factor(reu_id) + Year"

pModelMax_E_fit <- Stage2PSM(pModelMax_E ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))


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

