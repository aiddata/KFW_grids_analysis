
#-------------------------
# Build panel dataset
# All 151 PPTAL communities, no trimming
# Includes high pressure covariates and pre-trends
#------------------------

library(maptools)
library(reshape)
library(splitstackshape)
library(ggplot2)

#-----------------------------------------------
#Read in merged shapefile with pre-trends and covariates and rename pre-trends
#-----------------------------------------------

shpfile = "/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/GridDataProcessed/OhFive_gridanalysis_inputs_pretrends_ALL151.shp"
dta_Shp = readShapePoly(shpfile)

names(dta_Shp@data)[names(dta_Shp@data)=="pt_NDVI_ma"] <- "pre_trend_NDVI_max"
names(dta_Shp@data)[names(dta_Shp@data)=="pt_temp_me"] <- "pre_trend_temp_mean"
names(dta_Shp@data)[names(dta_Shp@data)=="pt_temp_ma"] <- "pre_trend_temp_max"
names(dta_Shp@data)[names(dta_Shp@data)=="pt_temp_mi"] <- "pre_trend_temp_min"
names(dta_Shp@data)[names(dta_Shp@data)=="pt_prec_me"] <- "pre_trend_precip_mean"
names(dta_Shp@data)[names(dta_Shp@data)=="pt_prec_ma"] <- "pre_trend_precip_max"
names(dta_Shp@data)[names(dta_Shp@data)=="pt_prec_mi"] <- "pre_trend_precip_min"

names(dta_Shp@data)[names(dta_Shp@data)=="pt_ntl"] <- "pre_trend_ntl"
names(dta_Shp@data)[names(dta_Shp@data)=="pt_cv"] <- "pre_trend_cv"
names(dta_Shp@data)[names(dta_Shp@data)=="pt_cy"] <- "pre_trend_cy"
names(dta_Shp@data)[names(dta_Shp@data)=="pt_rv"] <- "pre_trend_rv"
names(dta_Shp@data)[names(dta_Shp@data)=="pt_ry"] <- "pre_trend_ry"
names(dta_Shp@data)[names(dta_Shp@data)=="pt_sov"] <- "pre_trend_sov"
names(dta_Shp@data)[names(dta_Shp@data)=="pt_soy"] <- "pre_trend_soy"
names(dta_Shp@data)[names(dta_Shp@data)=="pt_suv"] <- "pre_trend_suv"
names(dta_Shp@data)[names(dta_Shp@data)=="pt_suy"] <- "pre_trend_suy"
names(dta_Shp@data)[names(dta_Shp@data)=="pt_wv"] <- "pre_trend_wv"


#-----------------------------------------------
#Run high pressure model
#-----------------------------------------------

#Creating averages for 1982-1995 pre-period, to use in high pressure model regressions

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

#High Pressure Model
#experimented with other options for high pressure model, now moved to "scratch" section

HPModel = lm(pre_trend_NDVI_max ~ Pop_1995 + pre_trend_temp_mean + pre_trend_temp_min + 
               pre_trend_temp_max + pre_trend_precip_min + pre_trend_precip_mean + pre_trend_precip_max + 
               pre_trend_ntl + Slope + Elevation + Riv_Dist + urbtravtim + Road_dist + 
               log_dist + mine_dist + fedcon_dis + stcon_dist + rail_dist + 
               pre_trend_cv + pre_trend_cy + pre_trend_ry + pre_trend_sov  +
               pre_trend_suv + pre_trend_suy + 
               MeanT_82_95 + MinT_82_95 + MaxT_82_95 + MinP_82_95 + MeanP_82_95 + MaxP_82_95 + ntl_92_95 +
               cv_94_95 + cy_91_95 + ry_90_95 + suv_94_95 + suy_91_95,
             data=dta_Shp@data)

#Running the model with cmreg
HPModel$Id <- cluster.vcov(HPModel,c(dta_Shp@data$Id))
CMREG <- coeftest(HPModel, HPModel$Id)
print(CMREG)
summary(HPModel)

#------------------------------------------------
#Identify high pressure cells
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

#Creating categorical measure of pre_trend_NDVI_max to interact with treatment binary
#(where lower than median value of pre_trend_NDVI_max is considered high pressure)

#create categorical measure using pre-trends
pretrend_NDVI_median<-fivenum(dta_Shp$pre_trend_NDVI_max)[3]
dta_Shp$pre_trend_NDVI_max_cat <- NA
dta_Shp$pre_trend_NDVI_max_cat <-ifelse(dta_Shp$pre_trend_NDVI_max<pretrend_NDVI_median,1,0)

#create categorical measure using predicted pre-trend from high pressure model
predict_NDVI_median<-fivenum(dta_Shp$predict_NDVI_max_pre)[3]
dta_Shp$predict_NDVI_max_pre_cat <- NA
dta_Shp$predict_NDVI_max_pre_cat <-ifelse(dta_Shp$predict_NDVI_max_pre<predict_NDVI_median,1,0)

#-------------------------------------------------
#Convert from a wide-form dataset for the Cross-sectional 
#to a long-form dataset for the panel model.
#-------------------------------------------------

#Drop data for unused years (after 2010) in order to allow "reshape" to work
kfw_grid_reshape <- dta_Shp
kfw_grid_reshape1<-kfw_grid_reshape[,-grep("(2011)",names(kfw_grid_reshape))]
kfw_grid_reshape2<-kfw_grid_reshape1[,-grep("(2012)",names(kfw_grid_reshape1))]
kfw_grid_reshape3<-kfw_grid_reshape2[,-grep("(2013)",names(kfw_grid_reshape2))]
kfw_grid_reshape4<-kfw_grid_reshape3[,-grep("(2014)",names(kfw_grid_reshape3))]
kfw_grid_reshape5<-kfw_grid_reshape4[,-grep("(1981)",names(kfw_grid_reshape4))]

#Eliminate data for vars we don't need
kfw_grid_reshape_a <- kfw_grid_reshape5
kfw_grid_reshape_b<-kfw_grid_reshape_a[,-grep("(cv)",names(kfw_grid_reshape_a))]
kfw_grid_reshape_c<-kfw_grid_reshape_b[,-grep("(cy)",names(kfw_grid_reshape_b))]
kfw_grid_reshape_d<-kfw_grid_reshape_c[,-grep("(rv)",names(kfw_grid_reshape_c))]
kfw_grid_reshape_e<-kfw_grid_reshape_d[,-grep("(ry)",names(kfw_grid_reshape_d))]
kfw_grid_reshape_f<-kfw_grid_reshape_e[,-grep("(sov)",names(kfw_grid_reshape_e))]
kfw_grid_reshape_g<-kfw_grid_reshape_f[,-grep("(soy)",names(kfw_grid_reshape_f))]
kfw_grid_reshape_h<-kfw_grid_reshape_g[,-grep("(suv)",names(kfw_grid_reshape_g))]
kfw_grid_reshape_i<-kfw_grid_reshape_h[,-grep("(suy)",names(kfw_grid_reshape_h))]
kfw_grid_reshape_j<-kfw_grid_reshape_i[,-grep("(wv)",names(kfw_grid_reshape_i))]
kfw_grid_reshape_k<-kfw_grid_reshape_j[,-grep("(model_int)",names(kfw_grid_reshape_j))]
kfw_grid_reshape_l<-kfw_grid_reshape_k[,-grep("(_95)",names(kfw_grid_reshape_k))]

# *** Recreate kfw_grid_wide if drop/add new variables using kfw_grid_reshape***

kfw_grid_wide<-kfw_grid_reshape_l

#Order variables chronologically to allow reshape to work properly
kfw_grid_wide<-kfw_grid_wide[,order(names(kfw_grid_wide))]

#Identify variables where values will change yearly in panel dataset
MeanT<-grep("MeanT_",names(kfw_grid_wide))
MeanP<-grep("MeanP_",names(kfw_grid_wide))
MinT<-grep("MinT_",names(kfw_grid_wide))
MaxT<-grep("MaxT_",names(kfw_grid_wide))
MinP<-grep("MinP_",names(kfw_grid_wide))
MaxP<-grep("MaxP_",names(kfw_grid_wide))
MaxL<-grep("MaxL_",names(kfw_grid_wide))
Pop<-grep("Pop_",names(kfw_grid_wide))

#Convert to long form panel dataset
all_reshape <- c(MeanT,MeanP,MaxT,MaxP,MinP,MinT,MaxL,Pop)
psm_Long <- reshape(kfw_grid_wide@data, varying=all_reshape, direction="long",idvar="GridID",sep="_",timevar="Year")

#Create years to demarcation
psm_Long$yrtodem <- NA
psm_Long$yrtodem=psm_Long$Year - psm_Long$demend_y

#Create demarcation treatment variable, using demend_y
#0 in years prior to demarcation, turns to 1 in year of demarcation
psmtest3 <- psm_Long
psmtest3$trtdem <- 0
psmtest3$trtdem[which(psmtest3$Year<psmtest3$demend_y)]<-0
psmtest3$trtdem[which(psmtest3$Year>=psmtest3$demend_y)]<-1

psm_Long <- psmtest3

#Correct communities with enforcement start date prior to demarcation
#Only 1 community, reu_id==84
psmtest5 <- psm_Long
psmtest5$enfdiff= psmtest5$enforce_st - psmtest5$demend_y
summary(psmtest5$enfdiff)
psmenf <- subset(psmtest5, psmtest5$enfdiff<0)
table(psmenf$reu_id)

psmtest5$enforce_st[which(psmtest5$reu_id==84)]<-2007

#Create enforcement treatment var
#0 before enforcement support begins, turns to 1 in the year that enforcement support starts
psmtest5$trtenf <- 0
psmtest5$trtenf[which(psmtest5$Year>=psmtest5$enforce_st)]<-1
table(psmtest5$trtenf,psmtest5$trtdem)
#Should only apply to communities that have been demarcated
# For non-demarcated communities, set all values back to 0
psmtest5$trtenf[is.na(psmtest5$demend_y)]<-0
table(psmtest5$trtenf,psmtest5$trtdem)

psm_Long <- psmtest5

write.csv(psm_Long,file="/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/GridDataProcessed/psm_Long_151_highpressure.csv")
psm_Long<-read.csv("/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/GridDataProcessed/psm_Long_151_highpressure.csv")

##Scratch

#Use to double check that panel dataset was built correctly
sub29<-kfw_grid_wide[kfw_grid_wide@data$GridID==225929,]
View(sub29@data)
View(sub29@data[,(100:200)])
View(sub29@data[,(200:296)])

#*** Alternative High Pressure Models ***

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
#everything
HPModel = lm(MaxL_levchange_95_82 ~ Pop_1995 + pre_trend_temp_mean + pre_trend_temp_min +
pre_trend_temp_max + pre_trend_precip_min + pre_trend_precip_mean + pre_trend_precip_max + 
  pre_trend_ntl + Slope + Elevation + Riv_Dist + urbtravtim + Road_dist + 
  log_dist + mine_dist + fedcon_dis + stcon_dist + rail_dist + 
  pre_trend_cv + pre_trend_cy + pre_trend_rv + pre_trend_ry + pre_trend_sov + pre_trend_soy +
  pre_trend_suv + pre_trend_suy + pre_trend_wv + 
  MeanT_82_95 + MinT_82_95 + MaxT_82_95 + MinP_82_95 + MeanP_82_95 + MaxP_82_95 +
  ntl_92_95 + cv_94_95 + cy_91_95 + rv_94_95 + ry_90_95 + sov_94_95 + soy_91_95 + suv_94_95 + 
  suy_91_95 + wv_94_95,data=dta_Shp@data)

