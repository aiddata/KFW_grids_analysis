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


#------
#Build Dataset

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

varList=c("MaxL_")
psm_Long_Untrimmed <- BuildTimeSeries(dta=dta_Shp,idField="GridID",varList_pre=varList,1982,2010,colYears=c("demend_y","enforce_st"),
                                      interpYears=c("Slope","Road_dist","Riv_Dist","UF","Elevation","terrai_are","Pop_","MeanT_","MeanP_","MaxT_",
                                                    "MaxP_","MinP_","MinT_", "reu_id", "Id" ))
psm_Long_Untrimmed$Year <- as.numeric(psm_Long_Untrimmed$Year)

#merge in demend_y and enforce_st to create correct treatment variables
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

#write.csv(psm_Long_Untrimmed,file="/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/GridDataProcessed/psm_Long_Untrimmed.csv")

#Done with dataset build
#------------

#Read in already created dataset
psm_Long_Untrimmed <- read.csv("/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/GridDataProcessed/psm_Long_Untrimmed.csv")


#Create subset that only includes years within -5 and +5 years of demarcation
# psm_Long_5yr <- psm_Long
# test <- psm_Long_5yr[psm_Long_5yr$yrtodem>=-5,]
# test1 <- test[test$yrtodem<=5,]
# psm_Long <- test
#write.csv(psm_Long_5yr,file="/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/GridDataProcessed/psm_Long_Untrimmed_5Yr.csv")

#Create subset that only includes reu_ids for the pairs made from 1st stage PSM (at community level) with 
# only the predictors that were significant (slope, elevation, road distance, pre-trend in min annual temp) at the COMMUNITY level

#psm_pairs (37 pairs, 74 communities)

pairs_id <- c(131,114,118,142,117, 121, 105, 148,  
  93, 107, 152, 150, 154, 112, 158, 159, 
  160, 161, 162, 163, 164, 146, 110, 180, 
  168, 151, 157,173, 176, 115,  80,  92,74,
  119, 132,  88, 128, 155, 129, 156, 100, 123,
  106, 172,  87,  78,  73, 122, 169, 144, 133, 111,  85,  89,
  79, 86,  91, 175,  81,  82, 125, 126, 141,  96, 109, 
  103, 143, 137, 135, 136, 134, 179, 178,  95)

#psmRes (no pairs, 92 communities out of 106, 14 dropped for common support)

res_id <- c(80, 131, 120, 114, 118,  96,  92, 109, 142,  74, 143, 167, 116, 117, 119, 132,  88, 128, 155, 129, 156, 100, 121, 123, 102,  98, 135,
103, 130, 105, 153, 127, 104, 134, 125, 126, 141, 124, 106, 148,  93, 172, 107,  87,  94,  78,  73, 152, 150, 154, 112, 138, 178, 136,
122, 169, 144, 158, 159, 145, 133, 160, 161, 162, 163, 164, 165, 146, 113, 137, 110, 180,  95, 111, 168, 179, 151,  85, 157,  89, 170,
79,  86,  91, 173, 174, 175, 176, 177,  81,  82, 115)

#subsetting the data
psm_Long_sigpsm <- psm_Long_Untrimmed
psm_Long_sigpsm1 <- psm_Long_sigpsm[psm_Long_sigpsm$reu_id %in% res_id,]
#write.csv(psm_Long_sigpsm1, file="/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/GridDataProcessed/psm_Long_Untrimmed_sigpsm.csv")

psm_Long <- psm_Long_sigpsm1

## Run Models

pModelMax_A <- "MaxL_ ~ trtdem + trtenf + factor(reu_id)"
pModelMax_B <- "MaxL_ ~ trtdem + trtenf  + Pop_ + MeanT_ + MeanP_ + MaxT_ + MaxP_ + MinT_ + MinP_  + factor(reu_id) "
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
          keep=c("trt","Pop","Mean","Max","Min","Year"),
          covariate.labels=c("Treatment (Demarcation)","Treatment (Demarcation + Enforcement Support)","Population","Mean Temp",
                           "Mean Precip","Max Temp","Max Precip","Min Temp","Min Precip","Year"),
          omit.stat=c("f","ser"),
          add.lines=list(c("Observations","404,405","404,405","404,405","404,405","404,405"),
                         c("Community Fixed Effects?","Yes","Yes","Yes","Yes","Yes"),
                         c("Year Fixed Effects?","No","No","No","Yes","Yes")),
          title="Regression Results",
          dep.var.labels=c("Max NDVI"))


##Workspace##

reg=lm(MaxL_ ~ factor(Year), data=psm_Long)
resid <- residuals(reg)
summary(resid)
plot(resid)

ViewTimeSeries(dta=psm_Long,IDfield="reu_id",TrtField="TrtBin",idPre="MaxL_")

ggplot(data = psm_Pairs, aes(x=variable, y=value, group="reu_id",colour=factor("TrtBin")),  
       geom_line(size=.5, linetype=3), 
       stat_summary(fun.y=mean,aes(x=variable, y=value, group="TrtBin",colour=factor("TrtBin")),data=psm_Long,geom='line',size=1.5),
       theme(axis.text.x=element_text(angle=90,hjust=1)))









