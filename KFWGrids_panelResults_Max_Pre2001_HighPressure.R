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
#Load in Processed Data - produced from script build_panel_ALL151_highpressure.R
#-------------------------------------------------
#-------------------------------------------------

psm_Long<-read.csv("/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/GridDataProcessed/psm_Long_151_highpressure.csv")


#Create post-2004 and post-2008 dummy
psm_Long$Post2004 <- 0
psm_Long$Post2004[psm_Long$Year > 2004] <- 1

psm_Long$Post2008<- 0
psm_Long$Post2008[psm_Long$Year>2008]<-1

## Run Models

pModelMax_A <- "MaxL ~ trtdem + factor(reu_id)"
pModelMax_B <- "MaxL ~ trtdem + Pop+ MeanT + MeanP +MaxT + MaxP + MinT + MinP + factor(reu_id) "
pModelMax_C <- "MaxL ~ trtdem + Pop +MeanT + MeanP +MaxT + MaxP + MinT + MinP + Year + factor(reu_id)"

pModelMax_C2004 <- "MaxL ~ trtdem + Pop + MeanT + MeanP +MaxT + MaxP + MinT + MinP + Year + Post2004*trtdem + factor(reu_id)"
pModelMax_C1_2004 <- "MaxL ~ trtdem + trtenf + Pop + MeanT + MeanP +MaxT + MaxP + MinT + MinP + Year + Post2004*trtenf + factor(reu_id)"

pModelMax_D <- "MaxL ~ trtdem  + Pop + MeanT + MeanP +MaxT + MaxP + MinT + MinP + factor(Year) + factor(reu_id)"
pModelMax_D1 <- "MaxL ~ trtdem + trtenf +Pop +MeanT + MeanP +MaxT + MaxP + MinT + MinP + factor(Year) + factor(reu_id)"
pModelMax_D2004 <- "MaxL ~ trtdem + Pop +MeanT + MeanP +MaxT + MaxP + MinT + MinP + Post2004*trtdem + factor(Year) + factor(reu_id)"
pModelMax_D1_2004 <- "MaxL ~ trtdem + trtenf + Pop +MeanT + MeanP +MaxT + MaxP + MinT + MinP + Post2004*trtenf + factor(Year) + factor(reu_id)"


pModelMax_A_fit <- Stage2PSM(pModelMax_A ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMax_B_fit <- Stage2PSM(pModelMax_B ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))

pModelMax_C_fit <- Stage2PSM(pModelMax_C ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))

pModelMax_C2004_fit <- Stage2PSM(pModelMax_C2004 ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMax_C1_2004_fit <- Stage2PSM(pModelMax_C1_2004 ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))

pModelMax_D_fit <- Stage2PSM(pModelMax_D ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMax_D1_fit <- Stage2PSM(pModelMax_D1 ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))

pModelMax_D2004_fit <- Stage2PSM(pModelMax_D2004 ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMax_D1_2004_fit <- Stage2PSM(pModelMax_D1_2004 ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))

#Create new interaction terms to assist in stargazer formatting later
psm_Long$predict_NDVI_max_pre_cat_int<-psm_Long$predict_NDVI_max_pre_cat*psm_Long$trtdem
psm_Long$predict_NDVI_max_pre_int <- psm_Long$predict_NDVI_max_pre*psm_Long$trtdem
psm_Long$pre_trend_NDVI_max_cat_int <- psm_Long$pre_trend_NDVI_max_cat*psm_Long$trtdem
psm_Long$pre_trend_NDVI_max_int <- psm_Long$pre_trend_NDVI_max*psm_Long$trtdem

pModelMax_E <- "MaxL ~ trtdem + Pop + MeanT + MeanP  + MaxT + MaxP + MinT + MinP + 
                predict_NDVI_max_pre_cat + predict_NDVI_max_pre_cat_int + factor(reu_id) + factor(Year)"

pModelMax_F <- "MaxL ~ trtdem +Pop +MeanT + MeanP +  MaxT + MaxP + MinT + MinP + 
                predict_NDVI_max_pre + predict_NDVI_max_pre_int + factor(reu_id) + factor(Year)"

pModelMax_G <- "MaxL ~ trtdem  + Pop +MeanT + MeanP + MaxT + MaxP + MinT + MinP + 
                pre_trend_NDVI_max_cat + pre_trend_NDVI_max_cat_int + factor(reu_id) + factor(Year)"
pModelMax_H <- "MaxL ~ trtdem +Pop +MeanT + MeanP + MaxT + MaxP + MinT + MinP + 
                pre_trend_NDVI_max + pre_trend_NDVI_max_int + factor(reu_id) + factor(Year)"

pModelMax_E_fit <- Stage2PSM(pModelMax_E,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMax_F_fit <- Stage2PSM(pModelMax_F,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))

pModelMax_G_fit <- Stage2PSM(pModelMax_G,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMax_H_fit <- Stage2PSM(pModelMax_H ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))


## Stargazer Output

#Used for JEEM resubmission
stargazer(pModelMax_A_fit $cmreg,pModelMax_B_fit $cmreg,pModelMax_C_fit $cmreg,pModelMax_D_fit $cmreg,
          pModelMax_D1_fit $cmreg,
          pModelMax_E_fit $cmreg,pModelMax_F_fit $cmreg,pModelMax_G_fit$cmreg, pModelMax_H_fit$cmreg,
          type="html",align=TRUE,
          keep=c("trt","Pop", "MeanT","MeanP","MaxT","MaxP","MinT","MinP",
                 "predict_NDVI_max_pre_cat","predict_NDVI_max_pre_cat_int", 
                 "predict_NDVI_max_pre","predict_NDVI_max_pre_int",
                 "pre_trend_NDVI_max_cat","pre_trend_NDVI_max_cat_int",
                 "pre_trend_NDVI_max","pre_trend_NDVI_max_int","Year"),
          covariate.labels=c("Treatment (Demarcation)","Treatment (Demarcation + Enforcement Support)", "Population", "Mean Temp","Mean Precip","Max Temp","Max Precip","Min Temp","Min Precip",
                              "Predicted NDVI Pre-Trend (Cat)","Predicted NDVI Pre-Trend(Cat)*Treatment","Predicted NDVI Pre-Trend",
                              "Predicted NDVI Pre-Trend * Treatment", "NDVI Pre-Trend (Cat)", "NDVI Pre-Trend(Cat)*Treatment",
                              "NDVI Pre-Trend","NDVI Pre-Trend*Treatment", "Year"),
          order=c("trt","Pop","MeanT","MeanP","MaxT","MaxP","MinT","MinP","NDVI","Year"),
          keep.stat=c("n"),
          add.lines=list(c("Observations","246,007","246,007","246,007","246,007","246,007","246,007","246,007","246,007","246,007"),
               c("Community Fixed Effects?","Yes","Yes","Yes","Yes","Yes","Yes","Yes","Yes","Yes"),
               c("Year Fixed Effects?","No","No","No","Yes","Yes","Yes","Yes","Yes","Yes")),
          title="Regression Results",
          dep.var.labels=c("Max NDVI")
)
#Used for JEEM resubmission
stargazer(pModelMax_C2004_fit $cmreg,pModelMax_D2004_fit $cmreg,pModelMax_C1_2004_fit $cmreg,
          pModelMax_D1_2004_fit $cmreg,
          type="html",align=TRUE,
          keep=c("trt","Pop", "MeanT","MeanP","MaxT","MaxP","MinT","MinP","Year"),
          covariate.labels=c("Treatment (Demarcation)","Treatment (Demarcation + Enforcement Support)", 
                             "Treatment (Demarcation)*Post2004", "Treatment (Dem + Enf)*Post2004",
                             "Population", "Mean Temp","Mean Precip","Max Temp","Max Precip",
                             "Min Temp","Min Precip", "Year"),
          order=c("trt","Pop","MeanT","MeanP","MaxT","MaxP","MinT","MinP","NDVI","Year"),
          keep.stat=c("n"),
          add.lines=list(c("Observations","246,007","246,007","246,007","246,007"),
                         c("Community Fixed Effects?","Yes","Yes","Yes","Yes"),
                         c("Year Fixed Effects?","No","No","Yes","Yes")),
          title="Regression Results",
          dep.var.labels=c("Max NDVI")
)


##STATS TABLE
stargazer(psm_Long, type="html",
          keep=c("MaxL","Slope","Road","Riv","Elevation","terrai_are","Pop","Mean","Min","MaxT",
          "MaxP","pre_trend_NDVI_max","predict_NDVI_max_pre"),
          # covariate.labels=c("NDVI","Slope (degree)","Distance to Road (m)","Distance to River (m)","Elevation (m)",
          #                    "Area (hectares)","Population Density","Mean Temperature","Mean Precipitation",
          #                    "Min Temperature","Min Precipitation","Max Temperature","Max Precipitation",
          #                    "NDVI Pre Trend","Predicted NDVI Pre Trend"),
          omit.summary.stat=c("n"))


#---------------------
#Workspace
#---------------------

## trying to implement lag function

psm_Long_lag <- TimeSeriesLag(psm_Long,"Year","GridID",1,"MaxL_","MaxL_lag",1983,2010)
psm_Long_lag_test <- psm_Long_lag[psm_Long_lag["GridID"]==319588,]

psm_Long_85 <- psm_Long[psm_Long["reu_id"]==85,]
psm_Long_143 <- psm_Long[psm_Long["reu_id"]==143,]

pModelMax_I <- "MaxL_ ~ MaxL_lag"
pModelMax_J <- "MaxL_ ~ MaxL_lag + factor(reu_id)"
pModelMax_K <- "MaxL_ ~ pre_trend_NDVI_max + factor(reu_id)"

pModelMax_I_fit <- Stage2PSM(pModelMax_I,psm_Long_lag,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMax_J_fit <- Stage2PSM(pModelMax_J,psm_Long_lag,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMax_K_fit <- Stage2PSM(pModelMax_K,psm_Long_lag,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))

## trying to weight by size of community for summary stats

psm_Long$commwt <- 1/psm_Long$terrai_are
summary(psm_Long$commwt)

library(SDMTools)
wt.mean(psm_Long$terrai_are,psm_Long$commwt)
wt.sd(psm_Long$terrai,psm_Long$commwt)

test<-psm_Long[!is.na(psm_Long$demend_y),]

##creating correlation coefficients with outcome/control vars and year of demarcation for summary stats table
#produced for JEEM second revision

cor(psm_Long$Slope,psm_Long$demend_y)

##learning how to output tables
test<-function(mean) {
  tabl<-mean(psm_Long$terrai_are)
  return(tabl)
  }

apply(psmenf,2,function(x) length (x[x<0]))
y<-psm_Long$demend_y
test<-psm_Long[,1:4]
tbl<-sapply(test,function(x) cor(x,y))

tbl<-sapply(test,function(x) mean(x))
xtable(tbl)

design.matrix <- model.matrix(mean(x), data = test)
xtable(design.matrix, digits = 0)

dat <- psm_Long[1:3, 1:4]
z <- xtable(bubba)
ztab<-xtable(z)
print.xtable(ztab,type="html")

bubba<-data.frame(row.names<-names(tbl),matrix(unlist(tbl)))
bub<-data.frame(row.names<-names(blah),matrix(unlist(blah)))

blah<-sapply(dat,function(x) mean(x))
bubba_large<-cbind(bubba,bub)
bubba_large<-merge(bubba,bub,by.x=bubba$row.names....names.tbl.,by.y=bub$row.names....names.blah.)

sapply(psmenf, function(x) cor(psmenf$demend_y,psmenf$terrai_are))
dtf <- sapply(psmenf, each(min, max, mean, sd, var, median, IQR))
dtf

#------
#Archived
#------

pModelMax_C_reuid <- "MaxL_ ~ trtdem + Pop_ + MeanT_ + MeanP_ +MaxT_ + MaxP_ + MinT_ + MinP_ + Year*factor(reu_id) + factor(reu_id)"
pModelMax_C_reuidenf <- "MaxL_ ~ trtdem + trtenf+ Pop_ + MeanT_ + MeanP_ +MaxT_ + MaxP_ + MinT_ + MinP_ + Year*factor(reu_id) + factor(reu_id)"
pModelMax_C_reuid_fe <- "MaxL_ ~ trtdem + Pop_ + MeanT_ + MeanP_ +MaxT_ + MaxP_ + MinT_ + MinP_ + Year*factor(reu_id) + factor(Year) + factor(reu_id)"
pModelMax_C_reuidenf_fe <- "MaxL_ ~ trtdem + trtenf+ Pop_ + MeanT_ + MeanP_ +MaxT_ + MaxP_ + MinT_ + MinP_ + Year*factor(reu_id) + factor(Year) + factor(reu_id)"

pModelMax_C_reuid_fit <- Stage2PSM(pModelMax_C_reuid ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMax_C_reuidenf_fit <- Stage2PSM(pModelMax_C_reuidenf ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMax_C_reuid_fe_fit <- Stage2PSM(pModelMax_C_reuid_fe ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMax_C_reuidenf_fe_fit <- Stage2PSM(pModelMax_C_reuidenf_fe ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))

