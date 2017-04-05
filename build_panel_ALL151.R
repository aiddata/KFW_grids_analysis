
#-------------------------
# Build panel dataset for analysis
# All 151 PPTAL communities
# Simple version -- no high pressure covariates or pre-trends
#------------------------

library(maptools)
library(reshape)
library(splitstackshape)
library(ggplot2)



kfw_grid100 = readShapePoly("OhFive_gridanalysis_inputs_ALL151.shp")


#Drop data for unused years (after 2010) in order to allow "reshape" to work
kfw_grid_reshape <- kfw_grid20
kfw_grid_reshape1<-kfw_grid_reshape[,-grep("(2011)",names(kfw_grid_reshape))]
kfw_grid_reshape2<-kfw_grid_reshape1[,-grep("(2012)",names(kfw_grid_reshape1))]
kfw_grid_reshape3<-kfw_grid_reshape2[,-grep("(2013)",names(kfw_grid_reshape2))]
kfw_grid_reshape4<-kfw_grid_reshape3[,-grep("(2014)",names(kfw_grid_reshape3))]
kfw_grid_reshape5<-kfw_grid_reshape4[,-grep("(1981)",names(kfw_grid_reshape4))]

kfw_grid_reshape5<-kfw_grid_reshape5[,order(names(kfw_grid_reshape5))]

MeanT<-grep("MeanT_",names(kfw_grid_reshape5))
MeanP<-grep("MeanP_",names(kfw_grid_reshape5))
MinT<-grep("MinT_",names(kfw_grid_reshape5))
MaxT<-grep("MaxT_",names(kfw_grid_reshape5))
MinP<-grep("MinP_",names(kfw_grid_reshape5))
MaxP<-grep("MaxP_",names(kfw_grid_reshape5))
MaxL<-grep("MaxL_",names(kfw_grid_reshape5))
Pop<-grep("Pop_",names(kfw_grid_reshape5))

all_reshape <- c(MeanT,MeanP,MaxT,MaxP,MinP,MinT,MaxL,Pop)
psm_Long <- reshape(kfw_grid_reshape5@data, varying=all_reshape, direction="long",idvar="GridID",sep="_",timevar="Year")

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

#Add enforcement start year

#correct communities with enforcement prior to demarcation
psmtest5 <- psm_Long
psmtest5$enfdiff= psmtest5$enforce_st - psmtest5$demend_y
summary(psmtest5$enfdiff)
psmenf <- subset(psmtest5, psmtest5$enfdiff<0)
table(psmenf$reu_id)

psmtest5$enforce_st[which(psmtest5$reu_id==84)]<-2007

#create enforcement treatment var
psmtest5$trtenf <- 0
psmtest5$trtenf[which(psmtest5$Year>=psmtest5$enforce_st)]<-1
table(psmtest5$trtenf)

psm_Long <- psmtest5


write.csv(psm_Long,file="psm_Long_ALL151.csv")

##Scratch

sub29<-kfw_grid_reshape5[kfw_grid_reshape5@data$GridID==225929,]
View(sub29@data)
View(sub29@data[,(100:200)])
View(sub29@data[,(200:296)])

#----
#Old code to build panel data using SCI
##****WRONG!!****** 

# varList=c("MaxL_")
# psm_Long <- BuildTimeSeries(dta=kfw_grid20,idField="GridID",varList_pre=varList,1982,2010,colYears=c("demend_y","apprend_y","regend_y"),
#                             interpYears=c("Slope","Road_dist","Riv_Dist","UF","Elevation","terrai_are","Pop_","MeanT_","MeanP_","MaxT_",
#                                           "MaxP_","MinP_","MinT_","ntl_",
#                                           "urbtravtim","reu_id", "Id" ))
# 
# all_reshape <- c(PCloss, mean_ln, minairTemp, maxairTemp, meanairTemp, minPre, maxPre, meanPre, MinDist, DecayDist, ProjCount, DecayDist100, DecayDist25)
# DFa4 <- reshape(DFa3, varying=all_reshape,direction="long", idvar="ID", sep="_", timevar="Year")
# 
# 
# psm_Long$Year <- as.numeric(psm_Long$Year)
