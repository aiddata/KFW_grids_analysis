
#-------------------------
# Join extracted covariate data at grid cell level
# All 151 PPTAL communities
# Add covariates and create pre-trends needed to determine high pressure communities 

library(maptools)
library(reshape)
library(splitstackshape)
library(ggplot2)
library(plyr)

library(devtools)
devtools::install_github("itpir/SCI@master")
library(SCI)
# library(stargazer)
# library(lmtest)
# library(multiwayvcov)
loadLibs()

#Obtain grid data with community level info
kfw_grid = readShapePoly("/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/Grid Data Extracts/KFW_Grids/shps/OhFive_wComm/CommunityFactors_Grid.shp")

#Drop Unused Variables at community level 
#(Pop, slope, elevation, urban travel time, distance to river/road, NDVI, temp and precip because we're using at cell level)
#keep treatment info and community identifiers
kfw_grid@data <- kfw_grid@data[,-(4:744),drop=FALSE]

#----
#correct information to include 2 PPTAL indigenous lands that were wrongly excluded from the grid level extract
#reu_id=90 and reu_id=123
# if want to produce original shapefile that miscorrectly attributed the cells and communities, run code without this section
#load in shapefile at community level from KFW_Amazon git repo
kfw_comm=readShapePoly("/Users/rbtrichler/Documents/AidData/Git Repos/KFW_Amazon/processed_data/kfw_analysis_inputs.shp")
#drop out unused variables at community level
kfw_comm<-kfw_comm[,-(3:277),drop=FALSE]
kfw_comm<-kfw_comm[,-(36:99),drop=FALSE]
# exclude everything but reu_id=90 & reu_id=123
comm<-c(90, 123)
kfw_comm1<-kfw_comm[kfw_comm@data$reu_id %in% comm,]
#create new dataframe with appropriate GridIDs and merge with comm level info
GridID<-c(393244,392458,443485,442699)
reu_id<-c(90,90,123,123)
df=data.frame(GridID, reu_id)
kfw_comm2<-join(kfw_comm1@data,df,by='reu_id',type='left',match='all')
#subset grids shapefile to 4 grid cells of interest, drop incorrect vars, merge correct vars
kfw_grid_add<-kfw_grid[(kfw_grid$GridID %in% GridID),]
kfw_grid_add<-kfw_grid_add[,-(3:36),drop=FALSE]
kfw_grid_add$reu_id<-NA
kfw_grid_add$reu_id[kfw_grid_add$GridID==393244]<-90
kfw_grid_add$reu_id[kfw_grid_add$GridID==392458]<-90
kfw_grid_add$reu_id[kfw_grid_add$GridID==443485]<-123
kfw_grid_add$reu_id[kfw_grid_add$GridID==442699]<-123
kfw_grid_add1<-merge(kfw_grid_add,kfw_comm2,by="GridID")
#drop and rename vars from merge
kfw_grid_add2<-kfw_grid_add1[,-grep("(Id)",names(kfw_grid_add1))]
kfw_grid_add3<-kfw_grid_add2[,-grep("(reu_id.x)",names(kfw_grid_add2))]
kfw_grid_add3<-rename(kfw_grid_add3,c(id="Id"))
kfw_grid_add3<-rename(kfw_grid_add3,c(reu_id.y="reu_id"))
#add information to kfw_grid, dropping out old GridID rows and adding in correct rows from kfw_comm
kfw_grid_x<-kfw_grid[!(kfw_grid$GridID %in% GridID),]
kfw_grid_add3<-kfw_grid_add3[,order(names(kfw_grid_add3))]
kfw_grid_x<-kfw_grid_x[,order(names(kfw_grid_x))]
kfw_grid_add3 <- spChFIDs(kfw_grid_add3, paste("b", row.names(kfw_grid_add3), sep="."))

kfw_grid_new=rbind(kfw_grid_x,kfw_grid_add3)

kfw_grid<-kfw_grid_new
#-----

#Drop non-PPTAL communities (but leave in those that were never demarcated)
kfw_grid$NA_check <- 0
kfw_grid$NA_check[is.na(kfw_grid$reu_id)] <- 1
kfw_grid0 <- kfw_grid[kfw_grid$NA_check != 1,]

#Merge grid-level covariate files
#elevation
#since dropped out non-PPTAL communities, there will be many records that can't be matched in merge (81363 records)
kfw_grid_elevation <- read.csv("/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/Grid Data Extracts/KFW_Grids/extracted_data/elevation/extract.csv")
names(kfw_grid_elevation)[2] = "Elevation"
kfw_grid1= merge(kfw_grid0, kfw_grid_elevation, by.x="GridID", by.y="Id")

#nighttime lights
kfw_grid_lights <- read.csv("/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/Grid Data Extracts/KFW_Grids/extracted_data/v4avg_lights_x_pct/extract_merge.csv")
names(kfw_grid_lights)=gsub("ad","ntl",names(kfw_grid_lights), fixed=TRUE)
kfw_grid2= merge(kfw_grid1, kfw_grid_lights, by.x="GridID", by.y="Id")

#access (urban travel time, to urban centers with pop>50,000)
kfw_grid_access <- read.csv("/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/Grid Data Extracts/KFW_Grids/extracted_data/access/extract.csv")
names(kfw_grid_access)[2]="urbtravtime"
kfw_grid3= merge (kfw_grid2, kfw_grid_access, by.x="GridID", by.y="Id")

#NDVI
kfw_grid_NDVI <- read.csv("/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/Grid Data Extracts/KFW_Grids/extracted_data/ltdr_yearly_max_mask_lt6k/extract_merge.csv")

kfw_grid_NDVI_dec<- kfw_grid_NDVI[,2:35]/10000
kfw_grid_NDVI_all <- cbind(kfw_grid_NDVI[,c(-2:-35)], kfw_grid_NDVI_dec)
names(kfw_grid_NDVI_all)=gsub("ad", "MaxL", names(kfw_grid_NDVI_all), fixed=TRUE)
names(kfw_grid_NDVI_all)[1]="Id"
kfw_grid4=merge (kfw_grid3, kfw_grid_NDVI_all, by.x="GridID", by.y="Id")

#Distance to River
kfw_grid_riv <- read.csv("/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/Grid Data Extracts/KFW_Grids/extracted_data/rivers_dist/extract.csv")
names(kfw_grid_riv)[2]="Riv_Dist"
kfw_grid5=merge (kfw_grid4, kfw_grid_riv, by.x="GridID", by.y="Id")

#Distance to Road
kfw_grid_road <- read.csv("/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/Grid Data Extracts/KFW_Grids/extracted_data/roads_dist/extract.csv")
names(kfw_grid_road)[2]="Road_dist"
kfw_grid6=merge(kfw_grid5, kfw_grid_road, by.x="GridID", by.y="Id")

#Slope
kfw_grid_slope <- read.csv("/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/Grid Data Extracts/KFW_Grids/extracted_data/slope/extract.csv")
names(kfw_grid_slope)[2]="Slope"
kfw_grid7=merge(kfw_grid6, kfw_grid_slope, by.x="GridID", by.y="Id")

#Temp
air_temp <- read.csv("/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/Grid Data Extracts/updated_climate_data/temp_extract_merge.csv")

for (i in 2:length(air_temp))
{
  # splt <- strsplit(colnames(air_temp)[i],"_")
  # month = splt[[1]][3]
  # year = splt[[1]][2]
  year = substr(colnames(air_temp)[i], 6, 9)
  month = substr(colnames(air_temp)[i], 10, 11)
  dt = paste(year,"-",month,sep="")
  colnames(air_temp)[i] <- dt
}

air_temp_ts <- melt(air_temp,id="Id")
air_temp_ts <- cSplit(air_temp_ts, "variable", "-")
air_temp_ts_mean <- aggregate(value ~ variable_1 + Id, air_temp_ts, FUN=mean)
air_temp_ts_max <- aggregate(value ~ variable_1 + Id, air_temp_ts, FUN=max)
air_temp_ts_min <- aggregate(value ~ variable_1 + Id, air_temp_ts, FUN=min)

air_temp_mean <- reshape(air_temp_ts_mean, idvar=c("Id"), direction="wide", timevar="variable_1")
air_temp_max <- reshape(air_temp_ts_max, idvar=c("Id"), direction="wide", timevar="variable_1")
air_temp_min <- reshape(air_temp_ts_min, idvar=c("Id"), direction="wide", timevar="variable_1")

for (i in 2:length(air_temp_mean))
{
  colnames(air_temp_mean)[i] <- sub("value.","MeanT_",colnames(air_temp_mean)[i])
  colnames(air_temp_max)[i] <- sub("value.","MaxT_",colnames(air_temp_max)[i])
  colnames(air_temp_min)[i] <- sub("value.","MinT_",colnames(air_temp_min)[i])
}

air_temp_mean <- air_temp_mean[,-(35),drop=FALSE]
air_temp_max <- air_temp_max[,-(35),drop=FALSE]
air_temp_min <- air_temp_min[,-(35),drop=FALSE]

kfw_grid8=merge(kfw_grid7, air_temp_mean, by.x="GridID", by.y="Id")
kfw_grid9=merge(kfw_grid8, air_temp_max, by.x="GridID", by.y="Id")
kfw_grid10=merge(kfw_grid9, air_temp_min, by.x="GridID", by.y="Id")

#Precip
precip <- read.csv("/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/Grid Data Extracts/updated_climate_data/precip_extract_merge.csv")

# splt splits columns when it sees _, e.g. ad_1982_02 yields 3 columns, "ad", "1982", "02" 
# month defines that the first row, represented as [[1]], and the 3rd column, [3], is the month
# year defines that the first row, [[1]], and the 2nd column, [2], is the year
# dt then renames the variable year-month, so 1982-02, rather than ad_1982_02

for (i in 2:length(precip))
{
  # splt <- strsplit(colnames(precip)[i],"_")
  # month = splt[[1]][3]
  # year = splt[[1]][2]
  year = substr(colnames(precip)[i], 6, 9)
  month = substr(colnames(precip)[i], 10, 11)
  dt = paste(year,"-",month,sep="")
  colnames(precip)[i] <- dt
}

precip_ts <- melt(precip,id="Id")
precip_ts <- cSplit(precip_ts, "variable", "-")
precip_ts_mean <- aggregate(value ~ variable_1 + Id, precip_ts, FUN=mean)
precip_ts_max <- aggregate(value ~ variable_1 + Id, precip_ts, FUN=max)
precip_ts_min <- aggregate(value ~ variable_1 + Id, precip_ts, FUN=min)
precip_mean <- reshape(precip_ts_mean, idvar=c("Id"), direction="wide", timevar="variable_1")
precip_max <- reshape(precip_ts_max, idvar=c("Id"), direction="wide", timevar="variable_1")
precip_min <- reshape(precip_ts_min, idvar=c("Id"), direction="wide", timevar="variable_1")

#Rename vars
for (i in 2:length(precip_mean))
{
  colnames(precip_mean)[i] <- sub("value.","MeanP_",colnames(precip_mean)[i])
  colnames(precip_max)[i] <- sub("value.","MaxP_",colnames(precip_max)[i])
  colnames(precip_min)[i] <- sub("value.","MinP_",colnames(precip_min)[i])
}

precip_mean <- precip_mean[,-(35),drop=FALSE]
precip_max <- precip_max[,-(35),drop=FALSE]
precip_min <- precip_min[,-(35),drop=FALSE]

kfw_grid11=merge(kfw_grid10, precip_mean, by.x="GridID", by.y="Id")
kfw_grid12=merge(kfw_grid11, precip_max, by.x="GridID", by.y="Id")
kfw_grid13=merge(kfw_grid12, precip_min, by.x="GridID", by.y="Id")

#Population (GPW v3, v4)
#drop pop columns and rename to Pop_
kfw_grid_pop <- read.csv("/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/Grid Data Extracts/KFW_Grids/extracted_data/gpw_v3.v4/merge_kfw_grid_gpwv3.v4.csv")

kfw_grid_pop<-kfw_grid_pop[,-grep("(sum)",names(kfw_grid_pop))]
kfw_grid_pop<-kfw_grid_pop[,-grep("(SP)",names(kfw_grid_pop))]
kfw_grid_pop<-kfw_grid_pop[,-grep("(Id)",names(kfw_grid_pop))]
kfw_grid_pop<-kfw_grid_pop[,-grep("(gpw_v3_density.2000)",names(kfw_grid_pop))]
kfw_grid_pop<-kfw_grid_pop[,-grep("(gpw_v4_density.2015)",names(kfw_grid_pop))]
kfw_grid_pop<-kfw_grid_pop[,-grep("(gpw_v4_density.2020)",names(kfw_grid_pop))]

colnames(kfw_grid_pop)<-sub("gpw_v3_density.","Pop_",colnames(kfw_grid_pop))
colnames(kfw_grid_pop)<-sub("gpw_v4_density.","Pop_",colnames(kfw_grid_pop))
colnames(kfw_grid_pop)<-sub(".mean","",colnames(kfw_grid_pop))

#interpolate so a value exists for every year from 1982-2010
#use each 5 year segment to fill in between, and 1990-1995 change to fill in pre-1982
#1982-1995
kfw_grid_pop$Pop_90.95<-(kfw_grid_pop$Pop_1995-kfw_grid_pop$Pop_1990)/5
kfw_grid_pop$Pop_1991<-kfw_grid_pop$Pop_1990+kfw_grid_pop$Pop_90.95
kfw_grid_pop$Pop_1992<-kfw_grid_pop$Pop_1991+kfw_grid_pop$Pop_90.95
kfw_grid_pop$Pop_1993<-kfw_grid_pop$Pop_1992+kfw_grid_pop$Pop_90.95
kfw_grid_pop$Pop_1994<-kfw_grid_pop$Pop_1993+kfw_grid_pop$Pop_90.95

kfw_grid_pop$Pop_1989<-kfw_grid_pop$Pop_1990-kfw_grid_pop$Pop_90.95
kfw_grid_pop$Pop_1988<-kfw_grid_pop$Pop_1989-kfw_grid_pop$Pop_90.95
kfw_grid_pop$Pop_1987<-kfw_grid_pop$Pop_1988-kfw_grid_pop$Pop_90.95
kfw_grid_pop$Pop_1986<-kfw_grid_pop$Pop_1987-kfw_grid_pop$Pop_90.95
kfw_grid_pop$Pop_1985<-kfw_grid_pop$Pop_1986-kfw_grid_pop$Pop_90.95
kfw_grid_pop$Pop_1984<-kfw_grid_pop$Pop_1985-kfw_grid_pop$Pop_90.95
kfw_grid_pop$Pop_1983<-kfw_grid_pop$Pop_1984-kfw_grid_pop$Pop_90.95
kfw_grid_pop$Pop_1982<-kfw_grid_pop$Pop_1983-kfw_grid_pop$Pop_90.95

#1995-2000
kfw_grid_pop$Pop_95.2000<-(kfw_grid_pop$Pop_2000-kfw_grid_pop$Pop_1995)/5
kfw_grid_pop$Pop_1996<-kfw_grid_pop$Pop_1995+kfw_grid_pop$Pop_95.2000
kfw_grid_pop$Pop_1997<-kfw_grid_pop$Pop_1996+kfw_grid_pop$Pop_95.2000
kfw_grid_pop$Pop_1998<-kfw_grid_pop$Pop_1997+kfw_grid_pop$Pop_95.2000
kfw_grid_pop$Pop_1999<-kfw_grid_pop$Pop_1998+kfw_grid_pop$Pop_95.2000

#2000-2005
kfw_grid_pop$Pop_2000.05<-(kfw_grid_pop$Pop_2005-kfw_grid_pop$Pop_2000)/5
kfw_grid_pop$Pop_2001<-kfw_grid_pop$Pop_2000+kfw_grid_pop$Pop_2000.05
kfw_grid_pop$Pop_2002<-kfw_grid_pop$Pop_2001+kfw_grid_pop$Pop_2000.05
kfw_grid_pop$Pop_2003<-kfw_grid_pop$Pop_2002+kfw_grid_pop$Pop_2000.05
kfw_grid_pop$Pop_2004<-kfw_grid_pop$Pop_2003+kfw_grid_pop$Pop_2000.05

#2005-2010
kfw_grid_pop$Pop_2005.10<-(kfw_grid_pop$Pop_2010-kfw_grid_pop$Pop_2005)/5
kfw_grid_pop$Pop_2006<-kfw_grid_pop$Pop_2005+kfw_grid_pop$Pop_2005.10
kfw_grid_pop$Pop_2007<-kfw_grid_pop$Pop_2006+kfw_grid_pop$Pop_2005.10
kfw_grid_pop$Pop_2008<-kfw_grid_pop$Pop_2007+kfw_grid_pop$Pop_2005.10
kfw_grid_pop$Pop_2009<-kfw_grid_pop$Pop_2008+kfw_grid_pop$Pop_2005.10

#Drop unused columns for merge
kfw_grid_pop<-kfw_grid_pop[,-grep("(Pop_90.95)",names(kfw_grid_pop))]
kfw_grid_pop<-kfw_grid_pop[,-grep("(Pop_95.2000)",names(kfw_grid_pop))]
kfw_grid_pop<-kfw_grid_pop[,-grep("(Pop_2000.05)",names(kfw_grid_pop))]
kfw_grid_pop<-kfw_grid_pop[,-grep("(Pop_2005.10)",names(kfw_grid_pop))]

kfw_grid13.1=merge(kfw_grid13, kfw_grid_pop, by="GridID")

kfw_grid13<-kfw_grid13.1

##Write dataset prior to merging in high pressure covars and creating pre-trends
writePolyShape(kfw_grid13,"/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/GridDataProcessed/OhFive_gridanalysis_inputs_ALL151.shp")

##---
#Merge in high pressure covars (from Ash) at the community level: Vegetation, Soil Type Zones
# will drop out those to be merged in at grid level instead (distance to conservation units, logging, mining, railways)

kfw_grid_veg <- read.csv("/Users/rbtrichler/Google Drive/REU/KfW/Inputs/HighPressureCommCovars__Abbrev_Ash.csv")
kfw_grid_veg <- kfw_grid_veg[,-(10:27),drop=FALSE]
kfw_grid14=merge(kfw_grid13, kfw_grid_veg, by.x="Id", by.y="id")

## Merge in high pressure covars (from Ash) at the community level: Crop Value and Yield

kfw_grid_crop=readShapePoly("/Users/rbtrichler/Google Drive/REU/KfW/Inputs/HighPressureCovars_CropYield/Indig_ag.shp")
kfw_grid_crop@data<-kfw_grid_crop@data[,-(1:13),drop=FALSE]
kfw_grid_crop@data<-kfw_grid_crop@data[,-(2:21),drop=FALSE]
kfw_grid_crop@data<-kfw_grid_crop@data[,-(183:184),drop=FALSE]
kfw_grid15=merge(kfw_grid14, kfw_grid_crop@data, by.x="Id", by.y="id")

## Merge in high pressure covars (from Ash) at the GRID level: Distance to conservation units, mining, logging, railways
#Distance to Federal Conservation Unit
kfw_grid_fedcons<-read.csv("/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/Grid Data Extracts/KFW_Grids/extracted_data/dist_from_conservation_federal_units/dfcf_e.csv")
names(kfw_grid_fedcons)[2]="fedcon_dist"
kfw_grid16=merge(kfw_grid15, kfw_grid_fedcons, by.x="GridID", by.y="Id")

#Distance to State Conservation Unit
kfw_grid_stcons<-read.csv("/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/Grid Data Extracts/KFW_Grids/extracted_data/dist_from_conservation_state_units/dfcs_e.csv")
names(kfw_grid_stcons)[2]="stcon_dist"
kfw_grid17=merge(kfw_grid16, kfw_grid_stcons, by.x="GridID", by.y="Id")

#Distance to Logging Center
kfw_grid_log<-read.csv("/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/Grid Data Extracts/KFW_Grids/extracted_data/dist_from_logging_center/dflc_e.csv")
names(kfw_grid_log)[2]="log_dist"
kfw_grid18=merge(kfw_grid17, kfw_grid_log, by.x="GridID", by.y="Id")

#Distance to Mining Activity
kfw_grid_mine<- read.csv("/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/Grid Data Extracts/KFW_Grids/extracted_data/dist_from_mining_interests/dfmi_e.csv")
names(kfw_grid_mine)[2]="mine_dist"
kfw_grid19=merge(kfw_grid18, kfw_grid_mine, by.x="GridID", by.y="Id")

#Distance to Railway
kfw_grid_rail<-read.csv("/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/Grid Data Extracts/KFW_Grids/extracted_data/dist_from_railways/dfrw_e.csv")
names(kfw_grid_rail)[2]="rail_dist"
kfw_grid20=merge(kfw_grid19, kfw_grid_rail, by.x="GridID", by.y="Id")

## Eliminate non-PPTAL demarcated communities

# kfw_grid20$NA_check <- 0
# kfw_grid20$NA_check[is.na(kfw_grid20$demend_y)] <- 1
# kfw_grid21 <- kfw_grid20[kfw_grid20$NA_check != 1,]
# kfw_grid20 <- kfw_grid21

## Create pre-trends
kfw_grid20$pt_NDVI_max <- timeRangeTrend(kfw_grid20,"MaxL_[0-9][0-9][0-9][0-9]",1982,1995,"GridID")

kfw_grid20$pt_temp_mean <- timeRangeTrend(kfw_grid20,"MeanT_[0-9][0-9][0-9][0-9]",1982,1995,"GridID")
kfw_grid20$pt_temp_max <- timeRangeTrend(kfw_grid20,"MaxT_[0-9][0-9][0-9][0-9]",1982,1995,"GridID")
kfw_grid20$pt_temp_min <- timeRangeTrend(kfw_grid20,"MinT_[0-9][0-9][0-9][0-9]",1982,1995,"GridID")

kfw_grid20$pt_prec_mean <- timeRangeTrend(kfw_grid20,"MeanP_[0-9][0-9][0-9][0-9]",1982,1995,"GridID")
kfw_grid20$pt_prec_max <- timeRangeTrend(kfw_grid20,"MaxP_[0-9][0-9][0-9][0-9]",1982,1995,"GridID")
kfw_grid20$pt_prec_min <- timeRangeTrend(kfw_grid20,"MinP_[0-9][0-9][0-9][0-9]",1982,1995,"GridID")

kfw_grid20$pt_ntl <- timeRangeTrend(kfw_grid20,"ntl_[0-9][0-9][0-9][0-9]",1992,1995,"GridID")

kfw_grid20$pt_cv <- timeRangeTrend(kfw_grid20,"cv[0-9][0-9][0-9][0-9]",1994,1995,"GridID")
kfw_grid20$pt_cy <- timeRangeTrend(kfw_grid20,"cy[0-9][0-9][0-9][0-9]",1991,1995,"GridID")
kfw_grid20$pt_rv <- timeRangeTrend(kfw_grid20,"rv[0-9][0-9][0-9][0-9]",1994,1995,"GridID")
kfw_grid20$pt_ry <- timeRangeTrend(kfw_grid20,"ry[0-9][0-9][0-9][0-9]",1991,1995,"GridID")
kfw_grid20$pt_sov <- timeRangeTrend(kfw_grid20,"sov[0-9][0-9][0-9][0-9]",1994,1995,"GridID")
kfw_grid20$pt_soy <- timeRangeTrend(kfw_grid20,"soy[0-9][0-9][0-9][0-9]",1991,1995,"GridID")
kfw_grid20$pt_suv <- timeRangeTrend(kfw_grid20,"suv[0-9][0-9][0-9][0-9]",1994,1995,"GridID")
kfw_grid20$pt_suy <- timeRangeTrend(kfw_grid20,"suy[0-9][0-9][0-9][0-9]",1991,1995,"GridID")
kfw_grid20$pt_wv <- timeRangeTrend(kfw_grid20,"wv[0-9][0-9][0-9][0-9]",1994,1995,"GridID")


## Write Final Shapefile, with pre-trends
writePolyShape(kfw_grid20,"/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/GridDataProcessed/OhFive_gridanalysis_inputs_pretrends_ALL151.shp")

#kfw_grid20 = readShapePoly("OhFive_gridanalysis_inputs_pretrends_ALL151.shp")

##Scratch

test<-kfw_grid[kfw_grid@data$GridID==392458,]
