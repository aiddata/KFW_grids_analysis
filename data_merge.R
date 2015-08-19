
##Join extracted covar data at pixel level

library(maptools)
library(splitstackshape)

#Obtain grid data with community level info
kfw_grid = readShapePoly("/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/Grid Data Extracts/KFW_Grids/shps/OhFive_wComm/CommunityFactors_Grid.shp")

#Drop Unused Variables at community level (NDVI, temp and precip because we're using at cell level); keep pop, treatment info, and community identifiers
kfw_grid@data <- kfw_grid@data[,-(8:744),drop=FALSE]

#Merge grid-level covariate files
#elevation
kfw_grid_elevation <- read.csv("/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/Grid Data Extracts/KFW_Grids/extracted_data/elevation/extract.csv")
names(kfw_grid_elevation)[2] = "Elevation"
kfw_grid1= merge(kfw_grid, kfw_grid_elevation, by.x="GridID", by.y="Id")

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
names(kfw_grid_NDVI)=gsub("ad", "MaxL", names(kfw_grid_NDVI), fixed=TRUE)
kfw_grid4=merge (kfw_grid3, kfw_grid_NDVI, by.x="GridID", by.y="Id")

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
air_temp <- read.csv("/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/Grid Data Extracts/KFW_Grids/extracted_data/terrestrial_air_temperature/extract_merge.csv")

for (i in 2:length(air_temp))
{
  splt <- strsplit(colnames(air_temp)[i],"_")
  splt[[1]][1] <- sub("X","",splt[[1]][1])
  month = splt[[1]][2]
  year = splt[[1]][1]
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

names(air_temp_mean)=gsub("ad_","",names(air_temp_mean), fixed=TRUE)
names(air_temp_max)=gsub("ad_","",names(air_temp_max), fixed=TRUE)
names(air_temp_min)=gsub("ad_","",names(air_temp_min), fixed=TRUE)

kfw_grid8=merge(kfw_grid7, air_temp_mean, by.x="GridID", by.y="Id")
kfw_grid9=merge(kfw_grid8, air_temp_max, by.x="GridID", by.y="Id")
kfw_grid10=merge(kfw_grid9, air_temp_min, by.x="GridID", by.y="Id")

#Precip
precip <- read.csv("/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/Grid Data Extracts/KFW_Grids/extracted_data/terrestrial_air_temperature/extract_merge.csv")

for (i in 2:length(precip))
{
  splt <- strsplit(colnames(precip)[i],"_")
  splt[[1]][1] <- sub("X","",splt[[1]][1])
  month = splt[[1]][2]
  year = splt[[1]][1]
  dt = paste(year,"-",month,sep="")
  colnames(precip)[i] <- dt
}

precip_ts <- melt(precip,id="id")
precip_ts <- cSplit(precip_ts, "variable", "-")
precip_ts_mean <- aggregate(value ~ variable_1 + id, precip_ts, FUN=mean)
precip_ts_max <- aggregate(value ~ variable_1 + id, precip_ts, FUN=max)
precip_ts_min <- aggregate(value ~ variable_1 + id, precip_ts, FUN=min)
precip_mean <- reshape(precip_ts_mean, idvar=c("id"), direction="wide", timevar="variable_1")
precip_max <- reshape(precip_ts_max, idvar=c("id"), direction="wide", timevar="variable_1")
precip_min <- reshape(precip_ts_min, idvar=c("id"), direction="wide", timevar="variable_1")

#Rename vars
for (i in 2:length(air_temp_mean))
{
  colnames(precip_mean)[i] <- sub("value.","MeanP_",colnames(precip_mean)[i])
  colnames(precip_max)[i] <- sub("value.","MaxP_",colnames(precip_max)[i])
  colnames(precip_min)[i] <- sub("value.","MinP_",colnames(precip_min)[i])
}

names(precip_mean)=gsub("ad_","",names(precip_mean), fixed=TRUE)
names(precip_max)=gsub("ad_","",names(precip_max), fixed=TRUE)
names(precip_min)=gsub("ad_","",names(precip_min), fixed=TRUE)

#names(kfw_grid_precip)=gsub("ad","MeanP",names(kfw_grid_precip), fixed=TRUE)

kfw_grid11=merge(kfw_grid10, precip_mean, by.x="GridID", by.y="Id")
kfw_grid12=merge(kfw_grid11, precip_max, by.x="GridID", by.y="Id")
kfw_grid13=merge(kfw_grid12, precip_min, by.x="GridID", by.y="Id")

## Write Final Shapefile

writePolyShape(kfw_grid13,""/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/GridDataProcessed/OhFive_gridanalysis_inputs.shp"")

