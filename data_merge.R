
##Join extracted covar data at pixel level

library(maptools)

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
kfw_grid_temp <- read.csv("/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/Grid Data Extracts/KFW_Grids/extracted_data/terrestrial_air_temperature/extract_merge.csv")
names(kfw_grid_temp)=gsub("ad","MeanT",names(kfw_grid_temp), fixed=TRUE)
kfw_grid8=merge(kfw_grid7, kfw_grid_temp, by.x="GridID", by.y="Id")

#Precip
kfw_grid_precip <- read.csv("/Users/rbtrichler/Documents/AidData/KFW Brazil Eval/Grid Data Extracts/KFW_Grids/extracted_data/terrestrial_air_temperature/extract_merge.csv")
names(kfw_grid_precip)=gsub("ad","MeanP",names(kfw_grid_precip), fixed=TRUE)
kfw_grid9=merge(kfw_grid8, kfw_grid_precip, by.x="GridID", by.y="Id")

#yearly avg for temp and precip

#write shape poly

