#### Environmental Data Processing ####

## Most files were processed in previous iteration of project ##
## Therefore most are just loaded in as final data format ##
## Environmental data is combined with occurrence records and extracted into dataset to be used in modelling ##

## Data that is made for this project: ##
  ## THaW --> script: c01_thaw_processing.R ##

## Written by: Shari Mang ##
## Date: November 2023 ##
## Modified: December 2023 ##
## Based on code by: Regan Early, written in 2021 ##


## ~~ ## NOTE ## ~~ ##
## Climate data in env_100m7_dd.grd DOES NOT MATCH occurrence records - had zeros for most of the data points -- DO NOT USE climate variables! ##
  ## file clim100m_dd.grd has climate data at each survey location and for pseudo-absences --> use this 
  ## If modeling with survey absences, need to remake climate variables sampled to those points as survey data was MIA at time of processing. 



#### Set Up ####
# if (!require("pacman")) install.packages("pacman")
# pacman::p_install()
pacman::p_load(
  stars, # spatiotemporal data handling
  sf, 
  tidyverse, 
  here,
  terra,
  readr,
  raster,
  rgeos,
  rgdal, 
  conflicted
)
conflicted::conflict_prefer("here", "here")

bng_epsg <- "EPSG:27700"
bng <- "+init=epsg:27700"

wd.thaw <- here("data/thaw/v3/")
wd.env <- here("data/")
wd.clim <- here("data/climate/")
wd.dist <- here("data/distribution/")

env_all <- terra::rast(paste0(wd.env,"env_100m7_dd.grd")) ## the environmental data from the previous project - to be used as reference grid to ensure rasters align

devdor_buff <- terra::vect("gis/devondorset_buffer5km.shp")
crs(devdor_buff) # already bng

# devdor_buff <- st_read("gis/devondorset_buffer5km.shp")
# devdor_buff <- st_transform(devdor_buff, crs = st_crs(env_all))



#### Load environmental variables ####
# > Climate data ####
clim <- terra::project(terra::rast(paste0(wd.clim, "clim100m_dd.grd")), bng_epsg)

# > Landcover data ####
# >> Ancient woodland ####
# Note sea has value 0
anc.layers <- c("anc_rep_wood", "anc_sn_wood")
anc <- terra::rast(paste0(wd.env,"/",anc.layers,".tif")) 
anc <- terra::project(anc, bng_epsg) 
anc[["anc_wood"]] <- anc[["anc_rep_wood"]] + anc[["anc_sn_wood"]] ## the most negative forest change - used in model
#terra::writeRaster(anc, paste0(wd.env, "anc_combo.tif"), overwrite = TRUE)
anc <- terra::rast(paste0(wd.env, "anc_combo.tif"))

# >> Forest cover ####
# Broadleaf woodland
combo_broadl_c <- terra::project(rast(paste0(wd.env,"combo_broadl_c.tif")), bng_epsg)
combo_broadl_m <- terra::project(rast(paste0(wd.env,"combo_broadl_m.tif")), bng_epsg)
# Conifer
combo_conif_c <- terra::project(rast(paste0(wd.env,"combo_conif_c.tif")), bng_epsg)
combo_conif_m <- terra::project(rast(paste0(wd.env,"combo_conif_m.tif")), bng_epsg)

## Adjust values of mean coverage for forest data --> doubled due to combining them
combo_broadl_m <- combo_broadl_m/2
combo_conif_m <- combo_conif_m/2

# >> THaW Tree-hedge and scrub ####
scrub <- terra::project(terra::rast(paste0(wd.thaw, "thaw_scrub_100m.tif")), bng_epsg)
names(scrub) <- "scrub"

treehedge <- terra::project(terra::rast(paste0(wd.thaw, "thaw_treehedge_100m.tif")), bng_epsg)
names(treehedge) <- "treehedge"

# > Terrain ####
aspNS <- terra::project(terra::rast(paste0(wd.env,"OS_Terrain_100_NS.tif")), bng_epsg) ## Aspect. Extent to which a slope faces north (NS=1) or south (NS=-1), 
aspWE <- terra::project(terra::rast(paste0(wd.env,"OS_Terrain_100_WE.tif")), bng_epsg) ## Aspect. Extent to which a slope faces east (WE=1) or west (WE=-1), 
slope <- terra::project(terra::rast(paste0(wd.env,"OS_Terrain_100_slope_pct.tif")), bng_epsg) ## Slope in %



#### Combine Environmental Data ####
# Combine variables into list
env_list <- list(clim, anc, combo_broadl_c, combo_broadl_m, combo_conif_c, combo_conif_m, aspNS, aspWE, slope, scrub, treehedge)
# ensure they're in the same crs 
#env_list <- sapply(env_list, function(x) terra::project(x, bng_epsg))

# Crop to study area
env_list_devdor <- sapply(env_list, function (x) terra::crop(x, devdor_buff, mask = TRUE)) # mask = TRUE crops to exactly the lines rather than extent

# Resample so they have the same resolution --> use enviro data from previous project as the reference grid
env_list_devdor <- sapply(env_list_devdor, function(x) terra::resample(x, env_all))

# Turn into stack of rasters 
env_r <- terra::rast(env_list_devdor)
#terra::writeRaster(env_r, paste0(wd.env, "enviro_data_combined.tif"), overwrite = TRUE)
#env_r <- terra::rast(paste0(wd.env, "enviro_data_combined.tif"))



#### Occurrence data ####
# presence and pseudo-absence data compiled in previous project -> file = occ100m1990_dd; see previous project script code3_explore_100m_DD.R line 158
occ <- sf::st_read(paste0(wd.dist, "occ100m1990_dd.shp"))
crs(occ) 

# Clean data
# Absence points don't have a m100 value; there is one presence without m100 --> remove
occ <- occ[!(is.na(occ$m100) & occ$occ==1),]
occ$source[occ$occ==0] <- "pseudo-absence" # indicate that 0 are from pseudo-absences

# Save cleaned occurrence data
# sf::st_write(occ, paste0(wd.dist, "occ100m1990_dd_clean.shp"), delete_dsn = TRUE)
occ <- st_read(paste0(wd.dist, "occ100m1990_dd_clean.shp"))



#### Combine Environmental and Occurrence data ####
# create sp data frame and extract environmental data from each occ point
env_occ <- cbind(occ, terra::extract(env_r, occ))
head(env_occ)

# need to remove NAs -> Na's in occ records that don't affect model (in variables not used) so don't use drop_na() as unnecessary data will be removed
nrow(env_occ) # 1241 with NAs
# NA's are in the combo variables
env_occ <- env_occ[!is.na(env_occ$combo_broadl_c),]
colnames(env_occ)
any(is.na(env_occ[, 9:ncol(env_occ)])) # no more NA's in environmental data
# removed 75 rows
# remove random ID column 
env_occ <- env_occ %>%
  dplyr::select(-ID)

# Save out as geopackage 
sf::st_write(env_occ, here("data/env_occ_clean.gpkg"), delete_dsn = TRUE)

# Save as csv 
# write_csv() keeps the geometry column in tact 
readr::write_csv(env_occ, here("data/env_occ_clean.csv"), append = FALSE)

# This one will turn geometry column into separate x and y; if layer_options is not indicated, geometry column is removed.
# sf::st_write(env_occ, here("data/env_occ_clean.csv"), layer_options = "GEOMETRY=AS_XY", delete_dsn = TRUE)






