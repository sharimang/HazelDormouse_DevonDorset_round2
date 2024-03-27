#### Process Climate Data to Make Most Recent 5 Year Means ####

## Needed most recent climate data for survey site selection ##

## Data source: 
  ## Monthly data from: https://catalogue.ceda.ac.uk/uuid/46f8c1377f8849eeb8570b8ac9b26d86
  ## Data: HadUK-Grid Gridded Climate Observations on a 1km grid over the UK, v1.2.0.ceda (1836-2022)
  ## It is from "HadUK-Grid gridded and regional average climate observations for the UK" which is updated annually.
  ## Data is at 1km resolution 

## At time of writing, data only available until 2022, so averaging 2018 to 2022
## Climate variables needed: Spring Rainfall, Max Spring Temperature, Min Winter Temperature, Sunshine Spring


## Written by: Shari Mang ##
## Date: February 2024 ##
## Based on code by: Regan Early, written in 2021 ##




#### SET UP ####

library(terra)

bng_epsg <- "EPSG:27700"

wd.env <- here("data_raw/climate_hadUK/")
wd.out <- here("data_raw/climate_monthly/")
wd.v3 <- here("data/thaw/v3/")
wd.clim.100 <- here("data_raw/climate_100m/")

months <- c("jan","feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")
spring <- c("mar", "apr", "may")
winter <- c("dec", "jan","feb")

devdor_buff <- terra::vect("gis/devondorset_buffer5km.shp")

# 100m reference file for resampling
t100 <- terra::rast(paste0(wd.v3, "thaw_grid_template.tif"))
t100 <- terra::project(t100, bng_epsg)



#### Rainfall variables (sum) ####
for (y in 2018:2022) {
  b <- terra::project(terra::rast(paste0(wd.env, "rainfall_hadukgrid_uk_1km_mon_",y,"01-",y,"12.nc")), bng_epsg)
  names(b) <- months 
 
  spr <- sum(b[[spring]]) # March - May
  names(spr) <- paste0("rainfall_spr_",y)
  terra::writeRaster(spr, paste0(wd.out, "rainfall_spr_",y,".tif"), overwrite = TRUE)
}

#### Sunshine hours variables (mean) ####
for (y in 2018:2022) {
  b <-  terra::project(terra::rast(paste0(wd.env, "sun_hadukgrid_uk_1km_mon_",y,"01-",y,"12.nc")), bng_epsg)
  names(b) <- months 

  spr <- mean(b[[spring]]) ## March - May
  names(spr) <- paste0("sun_spr_",y)
  terra::writeRaster(spr, paste0(wd.out, "sun_spr_",y,".tif"), overwrite = TRUE)
}

#### Minimum temperature variables (mean) ####
for (y in 2018:2022) {
  b <- terra::project(terra::rast(paste0(wd.env, "tasmin_hadukgrid_uk_1km_mon_",y,"01-",y,"12.nc")), bng_epsg)
  names(b) <- months 
 
  win <- mean(b[[winter]]) ## December - February
  names(win) <- paste0("tasmin_win_",y)
  terra::writeRaster(win, paste0(wd.out, "tasmin_win_",y,".tif"), overwrite = TRUE)
}

#### Maximum temperature variables (mean) ####
for (y in 2018:2022) {
  b <- terra::project(terra::rast(paste0(wd.env, "tasmax_hadukgrid_uk_1km_mon_",y,"01-",y,"12.nc")), bng_epsg)
  names(b) <- months 
  
  spr <- mean(b[[spring]]) ## March - May
  names(spr) <- paste0("tasmax_spr_",y)
  terra::writeRaster(spr, paste0(wd.out, "tasmax_spr_",y,".tif"), overwrite = TRUE)
}


##### Calculate 5 year means #####
# NA value is > 3e20... remove this so it doesn't get included when resampled to 100m

# Rainfall
rain <- mean(terra::rast(paste0(wd.out, "rainfall_spr_", c(2018:2022),".tif")))
# crop to Devon and Dorset
rain <- terra::crop(rain, devdor_buff, mask = TRUE)
rain <- terra::clamp(rain, upper = 10000000, values = FALSE) # changes huge values to NA
names(rain) <- "rainfall_spr_5yrmn_2022"
terra::writeRaster(rain, paste0(wd.out, "rainfall_spr_5yrmn_2022.tif"), overwrite = TRUE)

# Sunshine hours
sun <- mean(terra::rast(paste0(wd.out, "sun_spr_", c(2018:2022),".tif")))
# crop to Devon and Dorset
sun <- terra::crop(sun, devdor_buff, mask = TRUE)
sun <- terra::clamp(sun, upper = 10000000, values = FALSE) 
names(sun) <- "sun_spr_5yrmn_2022"
terra::writeRaster(sun, paste0(wd.out, "sun_spr_5yrmn_2022.tif"), overwrite = TRUE)

# Minimum temperature variables (mean)
mint <- mean(terra::rast(paste0(wd.out, "tasmin_win_", c(2018:2022),".tif")))
# crop to Devon and Dorset
mint <- terra::crop(mint, devdor_buff, mask = TRUE)
mint <- terra::clamp(mint, upper = 10000000, values = FALSE) 
names(mint) <- "tasmin_win_5yrmn_2022"
terra::writeRaster(mint, paste0(wd.out, "tasmin_win_5yrmn_2022.tif"), overwrite = TRUE)

# Maximum temperature variables (mean)
maxt <- mean(terra::rast(paste0(wd.out, "tasmax_spr_", c(2018:2022),".tif")))
# crop to Devon and Dorset
maxt <- terra::crop(maxt, devdor_buff, mask = TRUE)
maxt <- terra::clamp(maxt, upper = 10000000, values = FALSE) 
names(maxt) <- "tasmax_spr_5yrmn_2022"
terra::writeRaster(maxt, paste0(wd.out, "tasmax_spr_5yrmn_2022.tif"), overwrite = TRUE)


#### Resample to 100m resolution ####
rast.list <- list(rain, sun, mint, maxt)

i <- 1
for(i in 1:length(rast.list)) {
  rs <- terra::resample(rast.list[[i]], t100)
  terra::writeRaster(rs, paste0(wd.clim.100, names(rs), "_100m.tif"), overwrite = T)
}






