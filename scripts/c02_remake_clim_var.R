#### Making New Version of Climate Data ####

## Existing climate files are missing data in some survey locations; need to remake data for entire study area so they can be reampled ##
## Raw climate data is 1km resolution --> need 100m resolution ##

## Written by: Shari Mang ##
## Date:  December 2023 ##
## Based on code by: Regan Early, written in 2021 ##


#### Set up ####
install.packages("Microsoft365R", dependencies=TRUE, repos='http://cran.rstudio.com/')
library(Microsoft365R)

pacman::p_load(terra,
               sf, 
               tidyverse, 
               here, 
               magrittr, 
               conflicted)
conflicted::conflict_prefer("here", "here")

bng_epsg <- "EPSG:27700"

wd.v3 <- here("data/thaw/v3/")
wd.clim.mean <- here("data_raw/climate_mean/")
wd.clim.100 <- here("data_raw/climate_100m/")
wd.dist <- here("data/distribution/")
wd.clim <- here("data/climate/")

# 100m reference file for resampling
t100 <- terra::rast(paste0(wd.v3, "thaw_grid_template.tif"))
t100 <- terra::project(t100, bng_epsg)


# > Connect to One Drive ####
# Access files saved on One Drive
od <- Microsoft365R::get_business_onedrive(auth_type="device_code")
od$list_items("Postdoc/Dormouse")
od.clim <- "Postdoc/Dormouse/Data/5yr_means/"

# Transfer files from one drive to server cloud
od$download_folder(od.clim, dest = wd.clim.mean, overwrite = TRUE) 



#### Resample to 100m resolution ####
# There are previously made rasters with 5 year means at 1km resolution; resample the variables of interest to 100m resolution
clim.vars <- c("rainfall_spr", "sun_spr", "tasmax_spr", "tasmin_win", "tasrng_win") # variables used in model

c <- 1
i <- 1
for(c in 1:length(clim.vars)) {
  s <- terra::rast(paste0(wd.clim.mean, clim.vars[c], "_5yrmn_",1990:2019,".tif")) # load given climate variable for all years in that range
  s <- terra::resample(s, t100) # Convert to 100m resolution
  
  for (i in 1:nlyr(s)) {
    terra::writeRaster(s[[i]], paste0(wd.clim.100, names(s[[i]]),"_100m.tif")) # save out new version
  }
}



#### Extract Climate Data ####
## Create mask of occupied grid-cells for each year, multiply each year raster by the mask to nullify cells with no records that year, and add all rasters together
# Extract climate data for occurrence records 
# Sample climate data from the year of the occurrence record --> these are the climate means centred on that year.

# Climate variables of interest
clim.vars <- c("rainfall_spr", "sun_spr", "tasmax_spr", "tasmin_win", "tasrng_win")

# Load template climate raster -> already reduced to study area with correct extent, resolution, etc.
r <- terra::rast(paste0(wd.clim.100, "/tasmax_spr_5yrmn_2001_100m.tif")) ## not in BNG -- keep this way
r <- (r/r) - 1 # make all values 0

# Load occurrence records to be sampled -> all presence records and pseudo-absence 
occ <- sf::st_read(paste0(wd.dist, "occ100m1990_dd.shp")) # No survey absences in this data.
# unproject from BNG so that it matches climate raster
occ <- sf::st_transform(occ, crs(r))

# Make a column with the coordinates of each point 
occ <- occ %>%
  dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2]) ## now has lat lon in each row
head(occ)

# Extract climate data for each occurrence point
i <- 1
for(i in min(occ$year):max(occ$year)) { 
  occ.yr <- occ[occ$year==i,] # occurrence records for a given year
  xy <- cbind(as.numeric(occ.yr$lon), as.numeric(occ.yr$lat)) # coordinates of those records
  
  r.yr <- r # copy the template
  r.yr[cellFromXY(r.yr, xy)] <- 1 # for all points on raster with an occurrence record, raster value = 1
  
  # read in climate variables for given year -> in correct crs and dimensions 
  s <- terra::rast(paste0(wd.clim.100, clim.vars, "_5yrmn_", i, "_100m.tif"))
  s.yr <- s * r.yr # climate data in s will be retained in all locations with 1 in r (with occurrence records)
  names(s.yr) <- substr(names(s), 1, nchar(names(s))-5) # adjust the name structure as needed.
  assign(paste0("clim.",i), s.yr) # create an environment object clim to which values are assigned.
}

# Adds the climate at each grid-cell where species is present or absent. 
clim <- get(paste0("clim.", min(occ$year)))
for (y in min(occ$year)+1:max(occ$year)) { 
  a <- get(paste0("clim.",y))
  # a <- crop(a, extent(w))
  clim <- clim + a
  print(y)
}

names(clim) <- substr(names(s), 1, nchar(names(s))-5)

terra::writeRaster(clim, paste0(wd.clim, "clim100m_pres_pseudabs_dd.tif")) # Climate data for all presence and pseudo-absence records. 
clim <- terra::rast(paste0(wd.clim, "clim100m_pres_pseudabs_dd.tif"))


