##### Processing THaW data for Devon and Dorset #####

### Written by: Shari Mang ###
### Date: November 2023 ###
## Based on code by: Regan Early, written in 2021 ##


### Summary ###
## Combining thaw files into single raster then clipping it to the area of interest (Devon and Dorset)
## Reclassify the data into two categories used in model -> trees/hedges and scrub; defined by Regan Early in previous Dormouse project
## Summarize these rasters at 100m resolution
## non-thaw locations should be zero 




#### Set Up ####
library(rgdal) 
library(sf) # for spatial data, shapefiles
library(rgeos)
library(here)
library(raster)
library(stars)
library(sf)
library(conflicted)
conflicted::conflict_prefer("here", "here")

bng_epsg <- "EPSG:27700"
bng <- "+init=epsg:27700"

wd.v3 <- here("data/thaw/v3/")

OSGB.proj <- '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs'

env <- raster::raster("data/enviro_data_combined.tif")



#### Part 1: Mosaic rasters together ####
# Run on THaW cropped to Devon and Dorset #
# These were cropped in QGIS as doing it in R took much longer and the server usually just crashed.

# Load thaw files
thaw_path <- here("data_raw/thaw/dev_dor/")
thaw_files <- list.files(thaw_path, full.names = TRUE) # make a list of all the files in the folder
thaw_list <- sapply(thaw_files, function (x) raster::raster(x)) # turn all those files into a raster
# These have values from 0 to 15
# Rasters are already cropped to Devon and Dorset and are in the same projection as the environmental data --> don't need to adjust anything

# Mosaic rasters together
names(thaw_list) <- NULL # remove names -- seems to be needed for mosaic to work.
thaw_list$fun <- max # select the maximum value from each raster anywhere they overlap
thaw_raw <- do.call(mosaic, thaw_list)
raster::writeRaster(thaw_raw, paste0(wd.v3, "thaw_mosaic_devdor.tif"), overwrite = TRUE)



#### Part 2: 100m reference ####
# Make 100m resolution reference raster
t <- thaw_raw
t100 <- raster::aggregate(t, fact=100, fun=sum)
raster::writeRaster(t100, paste0(wd.v3, "template_100m_res.tif"), overwrite = TRUE)
# resolution is 100m but cells contain nonsense vales
t100 <- raster::raster(paste0(wd.v3, "template_100m_res.tif"), overwrite = TRUE)

# Reclassify cell values to be 1 in the 100m resolution template
m <- matrix(c(0,Inf,1), ncol=3, byrow=T)
t100 <- raster::reclassify(t100, m) ## Make all values 1 (excluding 0)
raster::writeRaster(t100, paste0(wd.v3, "template_100m_res_reclass.tif"), overwrite = TRUE)
t100 <- raster::raster(paste0(wd.v3, "template_100m_res_reclass.tif"))

# Resample to align with base grid (the environemental data) -> env
t100_rs <- resample(t100, env)
raster::writeRaster(t100_rs, paste0(wd.v3, "thaw_grid_template.tif"), overwrite = TRUE)
  # 100m resolution template with values of 1 for all cells that contain data



#### Part 3: Reclassify THaW data ####
# Group the data into treehedge and scrub categories

# Load 100m grid template
t100 <- raster::raster(paste0(wd.v3, "thaw_grid_template.tif"))
# Load mosaiced data for dev-dor to be reclassified
thaw_raw <- raster::raster(paste0(wd.v3, "thaw_mosaic_devdor.tif")) 

# Crop THaW data (1m) to the grid template (100m res) --> to ensure dimensions and extent are correct
t <- raster::crop(thaw_raw, t100)
raster::writeRaster(t, paste0(wd.v3, "thaw_devdor_template.tif")) ## still has raw thaw values, cropped to correct region
# t <- raster::raster(paste0(wd.v3, "thaw_devdor_template.tif"))
# Max value in raster is 255 

# > Tree-hedge ####
# Classify categories 2-4 and 6 as trees and hedges
m <- matrix(c(0,1.5,0,
              1.5,4.5,1,
              4.5,5.5,0,
              5.5,6.5,1, 
              7,Inf,0), # also need to get rid of 15 -> make it zero of it gets included into the aggregation
            nrow=5, byrow=T)
treehedge <- raster::reclassify(t, m)

raster::writeRaster(treehedge, paste0(wd.v3, "thaw_treehedge_1m.tif"), overwrite = TRUE)
treehedge <- raster(paste0(wd.v3, "thaw_treehedge_1m.tif")) # has correct extent

# Calculate the sum of 1m tree-hedge cells in each 100m grid-cell.
treehedge_100 <- raster::aggregate(treehedge, fact=100, fun=sum)/10000 # convert to proportion.
raster::writeRaster(treehedge_100, paste0(wd.v3, "thaw_treehedge_100m.tif"))
# treehedge_100 <- raster::raster(paste0(wd.v3, "thaw_treehedge_100m.tif"))  # has correct extent and dimensions


# > Scrub data #####
# Classify categories 1 and 5 as scrub below 1.3m
m <- matrix(c(0.5,1.5,1,
              1.5,4.5,0,
              4.5,5.5,1,
              5.5,6.5,0,
              7,Inf,0), # also need to get rid of 15 -> make it zero of it gets included into the aggregation
            nrow=5, byrow=T)
scrub <- raster::reclassify(t, m)

raster::writeRaster(scrub, paste0(wd.v3, "thaw_scrub_1m.tif"), overwrite = TRUE)
# scrub <- raster::raster(paste0(wd.v3, "thaw_scrub_1m.tif")) # has correct extent

# Calculate the sum of 1m scrub cells in each 100m grid-cell.
scrub_100 <- raster::aggregate(scrub, fact=100, fun=sum)/10000 # convert to proportion.
raster::writeRaster(scrub_100, paste0(wd.v3, "thaw_scrub_100m.tif"))
# scrub_100 <- raster::raster(paste0(wd.v3, "thaw_scrub_100m.tif"))


#### Part 4: Categories 2 and 3 individually ####
# Load raw thaw values
t <- raster::raster(paste0(wd.v3, "thaw_devdor_template.tif"))

# > Category 2 ####
# Isolate category 2 and make a 100m grid aggregate with proportion of pixels that are category 2
m <- matrix(c(0,1.5,0,
              1.5,2.5,1,
              2.5,Inf,0),
            nrow=3, byrow=T)
thaw_2 <- raster::reclassify(t, m)
raster::writeRaster(thaw_2, paste0(wd.v3, "thaw_cat2_1m.tif"), overwrite = TRUE)

# Calculate the sum of 1m category 2 cells in each 100m grid-cell --> convert to proportion
thaw_2_100 <- raster::aggregate(thaw_2, fact=100, fun=sum)/10000 
# Resample to 100m grid template -- need it to align with all other enviro data 
thaw_2_100 <- raster::resample(thaw_2_100, t100)
raster::writeRaster(thaw_2_100, paste0(wd.v3, "thaw_cat2_100m.tif"), overwrite = TRUE)


# > Category 3 separately ####
m <- matrix(c(0,2.5,0,
              2.5,3.5,1,
              3.5,Inf,0),
            nrow=3, byrow=T)
thaw_3 <- raster::reclassify(t, m)
raster::writeRaster(thaw_3, paste0(wd.v3, "thaw_cat3_1m.tif"), overwrite = TRUE)

# Calculate the sum of 1m category 3 cells in each 100m grid-cell --> convert to proportion
thaw_3_100 <- raster::aggregate(thaw_3, fact=100, fun=sum)/10000 
# Resample to 100m grid template -- need it to align with all other enviro data 
thaw_3_100 <- raster::resample(thaw_3_100, t100)
raster::writeRaster(thaw_3_100, paste0(wd.v3, "thaw_cat3_100m.tif"), overwrite = TRUE)


# > Category 1,2,5 Combined  ####
# Isolate category 2 and make a 100m grid aggregate with proportion of pixels that are category 2
m <- matrix(c(0,0.5,0,
              0.5,2.5,1,
              2.5,4.5,0,
              4.5,5.5,1,
              5.5,Inf,0),
            nrow=5, byrow=T)
thaw_125 <- raster::reclassify(t, m)

raster::writeRaster(thaw_125, paste0(wd.v3, "thaw_cat125_1m.tif"), overwrite = TRUE)

# Calculate the sum of 1m category 2 cells in each 100m grid-cell --> convert to proportion
thaw_125_100 <- raster::aggregate(thaw_125, fact=100, fun=sum)/10000 
# Resample to 100m grid template -- need it to align with all other enviro data 
thaw_125_100 <- raster::resample(thaw_125_100, t100)
raster::writeRaster(thaw_125_100, paste0(wd.v3, "thaw_cat125_100m.tif"), overwrite = TRUE)





#### THaW data codes ####

# Grid codes in the included THaW GeoTIFF raster datasets correspond to the below classifications: 
#   0 Null/Non THaW habitats 
# 1 Below 1.3m (Scrub, Bushes or Misc) 
# 2 Below 3m (Managed Hedgerow, Large Bushes) 
# 3 3m to 15m (Tree Canopy, Mature Hedgerow) 
# 4 Above 15m (Contiguous Tree Canopy) 
# 5 Below 1.3m (Fragmented - Scrub, Bushes or Misc e.g. fence/dry stone wall <10m2) 
# 6 Above 15m (emergent or isolated trees) 
# 15 Open water/Non THaW habitats


