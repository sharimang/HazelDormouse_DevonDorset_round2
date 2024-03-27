#### SURVEY POINT SAMPLING ####
## Generating survey locations within the reference region ##
## Sample points generated for each variable of interest while prioritizing sites with the most variable ranges represented ##
## Minimum distance between points is 3km ##



## Written by: Shari Mang ##
## Written on: February 2024 ##


# sample with distance condition: https://stackoverflow.com/questions/68757629/sf-generate-random-points-with-maximal-distance-condition/68760620#68760620
# spatstat package: https://stackoverflow.com/questions/64847597/how-do-i-generate-data-where-points-are-repelled-if-they-land-within-a-certain/64851194#64851194
# https://rdrr.io/cran/terra/man/sample.html
# sample and prune by distance: https://www.jla-data.net/eng/creating-and-pruning-random-points-and-polygons/
# spatialEco package example with subsample.distance 
  # https://gis.stackexchange.com/questions/460972/stratified-random-sampling-in-r
  # another https://gis.stackexchange.com/questions/397561/sampling-points-in-r-with-both-distance-constraint-and-group
  # stratified and controlling for distance; would need a third element to prioritize cells with max layers



#### SET UP ####
# remotes::install_github("jeffreyevans/spatialEco")

library(terra)
library(stars)
library(sp)
library(sf)
library(spatialEco)
library(dplyr)
library(tidyverse)
library(tibble)
library(nngeo)
library(janitor)
library(units)
library(here)
library(conflicted)
conflicted::conflict_prefer("here", "here")

wd.survey <- here("survey/")
wd.focus <- here("survey/focus_regions/")
wd.out <- here("model/thaw_pseudo_abs/")
wd.woodland <- here("data/UKCEH_Land_Cover/")
wd.thaw <- here("data/thaw/v3/")

bng_epsg <- "EPSG:27700"
mod.name <- "model_pseudo_allpres_01_"

# Occurrence data
occ_all <- sf::st_read(paste0(wd.survey, "occurrence_all_combined.gpkg"))
# Reference region from which to sample survey sites
ref_region_survey <- terra::rast(paste0(wd.survey, "reference_region_survey.tif"))

 

#### LOAD AND MAKE DATA TO BE SAMPLED ####

#### All regions of interest ####
# Mosaic individual variable survey areas together, summing the number of variables present in each cell
var_files <- list.files(wd.focus, full.names = TRUE) ## 22 layers
var_rasters <- sapply(var_files, function (x) terra::rast(x))
# Recrop to reference region ---> removing buffer happened after all files were written out 
var_rasters <- sapply(var_rasters, function (x) terra::mask(x, ref_region_survey))

# Mosaic rasters together
names(var_rasters) <- NULL # remove names -- seems to be needed for mosaic to work.
var_rasters$fun <- sum # sum the values everywhere rasters overlap
vars_survey <- do.call(mosaic, var_rasters)
names(vars_survey) <- "variables_per_cell"

# get variable names
var_rast_names <- sapply(var_files, function (x) strsplit((strsplit(x, paste0(wd.focus, "/"))[[1]][2]), ".tif")[[1]][1])
names(var_rasters) <- var_rast_names
var_rasters <- var_rasters[1:(length(var_rasters)-1)] # remove function from the end

saveRDS(var_rasters, file = paste0(wd.survey, "vars_survey_raster_list.RData")) # save list as .RData not .rds as it contains multiple objects
terra::writeRaster(vars_survey, paste0(wd.survey, "vars_survey_combined.tif"), overwrite = TRUE)
# vars_survey <- terra::rast(paste0(wd.survey, "vars_survey_combined.tif"))
terra::writeRaster(vars_survey, paste0(wd.survey, "vars_survey_combined_plus.tif"), overwrite = TRUE) # with thaw2_broadl raster included mosaicing 

# Calculate the area of each variable or interaction set
areas <- lapply(var_rasters, function (x) terra::expanse(x, unit = "km"))
areas_df <- do.call(rbind, areas)
areas_df <- areas_df %>%
  tibble::rownames_to_column(., "variable") %>%
  dplyr::arrange(., area) # put in order of increasing area
write.csv(areas_df, paste0(wd.survey, "variable_areas_df.csv"), row.names = FALSE)





#### Make Data ####
# Make a dataframe with the cells included in possible sample region for each variable as a column; append column with total number of variables represented in each cell 
# Turn list of rasters into a single raster stack
vars_survey_stack <- terra::rast(var_rasters)
# terra::writeRaster(vars_survey_stack, paste0(wd.survey, "vars_survey_stack.tif"), overwrite = TRUE)
vars_survey_stack <- terra::rast(paste0(wd.survey, "vars_survey_stack.tif")) 

# Extract raster values into a dataframe
vars_survey_region_df <- as.data.frame(terra::extract(vars_survey_stack, 1:ncell(vars_survey_stack)))
vars_survey_region_df$cell_num <- 1:ncell(vars_survey_stack)
vars_survey_region_df <- vars_survey_region_df %>%
  relocate(cell_num, .before = anc_wood_survey)

# Get coordinates of each pixel and number of layers in each pixel -- run as a stars object
vars_survey <- read_stars(paste0(wd.survey, "vars_survey_combined.tif")) %>%
  sf::st_set_crs(bng_epsg)

# Create points at each raster cell location 
vars_survey_points <- sf::st_as_sf(vars_survey, as_points = TRUE, na.rm = FALSE) # ordered by cell number
# Get lat/lon columns and cell_number
vars_survey_points <- vars_survey_points  %>%
  dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2]) %>%
  dplyr::mutate(cell_num = row_number()) %>%
  dplyr::rename(num_layers = vars_survey_combined.tif)

## Make a version with all NA cells removed 
vars_survey_points_NArm <- vars_survey_points %>%
  drop_na()

# Append to variable information
vars_survey_region_df <- merge(vars_survey_region_df, vars_survey_points, by = "cell_num")
# write.csv(vars_survey_region_df, paste0(wd.survey, "variables_survey_region_df.csv"), row.names = F)
vars_survey_region_df <- read.csv(paste0(wd.survey, "variables_survey_region_df.csv"))


#### SAMPLING ####

# > Set up data for loop ####
# Sampling with terra::spatSample -- requires raster format
set.seed(1000)

# Transform data to long format
# Keep lat and lon and they're needed later; don't need geometry in the loop
var_survey_l <- vars_survey_region_df %>%
  dplyr::select(-geometry) %>%
  gather(key = "variable", value = "include" , -c(cell_num, lon, lat, num_layers)) %>%
  drop_na()

# Calculate weighted selection variable - based on number of layers (variables) represented in each cell
var_survey_l$weighting <- exp(var_survey_l$num_layers)
saveRDS(var_survey_l, file = paste0(wd.survey, "variables_survey_region_long.rds")) # too big to reasonably save a csv
var_survey_l <- readRDS(paste0(wd.survey, "variables_survey_region_long.rds"))

# Order of variables for point selection -> sample variables with least area first
vars <- areas_df$variable

# Create blank raster template --> will be populated with weighting values
env_template <- terra::rast(paste0(wd.survey, "template_100m.tif"))
values(env_template) <- NA
names(env_template) <- "weighting" 

# Blank list where cell number of sampled points will be recorded 
sampled_cells <- list()

# Output for sample points --> saved into a list
sample_data_list <- list()

# Data frame of points for checking distance
pnt_cols <- c("cell", "weighting", "geometry")
prev_points <- as.data.frame(matrix(nrow=0, ncol = length(pnt_cols)))
colnames(prev_points) <- pnt_cols

# Set minimum distance 
min_dist <- 3000 
# Set number of samples
n_sample <- 250

# Data to use in loop
var_options <- var_survey_l


# > Run Sampling Loop ####
i <- 1
for(i in 1:length(vars)) {
  
  ## ~~~ Sample points for variable ~~~ ##
  # Isolate to variable of interest
  var_dat <- var_options[var_options$variable == vars[i], ]
  
  # Convert to raster
  sample_raster <- env_template # reset blank raster for each variable
  xy <- cbind(as.numeric(var_dat$lon), as.numeric(var_dat$lat)) # coordinates of the current data
  sample_raster[cellFromXY(sample_raster, xy)] <- var_dat$weighting # populate raster with weighting
  
  # Sample points weighted by number of layers
  smpl <- terra::spatSample(sample_raster, size= n_sample, method="weights", replace=FALSE, cells = TRUE, xy = TRUE)
  
  
  ## ~~~ Distance between points ~~~ ##
  # Only for currently sampled points 
  # convert to sf object
  smpl_sf <- st_as_sf(smpl, coords = c("x", "y"), crs = bng_epsg)
  # distance matrix
  smpl_dist <- st_distance(smpl_sf) %>%
    units::drop_units()
  
  # remove all values above and including the 0 diagonal line --> convert to NA
  r <- 1
  for(r in 1:nrow(smpl_dist)) {
    smpl_dist[r, (r):ncol(smpl_dist)] <- NA
  }
  
  # make a list of all locations with a distance < minimum distance
  r_rm <- list() # reset list
  c <- 1
  for(c in 1:ncol(smpl_dist)){
    if(any(smpl_dist[,c]<min_dist, na.rm = T)) {
      r_rm[[length(r_rm)+1]] <- c
    } 
  } # the output value is the row number from data used to make the distance matrix --> smpl_sf
  
  # remove proximate sample points from data --> smpl_sf
  row_rm <- unlist(r_rm, use.names = FALSE)
  if(!is.null(row_rm)){ # if there are rows to remove (i.e. it is not NULL)
    smpl_sf_update <- smpl_sf[-row_rm,] # remove those rows from data 
    # columns: cell, weighting, geometry; sf object
  }
  # if NULL make no changes -- no rows to remove.
  
  
  ## ~~~ Distance to Previous Variables' Points ~~~ ##
  # Exclude first variable -- no previous points to compare it with
  if(i != 1) {
    # Combine current and previous sample points
    check <- as.data.frame(smpl_sf_update) # convert to df to bind to template --> can't be in sf format for round 1 to work
    check_dist <- rbind(check, prev_points) 
    
    # convert back to sf to do distance matrix
    check_dist <- st_as_sf(check_dist, crs = bng_epsg)
    # calculate distance to previously selected points 
    prev_dist <- st_distance(check_dist) %>%
      units::drop_units()
    
    # remove all values above and including the 0 diagonal line --> convert to NA
    r <- 1
    for(r in 1:nrow(prev_dist)) {
      prev_dist[r, (r):ncol(prev_dist)] <- NA
    }
    
    # make a list of all locations with a distance < minimum distance
    r_rm <- list() # reset list
    c <- 1
    for(c in 1:ncol(prev_dist)){
      if(any(prev_dist[,c]<min_dist, na.rm = T)) {
        r_rm[[length(r_rm)+1]] <- c
      } 
    } # the output value is the row number from data used to make the distance matrix --> check_dist
    
    # remove proximate sample points from data --> smpl_sf_update
    row_rm <- unlist(r_rm, use.names = FALSE)
    if(!is.null(row_rm)){ # if there are rows to remove (i.e. it is not NULL)
      smpl_sf_update <- smpl_sf_update[-row_rm,] # remove those rows from data
      # columns: cell, weighting, geometry;  sf object
    } 
    # if NULL make no changes -- no rows to remove.
  } 
  
  
  ## ~~~ Select Final Points for Variable i ~~~ ###
  # Order data by highest to lowest number of layers per pixel
  smpl_sf_update <- smpl_sf_update[order(smpl_sf_update$weighting, decreasing = T),]
  
  # Keep the first n rows -- points with most layers will be at the top
  smpl_sf_update <- smpl_sf_update[1:15,] %>% # keep 15 points per variable
    drop_na() 
  ## for some, there may not be 15 points, so it will populate the dataframe with NA rows -- need to be removed so data is appended correctly
  
  # Add new points to point data that is used when checking distance
  prev_points <- rbind(prev_points, smpl_sf_update)
  
  
  ## ~~~ Output Data ~~~ ##
  # Add variable name to dataset
  smpl_sf_update$variable <- vars[i] # columns: cell, weighting, geometry, variable
  
  ### Sample Point Dataframe ###
  sample_data_list[[length(sample_data_list)+1]] <- smpl_sf_update # append to list of data sets for each variable
  
  ### Sample Point Cell Numbers ###
  sampled_cells <- append(sampled_cells, list(smpl_sf_update$cell)) # extract cell numbers and append to list of cell numbers
  
  
  ## ~~~ Remove Already Sampled Points From Data ~~~ ##
  # Remove sampled cells from dataframe before next round of sampling
  cell_rm <- unlist(sampled_cells, use.names = F)
  var_options <- var_options[!(var_options$cell_num %in% cell_rm),] # if cell_number is in this list, remove that row
  # dataframe updated for next loop through
  
}

# Convert the sample datasets into a single data frame
sample_points_data <- do.call(rbind,sample_data_list)
sample_points_data <- sample_points_data  %>%
  cbind(., st_coordinates(.))
#### DELETE ####
# # Append number of layers to output
# l <- data.frame(num_layers = unique(var_survey_l$num_layers), weighting = unique(var_survey_l$weighting))
# l <- l[order(l$num_layers),]
# round(l$weighting, 2)
# sample_points_data <- merge(sample_points_data, l, by = "weighting", all.x = TRUE) %>%
#   dplyr::relocate(weighting, .after = num_layers) %>%
#   dplyr::arrange(., variable) %>%
#   dplyr::mutate(point_id = row_number())

# Save out
sf::st_write(sample_points_data, paste0(wd.survey, "sample_points_survey.gpkg"), delete_dsn = TRUE)
readr::write_csv(sample_points_data, paste0(wd.survey, "sample_points_survey.csv"), append = FALSE)
sample_points_data <- st_read(paste0(wd.survey, "sample_points_survey.gpkg"))

# Save sampled cells 
saveRDS(sampled_cells, paste0(wd.survey, "sampled_cells_list.rds"))



#### Explore locations around Woodland patches ####
# Reference region as vector to crop woodland data to survey sample region
ref_region_vect <- terra::vect(paste0(wd.survey, "reference_region_survey_vect.shp"))

# > Conifer woodland ####
conif_shp <- terra::vect(paste0(wd.woodland, "lcm_2019_conif_woodland_devdor_dissolved.shp"))
# Shapefile made in QGIS from UKCEH Land Parcel 
  # data filtered to conifer cover (class 2)
  # "clip" function -> clipped to devdor buffer 
  # "dissolve" function -> adjacent polygons are dissolved into a single polygon; select option "keep disjoint features separate" 

# Crop to reference region
conif_shp <- terra::mask(conif_shp, ref_region_vect)
# convert to sf object
conif_shp <- st_as_sf(conif_shp)

# Calculate area 
conif_shp <- conif_shp %>%
  mutate(area_m2 = st_area(conif_shp),
         area_ha = units::set_units(st_area(conif_shp), "hectare")) %>%
  mutate(id = row_number())

sf::st_write(conif_shp, paste0(wd.woodland, "conif_poly_area_survey_region.gpkg"), delete_dsn = TRUE)
conif_shp <- st_read(paste0(wd.woodland, "conif_poly_area_survey_region.gpkg"))

# Histogram of patch size 
conif_df <- as.data.frame(conif_shp) %>%
  units::drop_units()

# Majority of points are within first bin
tmp <- classInt::classIntervals(conif_df$area_ha, n = 30, style = "equal", intervalClosure = "left")
# first bin is [0.5026214,53.52489) --> make histogram of that and then of remainder 
conif_df_sm <- conif_df %>%
  filter(., area_ha <53.52489)

ggplot(conif_df_sm, aes(x = area_ha)) + 
  geom_histogram(position="identity", fill = "grey", colour = "black") + 
  xlab("Conifer area (ha): max 53 ha")  +
  theme_classic() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))
ggsave(last_plot(), file = paste0(wd.survey, "conifer_size_small_hist.png"),  
       width = 30, height = 21, units = "cm")

# sizes outside 1st bin
conif_df_lrg <- conif_df %>%
  filter(., area_ha >53.52489)
ggplot(conif_df_lrg, aes(x = area_ha)) + 
  geom_histogram(position="identity", fill = "grey", colour = "black") + 
  xlab("Conifer area (ha): min 53 ha")  +
  theme_classic() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))
ggsave(last_plot(), file = paste0(wd.survey, "conifer_size_lrg_hist.png"),  
       width = 30, height = 21, units = "cm")



# > Broadleaf woodland ####
broadl_shp <- terra::vect(paste0(wd.woodland, "lcm_2019_broadleaf_woodland_devdor_dissolved.shp"))
# Crop to reference region
broadl_shp <- terra::mask(broadl_shp, ref_region_vect)
# convert to sf object
broadl_shp <- st_as_sf(broadl_shp)

# Calculate area 
broadl_shp <- broadl_shp %>%
  mutate(area_m2 = st_area(broadl_shp),
         area_ha = units::set_units(st_area(broadl_shp), "hectare")) %>%
  mutate(id = row_number())

sf::st_write(broadl_shp, paste0(wd.woodland, "broadl_poly_area_survey_region.gpkg"), delete_dsn = TRUE)
broadl_shp <- st_read(paste0(wd.woodland, "broadl_poly_area_survey_region.gpkg"))


# > Combine conifer and broadleaf ####
colnames(conif_shp); colnames(broadl_shp)
# Combine spatial information from both forest types into single object
woodland_shp <- rbind(conif_shp, broadl_shp)
woodland_shp <- sf::st_cast(woodland_shp, "MULTIPOLYGON") # currently is "Geometry" not polygon
sf::st_write(woodland_shp, paste0(wd.woodland, "conif_broadl_poly_area_survey_region.gpkg"), delete_dsn = TRUE)
woodland_shp <- sf::st_read(paste0(wd.woodland, "conif_broadl_poly_area_survey_region.gpkg"))

# Combine forest type polygons
# Get overall area of the combined woodland
woodland_diss <- spatialEco::sf_dissolve(woodland_shp) # dissolve polygons, regardless of forest type
woodland_diss <- woodland_diss %>%
  mutate(id = row_number(),
         area_ha = units::set_units(st_area(woodland_diss), "hectare")) %>%
  rename(geom = x) %>%
  relocate(geom, .after = last_col()) %>%
  sf::st_cast(., "MULTIPOLYGON")
sf::st_write(woodland_diss, paste0(wd.woodland, "conif_broadl_dissolved_area_survey_region.gpkg"), delete_dsn = TRUE)
woodland_diss <- sf::st_read(paste0(wd.woodland, "conif_broadl_dissolved_area_survey_region.gpkg"))

# Want to know the area of conifer within the larger woodland polygons
# intersect conifer with woodland and get the area of the intersected polygon.
x <- st_intersects(conif_shp, woodland_diss, sparse = TRUE)
# when sparse = TRUE, output is list with list element i an integer vector with all indices j for which x[i] intesects y[j]) is TRUE 
# get the area values from woodland for all intersecting locations; add that as a column to conif_shp
df <- as.data.frame(x)
df[df$row.id == 830,] ## some conifer polygons overlap 2+ woodland polygons so the row.id is duplicated
# also 877, 878
w <- df$col.id # Get the row identifier for the woodland data
df$wd_area <- woodland_diss$area_ha[w] # area values at each row listed in w - in order of w


# Get area of woodland, controlling for duplicated conifer polygon id
a_list <- list()
con_r <- unique(df$row.id)

i <- con_r[1]
for(i in con_r) {
  if(any(duplicated(df$row.id[df$row.id == i]))) # where the row.id is i, are there duplicates of the row.id 
     { # if true
    wd_area <- sum(df$wd_area[df$row.id == i])
  } else 
    { # if false
      wd_area <- df$wd_area[df$row.id == i]
    }
  tmp <- cbind(i, wd_area) # area is in hectares
  a_list[[length(a_list)+1]] <- tmp
}

wood_areas <- as.data.frame(do.call(rbind, a_list)) %>% # i is the corresponding conifer shp row
  rename(id = i,
         wd_area_ha = wd_area)
nrow(wood_areas); nrow(conif_shp) # Matches! Yay!

# Append woodland area to conifer spatial data; i = conif_shp$id
conif_shp <- merge(conif_shp, wood_areas, by = "id")
# Save out new version 
sf::st_write(conif_shp, paste0(wd.woodland, "conif_poly_area_survey_region.gpkg"), delete_dsn = TRUE)

# Histogram sizes of combined woodland 
ggplot(conif_shp, aes(x = wd_area_ha)) + 
  geom_histogram(position="identity", fill = "grey", colour = "black") + 
  xlab("Woodland Area (ha)")  +
  theme_classic() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))
ggsave(last_plot(), file = paste0(wd.survey, "woodland_area_hist.png"),  
       width = 30, height = 21, units = "cm")


#### THaW Category 1,2,5 ####
# New THaW data grouping that represents scrub and hedgerows

# Called thaw_2 in code even though it's categories 1, 2, and 5
thaw_2 <- terra::project(terra::rast(paste0(wd.thaw, "thaw_cat125_100m.tif")), bng_epsg)

# Extract values and make histogram
thaw_2_dat <- as.data.frame(terra::extract(thaw_2, 1:ncell(thaw_2)))
thaw_2_dat$cell_num <- 1:ncell(thaw_2)
thaw_2_dat <- thaw_2_dat %>% rename(thaw2_cover = thaw_cat125_100m)
thaw_2_dat$thaw2_cover[thaw_2_dat$thaw2_cover == 0] <- NA

ggplot(thaw_2_dat, aes(x = thaw2_cover)) + 
  geom_histogram(position="identity", fill = "grey", colour = "black") + 
  xlab("THaW Cat. 2 Cover")  +
  theme_classic() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))
ggsave(last_plot(), file = paste0(wd.survey, "thaw125_hist.png"),  
       width = 30, height = 21, units = "cm")

median(thaw_2_dat$thaw2_cover, na.rm = TRUE) # 0.01341644
mean(thaw_2_dat$thaw2_cover, na.rm = TRUE) # 0.01911895

thaw_2_thresh <- thaw_2_dat %>%
  dplyr::filter(thaw2_cover > 0.025) %>%
  dplyr::filter(thaw2_cover < 0.25) %>%
  drop_na()
head(thaw_2_thresh)
min(thaw_2_thresh$thaw2_cover)

ggplot(thaw_2_thresh, aes(x = thaw2_cover)) + 
  geom_histogram(position="identity", fill = "grey", colour = "black") + 
  xlab("THaW Cat. 2 [0.025, 0.25]")  +
  theme_classic() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))
ggsave(last_plot(), file = paste0(wd.survey, "thaw125_hist_025_25.png"),  
       width = 30, height = 21, units = "cm")



# > Occurrence data within THaW 2 cells ####
thaw_2 <- terra::project(terra::rast(paste0(wd.thaw, "thaw_cat125_100m.tif")), bng_epsg)

occ_sub <- occ_all[occ_all$source != "pseudo-absence", ]
occ_sub <- occ_sub %>%
  rename(X = easting, 
         Y = northing) %>%
  dplyr::select(-c(id)) %>%
  dplyr::relocate(source, .before = X)

occ_sub_pres <- occ_sub[occ_sub$occ==1,]
occ_sub_abs <- occ_sub[occ_sub$occ == 0,]

thaw_2_pres <- as.data.frame(terra::extract(thaw_2, occ_sub_pres))
thaw_2_abs <- as.data.frame(terra::extract(thaw_2, occ_sub_abs))
thaw_2_occ <- as.data.frame(terra::extract(thaw_2, occ_sub))
thaw_2_occ <- cbind(thaw_2_occ, occ_sub[, c("source", "occ")])
thaw_2_occ$occ <- as.factor(thaw_2_occ$occ)

ggplot(thaw_2_pres, aes(x = thaw_cat125_100m)) + 
  geom_histogram(position="identity", fill = "grey", colour = "black") + 
  xlab("Thaw 2: Presence locations")  +
  theme_classic() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))
ggsave(last_plot(), file = paste0(wd.survey, "thaw125_DM_presence_hist.png"),  
       width = 30, height = 21, units = "cm")

ggplot(thaw_2_abs, aes(x = thaw_cat125_100m)) + 
  geom_histogram(position="identity", fill = "grey", colour = "black") + 
  xlab("Thaw 2: Absence locations")  +
  theme_classic() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))
ggsave(last_plot(), file = paste0(wd.survey, "thaw125_DM_absence_hist.png"),  
       width = 30, height = 21, units = "cm")


ggplot(thaw_2_occ, aes(x = thaw_cat125_100m, color = occ)) + 
  geom_histogram(position="identity", fill = "white") + 
  xlab("Thaw 2: Occurrence")  +
  guides(colour=guide_legend(title="Occurrence")) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18), 
        legend.text = element_text(size = 14),
        legend.title = element_text(size=14),
        plot.title = element_text(size = 18), 
        strip.text = element_text(size = 14)) +
  facet_wrap(~ occ) 
ggsave(last_plot(), file = paste0(wd.survey, "thaw125_DM_occurrence_hist.png"),  
       width = 30, height = 21, units = "cm")


# > Sample points in THaW 2 near woodland ####
# Want to add more survey sites within THaW 2 habitat that is within a given distance of woodland patches
# >> Isolate regions of interest ####
# Crop to study region 
thaw_2 <- terra::mask(thaw_2, ref_region_survey)
terra::writeRaster(thaw_2, paste0(wd.survey, "thaw_cat125_ref_region.tif"), overwrite = TRUE)

# Threshold for thaw 2 region of interest --> minimum 5% cover
thaw_2_clamp <- terra::clamp(thaw_2, lower = 0.05, values = FALSE) 
terra::writeRaster(thaw_2_clamp, paste0(wd.survey, "thaw_cat125_above_05.tif"))

# convert to stars object
thaw_2_clamp <- st_as_stars(thaw_2_clamp)
# or read_stars()

# Convert raster to point data 
thaw_2_pnt <- sf::st_as_sf(thaw_2_clamp, as_points = TRUE, na.rm = FALSE) # Keep NA for now to assign cell numbers
thaw_2_pnt <- thaw_2_pnt %>%
  dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2]) %>%
  dplyr::mutate(cell_num = 1:ncell(thaw_2)) %>%
  dplyr::rename(thaw2_cover = thaw_cat125_100m) %>%
  tidyr::drop_na()


# Woodland data --> filter to large regions (min. 20ha)
# Broadleaf forest alone 
broadl_lrg <- st_read(paste0(wd.woodland, "broadl_poly_area.gpkg")) %>%
  dplyr::filter(., area_ha >= 20)

# Conifer alone
conif_lrg <- st_read(paste0(wd.woodland, "conif_poly_area.gpkg")) %>%
  dplyr::filter(., wd_area_ha >= 20)
min(conif_lrg$wd_area_ha)

# Woodland combined 
woodland_lrg <- st_read(paste0(wd.woodland, "conif_broadl_dissolved_area.gpkg")) %>%
  dplyr::filter(., area_ha >= 20)


# >> Calculate distance between THaW 2 cells and woodland ####
# Distance to broadleaf woodland and ONLY the broadleaf patch size is considered -- not contiguous conifer sections
dist_calc <- nngeo::st_nn(thaw_2_pnt, broadl_lrg, returnDist = TRUE)
# Don't set max distance as it then excludes points so can't match them back to point dataset
n <- unlist(dist_calc$nn)
length(n); nrow(thaw_2_pnt) # should be same length

# Put distance values for each point into a dataframe
dist_thaw2_broadl <- data.frame(cell_num = 1:length(thaw_2_pnt$cell_num),
                               sf::st_coordinates(thaw_2_pnt),
                               broadl_id = unlist(dist_calc[[1]]), # first item is the index id from the neighbouring vector
                               distance = round(unlist(dist_calc[[2]]), 2))#%>%
  #sf::st_as_sf(., coords = c("X","Y"), crs = bng_epsg)

# Put distance values into raster
env_template <- terra::rast(paste0(wd.survey, "thaw_cat125_above_05.tif"))
values(env_template) <- NA

dist_thaw2_broadl_r <- env_template 
xy <- cbind(as.numeric(dist_thaw2_broadl$X), as.numeric(dist_thaw2_broadl$Y)) # coordinates of the cells to be populated with values
dist_thaw2_broadl_r[cellFromXY(dist_thaw2_broadl_r, xy)] <- dist_thaw2_broadl$distance
names(dist_thaw2_broadl_r) <- "dist_thaw2_broadl"

# Save out 
terra::writeRaster(dist_thaw2_broadl_r, paste0(wd.survey, "distance_thaw125_broadl_woodland.tif"), overwrite = T)


# >> Isolate THaW 2 locations within distance range of interest ####
# Consider location between 400m to 1000m from woodland
dist_thaw2_broadl_smpl_rng <- terra::clamp(dist_thaw2_broadl_r, lower = 400, upper = 1000, values = FALSE)
terra::writeRaster(dist_thaw2_broadl_smpl_rng, paste0(wd.survey, "distance_thaw125_broadl_sample_range.tif"), overwrite = T)

# Filter point data to between 400m to 1000m
dist_thaw2_broadl_pnt <- dist_thaw2_broadl %>%
  sf::st_as_sf(., coords = c("X","Y"), crs = bng_epsg) %>%
  dplyr::filter(dplyr::between(distance, 400, 1000))
max(dist_thaw2_broadl_pnt$distance); min(dist_thaw2_broadl_pnt$distance)


# >> Exclude points that overlap with Woodland polygons ####
# Sample points OUTSIDE of woodland (conifer and broadleaf) -- accounts for overlap with woodland patches

# Identify which points overlap with combined woodland
pnt_inter <- sf::st_intersects(dist_thaw2_broadl_pnt, woodland_diss, sparse = T) %>% # outputs row.name of x if they intersect
  as.numeric() # populates with NA instead of empty for those that don't intersect
# Filter out points within polygons
finalpnt <- dist_thaw2_broadl_pnt[is.na(pnt_inter),] # keep only points (rows) that don't intersect --> assigned NA

# Turn possible sample locations back to raster 
finalpnt <- finalpnt %>%
  dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2])

env_template <- terra::rast(paste0(wd.survey, "thaw_cat125_above_05.tif"))
values(env_template) <- NA
sample_raster <- env_template 
names(sample_raster) <- "thaw125_broadl_sample_range"
xy <- cbind(as.numeric(finalpnt$lon), as.numeric(finalpnt$lat)) # coordinates of the current data
sample_raster[cellFromXY(sample_raster, xy)] <- finalpnt$distance
terra::writeRaster(sample_raster, paste0(wd.survey, "distance_thaw125_broadl_sample_range_rd.tif"), overwrite = T)

thaw2_broadl_region <- terra::rast(paste0(wd.survey, "distance_thaw125_broadl_sample_range_rd.tif"))

# Convert to 1 and 0 
m <- matrix(c(0,Inf,1), ncol=3, byrow=T)
thaw2_broadl_survey <- terra::classify(thaw2_broadl_region, m)
names(thaw2_broadl_survey) <- "thaw2_broadl_survey"
terra::writeRaster(thaw2_broadl_survey, paste0(wd.survey, "thaw2_broadl_survey.tif"), overwrite = TRUE)


# >> Sample THaW 2 points ####
# Sample from this list of points with a minimum distance of 3000m
thaw2_broadl_smpl <- spatialEco::subsample.distance(finalpnt, size = 100, d = 3000, replacement = F) # minimum 3km distance between sample points
# Save out full set
sf::st_write(thaw2_broadl_smpl, paste0(wd.survey, "thaw125_broadl_sample_sites_all.gpkg"))
# Remove points within min distance of existing sample sites 
sample_points_data <- st_read(paste0(wd.survey, "sample_points_survey.gpkg"))


# Combine all sample points 
# Make data formats the same 
tb_smp <- thaw2_broadl_smpl %>%
  dplyr::mutate(source = "thaw2_broadl") %>%
  dplyr::select(-c(broadl_id, distance)) %>%
  dplyr::rename(geom = geometry)

# st_as_sf(smpl, coords = c("x", "y"), crs = bng_epsg)

# Combine new sample points with previous sample points 
site <- sample_points_data %>%
  rename(source = variable, 
         cell_num = cell) %>%
  dplyr::select(-c(X, Y, num_layers, weighting, near_occ, point_id))
comb <- rbind(site, tb_smp)  

# Calculate distance between all points
comb_dist <- st_distance(comb) %>%
  units::drop_units()

# remove all values above and including the 0 diagonal line --> convert to NA
r <- 1
for(r in 1:nrow(comb_dist)) {
  comb_dist[r, (r):ncol(comb_dist)] <- NA
}
# Only want the output list to be for the thaw points which are the last 100 columns 
# Check by row instead of column, only keep rows for thaw-woodland points
# Remove sample point rows -- first 271 
comb_dist <- comb_dist[-c(1:271),]

# make a list of all locations with a distance < minimum distance
r_rm <- list() # reset list
c <- 1
for(c in 1:nrow(comb_dist)){
  if(any(comb_dist[c,]<min_dist, na.rm = T)) {
    r_rm[[length(r_rm)+1]] <- c
  } 
}
length(r_rm)

# remove proximate sample points from data thaw-woodland points
row_rm <- unlist(r_rm, use.names = FALSE)
thaw2_broadl_smpl_rm <- thaw2_broadl_smpl[-row_rm,] # remove those rows from data 
thaw2_broadl_smpl_rm <- thaw2_broadl_smpl_rm %>%
  mutate(point_id = row_number())
# Save out
sf::st_write(thaw2_broadl_smpl_rm, paste0(wd.survey, "thaw125_broadl_sample_sites.gpkg"), delete_dsn = T)
readr::write_csv(thaw2_broadl_smpl_rm, paste0(wd.survey, "thaw125_broadl_sample_sites.csv"), append = FALSE)

ggplot(thaw2_broadl_smpl_rm, 
       aes(x = distance)) + 
  geom_histogram(position="identity", fill = "grey", colour = "black") + 
  xlab("Distance from woodland")  +
  theme_classic() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))


#### Explore selected samples ####
# > Combine all sample points ####
# Load manually filtered thaw-broadleaf sample data
th_br_pt <- st_read(paste0(wd.survey, "thaw125_broadl_sample_sites_reduced.shp")) %>%
  dplyr::select(-c(fid, broadl_id, distance)) %>%
  dplyr::mutate(variable = "dist_thaw2_broadl") %>%
  dplyr::rename(geom = geometry) %>%
  cbind(., st_coordinates(.)) %>%
  dplyr::relocate(point_id, .before= geom)

# Add to other sample points 
sample_points_data <- st_read(paste0(wd.survey, "sample_points_survey.gpkg"))
survey_sites <- sample_points_data %>%
  dplyr::rename(cell_num = cell) %>%
  dplyr::select(-c(num_layers, weighting, near_occ)) %>%
  dplyr::relocate(geom, .after = last_col())

survey_sites <- rbind(survey_sites, th_br_pt)
survey_sites$point_id <- 1:nrow(survey_sites)

# >  Number of points within each variable/interactions range ####
vars_survey <- terra::rast(paste0(wd.survey, "vars_survey_combined_plus.tif"))
vars_survey_stack <- terra::rast(paste0(wd.survey, "vars_survey_stack.tif")) 
# Add thaw 2 to raster of variable ranges
thaw2_broadl_survey <- terra::rast(paste0(wd.survey, "thaw2_broadl_survey.tif"))
# combine with other survey regions 
vars_survey_stack <- c(vars_survey_stack, thaw2_broadl_survey)
names(vars_survey_stack)
# Save out updated survey stack 
terra::writeRaster(vars_survey_stack, paste0(wd.survey, "vars_survey_stack_plus.tif"), overwrite = TRUE)  


# loop through each raster and get value at each sample location 
sample_per_var <- lapply(vars_survey_stack, function (x) cbind(terra::extract(x, survey_sites), survey_sites$geom))
str(sample_per_var)
# Drop the "ID" column
sample_per_var <- lapply(sample_per_var, function(x) x[!(names(x) %in% c("ID"))])
# Combine into single dataset
sample_per_var <- sample_per_var %>%
  purrr::reduce(left_join, by = "geometry") %>%
  dplyr::mutate(point_id = row_number()) %>%
  dplyr::relocate(geometry, .after = last_col())

# For each sample site, calculate number of ranges of interest represented 
sample_per_var <- sample_per_var %>%
  dplyr::mutate(num_layers = rowSums(dplyr::select(., -c(point_id, geometry)), na.rm = T)) %>%
  dplyr::relocate(geometry, .after = last_col()) %>%
  dplyr::relocate(point_id, .before = anc_wood_survey) %>%
  st_as_sf(., crs = bng_epsg)%>%
  cbind(., st_coordinates(.))

# sf::st_write(sample_per_var, paste0(wd.survey, "sample_site_per_variable_survey.gpkg"), delete_dsn = TRUE)
# readr::write_csv(sample_per_var, paste0(wd.survey, "sample_site_per_variable_survey.csv"))
sample_per_var <- sf::st_read(paste0(wd.survey, "sample_site_per_variable_survey.gpkg"))

# Add number of layers info to survey sites data
m <- sample_per_var[, c("point_id", "num_layers")] %>%
  st_drop_geometry(.) # drop geom or merge doesn't work
survey_sites <- merge(survey_sites, m, by = "point_id")

# sf::st_write(survey_sites, paste0(wd.survey, "survey_sites.gpkg"), delete_dsn = TRUE)
# readr::write_csv(survey_sites, paste0(wd.survey, "survey_sites.csv"), append = FALSE)



# > Modelled suitability at each sampled point ####
pred <- terra::rast(paste0(wd.out, mod.name,"final_preds_thawfix.tif")) 

sample_pred <- cbind(survey_sites, terra::extract(pred$fit, survey_sites))
sample_pred$fit <- round(sample_pred$fit, 3)
summary(sample_pred$fit)
sf::st_write(sample_pred, paste0(wd.survey, "sample_points_survey_predict.gpkg"), delete_dsn = TRUE)

nrow(sample_pred[sample_pred$fit < 0.30,]) # 125 of 290
nrow(sample_pred)

# Histograms of sample point predictions
sample_df <- as.data.frame(sample_pred)

ggplot(sample_df, 
       aes(x = fit)) + 
  geom_histogram(position="identity", fill = "grey", colour = "black") + 
  xlab("Predicted Suitability")  +
  theme_classic() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))
ggsave(last_plot(), file = paste0(wd.survey, "sample_site_pred_suitability.png"),  
       width = 30, height = 21, units = "cm")


# > Distance to existing occurrence data ####
# Exclude pseudo-absences
unique(occ_all$source)
occ_sub <- occ_all[occ_all$source != "pseudo-absence", ]
occ_sub <- occ_sub %>%
  dplyr::rename(X = easting, 
         Y = northing) %>%
  dplyr::select(-c(id, occ)) %>%
  dplyr::relocate(source, .before = X)

# Combine occurrence data with sample points 
site <- survey_sites %>%
  dplyr::rename(source = variable,
         geom = geometry) %>%
  dplyr::select(-c(point_id, cell_num, num_layers)) # want same columns as in occ data
comb <- rbind(site, occ_sub)  

# Calculate distance between all points
comb_dist <- st_distance(comb) %>%
  units::drop_units()

# remove all values above and including the 0 diagonal line --> convert to NA
r <- 1
for(r in 1:nrow(comb_dist)) {
  comb_dist[r, (r):ncol(comb_dist)] <- NA
}
# only need first 290 columns -- these are sample sites; can exclude those rows as well
nrow(survey_sites) # 290
comb_dist <- comb_dist[-c(1:290), c(1:290)]

# make a list of all locations with a distance < minimum distance
r_rm <- list() # reset list
c <- 1
for(c in 1:ncol(comb_dist)){
  if(any(comb_dist[,c]<min_dist, na.rm = T)) {
    r_rm[[length(r_rm)+1]] <- c
  } 
}
length(r_rm)
# 173 of 290 sites are within 3km of previously surveyed area (60%)

# Add a column to sample point output that indicates if sample site is within 3km of existing occurrence record location
row_rm <- unlist(r_rm, use.names = FALSE)
survey_sites$near_occ <- NA
survey_sites$near_occ[row_rm] <- "Y" # is near an existing occurrence record
survey_sites$near_occ[is.na(survey_sites$near_occ)] <- "N"

# Save updates output --> near_occ = Y indicates that a existing occurrence record is within 3km of the given point
sf::st_write(survey_sites, paste0(wd.survey, "survey_sites.gpkg"), delete_dsn = TRUE)
readr::write_csv(survey_sites, paste0(wd.survey, "survey_sites.csv"), append = FALSE)



#### Refining Site Selection ####
# In QGIS, the data was cleaned with unsuitable points manually removed or shifted to more suitable locations 
# Need to recalculate number of variables represented at each site for the reduced survey site data 
# Exploring how to prioritize  points and reduce the number of sites to a reasonable number.

# Load cleaned, reduced survey site data
survey_sites <- st_read(paste0(wd.survey, "survey_sites_reduced.gpkg")) 

# > Calculate the number of sites within each variable/interactions range ####
# survey region for all variables and interactions
vars_survey_stack <- terra::rast(paste0(wd.survey, "vars_survey_stack_plus.tif"))  

# Loop through raster stack with updated points to get variables represented at each sample location 
sample_per_var_rd <- lapply(vars_survey_stack, function (x) cbind(terra::extract(x, survey_sites), survey_sites$geom))
str(sample_per_var_rd)
# Drop the "ID" column
sample_per_var_rd <- lapply(sample_per_var_rd, function(x) x[!(names(x) %in% c("ID"))])
# Combine into single dataset
sample_per_var_rd <- sample_per_var_rd %>%
  purrr::reduce(left_join, by = "geometry") %>%
  dplyr::relocate(geometry, .after = last_col())

# Append on original point id number of survey data 
j <- survey_sites %>%
  dplyr::select(c(point_id, geom)) %>%
  dplyr::rename(geometry = geom)

sample_per_var_rd <- sample_per_var_rd %>%
  dplyr::left_join(., j, by = "geometry")
# now has original point_id values for each point

# sum number of variables represented at each point
sample_per_var_rd <- sample_per_var_rd %>%
  dplyr::mutate(num_layers = rowSums(dplyr::select(., -c(point_id, geometry)), na.rm = T)) %>%
  dplyr::relocate(geometry, .after = last_col()) %>%
  dplyr::relocate(point_id, .before = anc_wood_survey) %>%
  st_as_sf(., crs = bng_epsg)%>%
  cbind(., st_coordinates(.))

## Save Reduced survey site version ##
sf::st_write(sample_per_var_rd, paste0(wd.survey, "sample_site_per_variable_survey_reduced.gpkg"), delete_dsn = TRUE)
readr::write_csv(sample_per_var_rd, paste0(wd.survey, "sample_site_per_variable_survey_reduced.csv"))

sample_per_var_rd <- sf::st_read(paste0(wd.survey, "sample_site_per_variable_survey_reduced.gpkg"))

# Summary of sites per variable
sample_per_var_simp <- sample_per_var_rd %>%
  dplyr::select(-c(X, Y, num_layers, point_id)) %>%
  sf::st_drop_geometry()

sample_summary <- colSums(sample_per_var_simp, na.rm = TRUE)
sample_summary <- as.data.frame(sample_summary) %>%
  tibble::rownames_to_column(., "variable") %>%
  dplyr::rename(represented_by_n_sites = sample_summary)
View(sample_summary)
write.csv(sample_summary, paste0(wd.survey, "sites_per_variable_reduced.csv"), row.names = F)

# Add number of layers info to survey sites data
m <- sample_per_var_rd[, c("point_id", "num_layers")] %>%
  st_drop_geometry(.) # drop geom or merge doesn't work
survey_sites <- survey_sites %>%
  dplyr::select(-num_layers) # remove old verion
survey_sites <- merge(survey_sites, m, by = "point_id")

## Save Reduced survey site version ##
sf::st_write(survey_sites, paste0(wd.survey, "survey_sites_reduced_update01.gpkg"), delete_dsn = TRUE) ## updating number of layers per site calculation
readr::write_csv(survey_sites, paste0(wd.survey, "survey_sites_reduced_update01.csv"), append = FALSE)


# > Reducing number of sites ####
# >250 sites, need to reduce to manageable amount.

# Add ranking from least to most number of sites representing variable. 
sample_summary <- sample_summary %>%
  dplyr::arrange(., represented_by_n_sites) %>%
  dplyr::mutate(priority = row_number())
readr::write_csv(sample_summary, paste0(wd.survey, "sites_per_variable_reduced.csv"), append = FALSE)

# Explore ways to remove sites
# Focus on keeping variables with fewest sites as they should also capture the more common ones 
  # ranking 1:11 captures all the interactions
# priority <- sample_summary$variable[1:11]
# priority <- c(priority, "min_win_temp_LOW_scrub_HIGH")

# Manual selection of priorities 
# Selected if there are only a few sites for that variable or if it's a variable/interaction we're really interested in 
  # selecting the less common ones will also capture the more common variables that are represented at the same sites 
  # But played with exactly what to target to get the best representation.
priority <- c("sun_spr_HIGH_conif_HIGH", "sun_spr_HIGH_scrub_HIGH", 
              "sun_spr_HIGH_max_spr_t_HIGH", "min_win_temp_HIGH_scrub_HIGH", 
              "min_win_temp_HIGH_treehedge_LOW", "rain_spr_HIGH_treehedge_LOW", 
              "rain_spr_HIGH_treehedge_HIGH", "thaw2_broadl_survey", 
              "min_win_temp_LOW_treehedge_LOW", "min_win_temp_LOW_scrub_HIGH", 
              "scrub_survey")

# Extract all point numbers where priority vars are represented. 
# Make long version of sample_per_var
sample_per_var_l <- sample_per_var_rd %>%
  dplyr::select(-c(X,Y,geometry)) %>%
  tidyr::pivot_longer(cols = anc_wood_survey:thaw2_broadl_survey, 
                      names_to = "variable", values_to = "site", values_drop_na = TRUE) %>%
  dplyr::select(-site) %>%
  dplyr::relocate(variable, .after = point_id)
  
# Extract point ID for priority variables 
priority_pnt <- sample_per_var_l$point_id[sample_per_var_l$variable %in% priority]
# remove duplicates 
priority_pnt <- unique(priority_pnt) 
length(priority_pnt)
# 145 points -> manual

# Isolate these points and look at representation across all variables. 
sample_per_var_priority <- sample_per_var_rd[sample_per_var_rd$point_id %in% priority_pnt,]

# Summary of sites per variable
sample_per_var_simp <- sample_per_var_priority%>%
  dplyr::select(-c(X, Y, num_layers, point_id)) %>%
  sf::st_drop_geometry()

sample_summary_priority <- colSums(sample_per_var_simp, na.rm = TRUE)
sample_summary_priority <- as.data.frame(sample_summary_priority) %>%
  tibble::rownames_to_column(., "variable") %>%
  dplyr::rename(represented_by_n_sites = sample_summary_priority)
#sample_summary_priority <- sample_summary_priority[order(sample_summary_priority$represented_by_n_sites),]
View(sample_summary_priority)

# Save out version 
# Manual selection 
readr::write_csv(sample_summary_priority, paste0(wd.survey, "sites_per_variable_prioritized_manual.csv"), append = F)

# Reduce spatial data to priority points 
survey_sites <- sf::st_read(paste0(wd.survey, "survey_sites_reduced_update01.gpkg")) 

# Keep points that are in priority dataset
priority_pnt <- sample_per_var_priority$point_id # points to keep
survey_sites_pr <- survey_sites[survey_sites$point_id %in% priority_pnt,]
# Save out 
# Manual selection
sf::st_write(survey_sites_pr, paste0(wd.survey, "survey_sites_reduced_priority_manual.gpkg"), delete_dsn = TRUE) 
readr::write_csv(survey_sites_pr, paste0(wd.survey, "survey_sites_reduced_priority_manual.csv"), append = FALSE)



# > Sample sites per variable summary ####
# Make the points per variable summary for the reduced data
survey_sites_pr <- sf::st_read(paste0(wd.survey, "survey_sites_reduced_priority_manual.gpkg")) 

# Extract point_ids from sample_site_per_variable_survey_reduced.gpkg --> the data priority points were selected from
sample_per_var_rd <- sf::st_read(paste0(wd.survey, "sample_site_per_variable_survey_reduced.gpkg"))
keep <- survey_sites_pr$point_id

# keep rows there point_id is a priority point
sample_per_var_pr_w <- sample_per_var_rd[sample_per_var_rd$point_id %in% keep,]
nrow(sample_per_var_pr_w)
## Save Priority sites version##
sf::st_write(sample_per_var_pr_w, paste0(wd.survey, "sample_site_per_variable_survey_priority_manual.gpkg"), delete_dsn = TRUE)
readr::write_csv(sample_per_var_pr_w, paste0(wd.survey, "sample_site_per_variable_survey_priority_manual.csv"), append = FALSE)



# > Habitat for selected sites ####
# Uploads notes with habitat assignment for survey points 
site_notes <- read.csv(paste0(wd.survey, "survey_sites_selection_notes.csv")) %>%
  janitor::clean_names()

# remove rows for manually removed points
site_notes <- site_notes[!site_notes$remove == "X",]
unique(site_notes$habitat)

# Habitat information for each point
site_hab <- site_notes %>%
  dplyr::select(c(point_id, habitat))
View(site_hab)

# Append habitat information to survey site sf for priority points
survey_sites_pr <- sf::st_read(paste0(wd.survey, "survey_sites_reduced_priority.gpkg"))
# left_join -> keeps all rows in x and adds info from y where the "by" matches
survey_sites_pr <- survey_sites_pr %>%
  dplyr::left_join(., site_hab, by = "point_id")

sf::st_write(survey_sites_pr, paste0(wd.survey, "survey_sites_reduced_priority.gpkg"), delete_dsn = TRUE) 
readr::write_csv(survey_sites_pr, paste0(wd.survey, "survey_sites_reduced_priority.csv"), append = FALSE)

# Append to sample per variable data 
survey_sites_pr_man <- sf::st_read(paste0(wd.survey, "survey_sites_reduced_priority_manual.gpkg")) 
survey_sites_pr_man <- survey_sites_pr_man %>%
  dplyr::left_join(., site_hab, by = "point_id")
sf::st_write(survey_sites_pr_man, paste0(wd.survey, "survey_sites_reduced_priority_manual.gpkg"), delete_dsn = TRUE) 
readr::write_csv(survey_sites_pr_man, paste0(wd.survey, "survey_sites_reduced_priority_manual.csv"), append = FALSE)

sample_per_var_pr_man <- sf::st_read(paste0(wd.survey, "sample_site_per_variable_survey_priority_manual.gpkg"))
sample_per_var_pr_man  <- sample_per_var_pr_man  %>%
  dplyr::left_join(., site_hab, by = "point_id")
sf::st_write(sample_per_var_pr_man, paste0(wd.survey, "sample_site_per_variable_survey_priority_manual.gpkg"), delete_dsn = TRUE)
readr::write_csv(sample_per_var_pr_man, paste0(wd.survey, "sample_site_per_variable_survey_priority_manual.csv"), append = FALSE)


# Look at distribution of habtiat types of each selection process 
ggplot(survey_sites_pr_man, aes(x = habitat)) +
  geom_bar() + 
  labs(title = "Manual prioritization") +
  theme_classic()
ggsave(last_plot(), file = paste0(wd.survey, "sites_priority_manual_habitat.png"),  
       width = 30, height = 21, units = "cm")


# > Points in Conifer ####
# Explore removing points in conifer plantations are there are loads of them. 
conif_shp <- st_read(paste0(wd.woodland, "conif_poly_area.gpkg")) # conifer woodland polygon
# Cut out all that have habitat = conifer 
conif_hab_pnt <- sf::st_intersection(survey_sites_pr, conif_shp) %>%
  dplyr::select(c(point_id, variable, habitat, area_ha, geom))
View(conif_hab_pnt)
# 57 points 
var_rm <- unique(conif_hab_pnt$variable)
int <- conif_hab_pnt$point_id

# Check how this effects individual variables 
sample_per_var_pr_w <- sf::st_read(paste0(wd.survey, "sample_site_per_variable_survey_priority.gpkg"))
# Remove rows 
s_per_var_conifRM <- sample_per_var_pr_w[!sample_per_var_pr_w$point_id %in% int,]%>%
  dplyr::select(-c(X, Y, num_layers, point_id, habitat)) %>%
  sf::st_drop_geometry()

summary_conifRM <- colSums(s_per_var_conifRM, na.rm = TRUE)
summary_conifRM <- as.data.frame(summary_conifRM) %>%
  tibble::rownames_to_column(., "variable") %>%
  dplyr::rename(represented_by_n_sites = summary_conifRM)
readr::write_csv(summary_conifRM, paste0(wd.survey, "sites_per_variable_conifRM.csv"), append = F)

# Check how this compares to overlap with broadleaf woodland 
broadl_shp <- st_read(paste0(wd.woodland, "broadl_poly_area.gpkg"))
broadl_hab_pnt <- sf::st_intersection(survey_sites_pr, broadl_shp) %>%
  dplyr::select(c(point_id, variable, habitat, area_ha, geom))
View(broadl_hab_pnt)

# cut out all from conifer survey range that are in habitat = conifer (need sample_per_var table)
conif_local <- sample_per_var_pr_w[!is.na(sample_per_var_pr_w$combo_conif_m_survey),] # 69 points
priority <- sample_summary$variable[1:11] # priority variables 

# If column is priority variable, extract point_id for rows where values = 1
pnt_list <- list()
i<- 1
for(i in 1:length(priority)) { 
  # isolate column of interest
  df <- conif_local[, c("point_id", priority[i])] %>%
    sf::st_drop_geometry()
  # fails if all rows are NA -> need if statement
  if(any(!is.na(df[[priority[i]]]))) { # are there any non-NA values
    # if true
    # remove NA rows
    df_true <- df %>%
      drop_na() %>%
      dplyr::rename(variable = priority[i])
    # put variable name in cells
    df_true[[2]] <- priority[i]
    # append to list of dataframe 
    pnt_list[[length(pnt_list)+1]] <- df_true
  } 
}

# Make into data frame 
conif_priority_overlap <- do.call(rbind,pnt_list)
View(conif_priority_overlap)
conif_priority_overlap_pnt <- unique(conif_priority_overlap$point_id)
# 69 points are within priority points --> 100% overlap. 

# count the points for each priority variable 
v <- unique(conif_priority_overlap$variable)
l <- list()
i<- 1
for(i in 1:length(v)) { 
  n <- nrow(conif_priority_overlap[conif_priority_overlap$variable == v[i],])
  l[[length(l)+1]] <- c(n, v[i])
 }

conif_pnt_per_var <- as.data.frame(do.call(rbind,l))
conif_pnt_per_var <- conif_pnt_per_var %>%
  dplyr::rename(num_points = V1,
                variable = V2)


#### Further Refinement after Manual Cleaning ####
# There were more conifer locations than we wanted; 
  # looked through any in conifer regions and decided which should be excluded from survey sites or which should just be deprioritized 
  # Also looked through points not included in selection of sites to select those that would be good backup locations. 

# Load updated survey sites that were manually adjusted 
survey_sites_pr_man <- sf::st_read(paste0(wd.survey, "survey_sites_reduced_priority_manual_edit.gpkg")) 

# Load points to consider that are not within current survey set. 
survey_sites_extra <- sf::st_read(paste0(wd.survey, "Survey_sites_extra.gpkg"))

# Add habtiat category to data 
# Uploads notes with habitat assignment for survey points 
site_notes <- read.csv(paste0(wd.survey, "survey_sites_selection_notes.csv")) %>%
  janitor::clean_names()

# remove rows for manually removed points
site_notes <- site_notes[!site_notes$remove == "X",]
unique(site_notes$habitat)

# Habitat information for each point
site_hab <- site_notes %>%
  dplyr::select(c(point_id, habitat))
View(site_hab)

# Join to extra site data
# left_join -> keeps all rows in x and adds info from y where the "by" matches
survey_sites_extra <- survey_sites_extra %>%
  dplyr::left_join(., site_hab, by = "point_id")


# Combine data sets together 
colnames(survey_sites_pr_man)
colnames(survey_sites_extra)
survey_sites_pr_man <- rbind(survey_sites_pr_man, survey_sites_extra)
# relcalcuate X and Y as some points moved
survey_sites_pr_man <- survey_sites_pr_man %>%
  dplyr::select(-c(X, Y)) %>%
  cbind(., st_coordinates(.))

sf::st_write(survey_sites_pr_man, paste0(wd.survey, "survey_sites_reduced_priority_manual.gpkg"), delete_dsn = TRUE) 
readr::write_csv(survey_sites_pr_man, paste0(wd.survey, "survey_sites_reduced_priority_manual.csv"), append = FALSE)

# > Recount the number of points for each variable of interest ####
survey_sites_pr_man <- sf::st_read(paste0(wd.survey, "survey_sites_reduced_priority_manual.gpkg"))
vars_survey_stack <- terra::rast(paste0(wd.survey, "vars_survey_stack_plus.tif")) 

# Loop through raster stack with updated points to get variables represented at each sample location 
survey_tmp <- lapply(vars_survey_stack, function (x) cbind(terra::extract(x, survey_sites_pr_man), survey_sites_pr_man$geom))
# Drop the "ID" column
survey_tmp <- lapply(survey_tmp, function(x) x[!(names(x) %in% c("ID"))])

# Combine into single dataset
survey_tmp <- survey_tmp %>%
  purrr::reduce(left_join, by = "geometry") %>%
  dplyr::relocate(geometry, .after = last_col())
View(survey_tmp)

# Append on original point id number of survey data 
j <- survey_sites_pr_man %>%
  dplyr::select(c(point_id, geom)) %>%
  dplyr::rename(geometry = geom) # to match new file

survey_tmp <- survey_tmp %>%
  dplyr::left_join(., j, by = "geometry")
# now has original point_id values for each point

# sum number of variables represented at each point
survey_tmp <- survey_tmp %>%
  dplyr::mutate(num_layers = rowSums(dplyr::select(., -c(point_id, geometry)), na.rm = T)) %>%
  dplyr::relocate(geometry, .after = last_col()) %>%
  dplyr::relocate(point_id, .before = anc_wood_survey) %>%
  st_as_sf(., crs = bng_epsg)%>%
  cbind(., st_coordinates(.))

survey_tmp <- survey_tmp %>%
  dplyr::arrange(., point_id)
sample_per_var_f <- survey_tmp

# Update the number of layers represented by each point
sf::st_write(sample_per_var_f, paste0(wd.survey, "sample_site_per_variable_survey_reduced_priority_manual_final.gpkg"), delete_dsn = TRUE)
readr::write_csv(sample_per_var_f, paste0(wd.survey, "sample_site_per_variable_survey_reduced_priority_manual_final.csv"))

# Summary of sites per variable
survey_tmp_simp <- survey_tmp %>%
  dplyr::select(-c(X, Y, num_layers, point_id)) %>%
  sf::st_drop_geometry()

sample_summary <- colSums(survey_tmp_simp, na.rm = TRUE)
sample_summary <- as.data.frame(sample_summary) %>%
  tibble::rownames_to_column(., "variable") %>%
  dplyr::rename(represented_by_n_sites = sample_summary)
View(sample_summary)
write.csv(sample_summary, paste0(wd.survey, "sites_per_variable_reduced_priority_manual_final.csv"), row.names = F)


# Add number of layers info to survey sites data
m <- survey_tmp[, c("point_id", "num_layers")] %>%
  st_drop_geometry(.) # drop geom or merge doesn't work
survey_sites_pr_man <- survey_sites_pr_man %>%
  dplyr::select(-num_layers) # remove old verion
survey_sites_pr_man <- merge(survey_sites_pr_man, m, by = "point_id")

## Save Reduced survey site version ##
sf::st_write(survey_sites_pr_man, paste0(wd.survey, "survey_sites_reduced_priority_manual_final.gpkg"), delete_dsn = TRUE) ## updating number of layers per site calculation
readr::write_csv(survey_sites_pr_man, paste0(wd.survey, "survey_sites_reduced_priority_manual_final.csv"), append = FALSE)



# > Assign priority of survey sites ####
survey_sites_pr_man <- sf::st_read(paste0(wd.survey, "survey_sites_reduced_priority_manual_final.gpkg"))
# If the target variable is any of the below, assign to Level 1 priority
priority <- c("sun_spr_HIGH_conif_HIGH", "sun_spr_HIGH_scrub_HIGH", 
              "sun_spr_HIGH_max_spr_t_HIGH", "min_win_temp_HIGH_scrub_HIGH", 
              "min_win_temp_HIGH_treehedge_LOW", "rain_spr_HIGH_treehedge_LOW", 
              "rain_spr_HIGH_treehedge_HIGH", "dist_thaw2_broadl", 
              "min_win_temp_LOW_treehedge_LOW", "min_win_temp_LOW_scrub_HIGH", 
              "scrub_survey")
survey_sites_pr_man$priority <- NA
survey_sites_pr_man$priority[survey_sites_pr_man$variable %in% priority] <- 1

length(survey_sites_pr_man$priority[!is.na(survey_sites_pr_man$priority)])
# 87 priority points
priority_point_id <- survey_sites_pr_man$point_id[!is.na(survey_sites_pr_man$priority)]

# Check how other variables are represented by the level 1 sites
p_check <- sample_per_var_f[sample_per_var_f$point_id %in% priority_point_id,]
View(p_check)
# no rain spring high- conif high
# no sun_spr_LOW_conif_HIGH

# Assign priority 1 to those variables that aren't represented 
extra <- c("rain_spr_HIGH_conif_HIGH", "sun_spr_LOW_conif_HIGH")
survey_sites_pr_man$priority[survey_sites_pr_man$variable %in% extra] <- 1
length(survey_sites_pr_man$priority[!is.na(survey_sites_pr_man$priority)]) 
# 96 points
# Now all variables are represented in the priority points.
x <- unique(survey_sites_pr_man$variable[is.na(survey_sites_pr_man$priority)])
tmp <- survey_sites_pr_man[is.na(survey_sites_pr_man$priority),]

# Assign priority level 2
priority2 <- c("anc_wood_survey", "rain_spr_LOW_conif_HIGH", 
               "sun_spr_LOW_max_spr_t_HIGH", "sun_spr_LOW_scrub_HIGH")
survey_sites_pr_man$priority[survey_sites_pr_man$variable %in% priority2] <- 2
# 25 points at level 2 
priority_point_id <- survey_sites_pr_man$point_id[!is.na(survey_sites_pr_man$priority)]
# Check how other variables are represented 
p_check <- sample_per_var_f[sample_per_var_f$point_id %in% priority_point_id,]

x <- unique(survey_sites_pr_man$variable[is.na(survey_sites_pr_man$priority)])
tmp <- survey_sites_pr_man[is.na(survey_sites_pr_man$priority),]

# Set remaining points to level 3 
survey_sites_pr_man$priority[is.na(survey_sites_pr_man$priority)] <- 3
length(survey_sites_pr_man$priority[survey_sites_pr_man$priority == 3]) # 37 least priority points

# Manually adjust priority after reviewing 
# more that should be level 1 priority
extra1 <- c(127, 228, 232, 247, 260, 266, 243, 246, 256, 260,
            244, 263, 257, 158, 269, 161, 162, 163)
survey_sites_pr_man$priority[survey_sites_pr_man$point_id %in% extra1] <- 1
# 113 sites are level 1; 25 sites level 2; 20 sites are level 3

# Remove Occ column 
survey_sites_pr_man <- survey_sites_pr_man %>%
  dplyr::select(-near_occ)

# > Save out final sites ####
sf::st_write(survey_sites_pr_man, paste0(wd.survey, "survey_sites_reduced_priority_manual_final.gpkg"), delete_dsn = TRUE) ## updating number of layers per site calculation
readr::write_csv(survey_sites_pr_man, paste0(wd.survey, "survey_sites_reduced_priority_manual_final.csv"), append = FALSE)

# Add priority value to cells of sites per variable 
sample_per_var_f <- sf::st_read(paste0(wd.survey, "sample_site_per_variable_survey_reduced_priority_manual_final.gpkg"))

# Subset point_id and priority level 
j <- survey_sites_pr_man %>% 
  dplyr::select(c(point_id, priority)) %>%
  sf::st_drop_geometry()

sample_per_var_f <- merge(sample_per_var_f, j, by = "point_id")
sf::st_write(sample_per_var_f, paste0(wd.survey, "sample_site_per_variable_survey_reduced_priority_manual_final.gpkg"), delete_dsn = TRUE)
readr::write_csv(sample_per_var_f, paste0(wd.survey, "sample_site_per_variable_survey_reduced_priority_manual_final.csv"), append = FALSE)










