#### Predict and Evaluate Model ####

## Generating predictions and evaluating the model ##
## Corresponds to model code in c05m01_sdm_pseudo-abs.R ##
## Model 01 --> All presence data and pseudo-absences, incorporating THaW data ##

## Written by: Shari Mang ##
## Date: November 2023 ##
## Based on code by: Regan Early, written in 2021 ##


# Notes to self: 
# my dat_stnd is Regan's dat.m
# Regan's "dat_30Nov2021_dd.csv" is the combined enviro and occurrence data used in model -> equal to my "data/env_occ_clean.csv"



#### Set Up ####
pacman::p_load(dplyr, MuMIn, car, snow, lme4, performance, see, glmmTMB, DHARMa, here, effects, raster, terra, pROC, dismo) 
# dismo -> MESS with rasters
# https://easystats.github.io/performance/ 

bng_epsg <- "EPSG:27700"

wd.dat <- here("data/")
wd.env <- here("data/")
wd.out <- here("model/thaw_pseudo_abs/")
wd.clim.100 <- here("data_raw/climate_100m/")
wd.dist <- here("data/distribution/")

mod.name <- "model_pseudo_allpres_01_"

destdize <- read.csv(paste0(wd.out, "model_pseudo_allpres_01_dat_means_sd_outlier_rm.csv"))
# row 1 = mean, row 2 = SD

# load final model 
load(file = paste0(wd.out, mod.name, "final_averaged_thawfix")) # name = final_m01



#### Test Robustness ####
# Run model with same variables using different subset of calibration data 
dat <- read.csv(here("data/env_occ_clean_scrub_outliers_rm.csv"))
dat_stnd <- read.csv(paste0(wd.out, mod.name, "dat_stnd_scrub_outliers_removed.csv")) # standardized data
dat.c <- read.csv(paste0(wd.out, mod.name, "dat_stnd_cal.csv")) # data used in model
dat.v <- read.csv(paste0(wd.out, mod.name, "dat_stnd_val.csv")) # data for validation

# Variables in the final model
vars <- as.vector(colnames(final_m01$coefficients))
vars <- vars[!grepl("Intercept|I|:", vars)] ## remove intercept, quadratic, and interaction terms from the variables

# Predict model with validation data
p.val <- predict(final_m01, dat.v, re.form=NA, type="response")

roc_obj <- pROC::roc(dat.v$occ, p.val)
auc(roc_obj) 
roc_obj$auc
## 0.7409



#### Predict Suitability for Devon and Dorset ####

# > load/make environmental data for whole study area ###

# ~~~ Have been updated with correct THaW rasters ~~~ # 
# ~~~ Both enviro_data_combined_full_clim_m01.tif and enviro_data_combined_full_clim_stnd_m01.tif ~~~ #

# Environmental data rasters 
env_r <- terra::rast(paste0(wd.dat, "enviro_data_combined.tif")) # Stack of enviro rasters
# full enviro coverage EXCEPT for climate variables which only have data for occurrence records --> need full coverage for prediction
names(env_r)
vars.all <- c('combo_broadl_m', 'combo_broadl_c', 'combo_conif_m', 'combo_conif_c', 'anc_wood', 
              'treehedge', 'scrub',
              'OS_Terrain_100_NS', 'OS_Terrain_100_WE','OS_Terrain_100_slope_pct', 
              "tasrng_win_5yrmn", "tasmax_spr_5yrmn", "tasmin_win_5yrmn", 
              "sun_spr_5yrmn", "rainfall_spr_5yrmn")
# reduce to the possible rasters to be included in models
env_r <- env_r[[vars.all]] # get data for all climate variables included initially


# > Make new climate rasters ####
  # want climate data for all cells, not just ones with occurrence data
clim.vars <- c("rainfall_spr", "sun_spr", "tasmax_spr", "tasmin_win", "tasrng_win")

i <- 2019 # Use the 5 year climate mean calculated for the most recent period
s <- terra::rast(paste0(wd.clim.100, clim.vars, "_5yrmn_", i, "_100m.tif"))
crs(s) <- crs(env_r)
sp <- terra::project(s, bng_epsg) # project so same crs as env data

nn <- names(s) # list of raster names

i <- 1
for (i in 1:length(nn)) {
  n2 <- substr(nn[i], start=1, stop=nchar(nn[i])-5) ## trim year from end so you have a name that corresponds to that in the envr file (e.g., sun_spr_5yrmn)
  env_r[[n2]] <- sp[[i]] # replace raster in env with equivalent variable for 2019
}
plot(env_r$tasrng_win_5yrmn) 

# terra::writeRaster(env_r, file=paste0(wd.dat, "enviro_data_combined_full_clim_m01.tif"), overwrite = TRUE)
env_r <- terra::rast(paste0(wd.dat, "enviro_data_combined_full_clim_m01.tif"))
names(env_r)


# > Standardize environmental data #### 
# same method as for data included in the model
vars.all
n <- vars.all[1]
for (n in vars.all) {
  std <- destdize[,n] # values of model input environment data to standardise raster values
  env_r[[n]] <- (env_r[[n]] - std[1]) / std[2] # subtract original mean and divide by original standard deviation
}

# Save out standardized version of raster
terra::writeRaster(env_r, paste0(wd.dat, "enviro_data_combined_full_clim_stnd_m01.tif"), overwrite = TRUE)

env_r <- terra::rast(paste0(wd.dat, "enviro_data_combined_full_clim_stnd_m01.tif"))
env_r <- env_r[[vars]] # reduce to just variables included in this model run
names(env_r)



#### Make the prediction for the averaged model ####
# Use standardized data to use for predictions
p <- predict(env_r, final_m01, re.form=NA, type="response", full=T, se.fit=T) 
# Note: type="response" option tells R to output the prediction is on the scale of the response variable; 
# for a default binomial model the default predictions are of log-odds (probabilities on logit scale) and type = "response" gives the predicted probabilities. 
writeRaster(p, file=paste0(wd.out, mod.name,"final_preds_thawfix.tif"), overwrite = TRUE)
p <- raster::raster(paste0(wd.out, mod.name,"final_preds_thawfix.tif"))

# Standard error of the prediction
p.mod01.se <- p$se.fit
writeRaster(p.mod01.se, file=paste0(wd.out, mod.name, "final_preds_thawfix_se.tiff"), overwrite = TRUE)


#### MESS ####
occ <- sf::st_read(paste0(wd.dist, "occ100m1990_dd_clean.shp")) # load occurrence records used on model
# Mess requires raster format
env_r <- as(env_r, "Raster")
env.pa.mat <- raster::extract(env_r, occ)

mess.out <- dismo::mess(x=env_r, v=env.pa.mat, full=TRUE)
names(mess.out) <- c(vars, "rmess")
plot(mess.out[["rmess"]])

for(i in names(mess.out)) {
  writeRaster(mess.out[[i]], paste0(wd.out, mod.name, i, "_mess.tif"), overwrite = TRUE) ## the MESS surfaces for all variables together and all individual variables
}

# Mask of negative values --> indicates dissimilar environments 
mess.out <- terra::rast(paste0(wd.out, "model_pseudo_allpres_01_rmess_mess.tif"))
mess_mask <- mess.out < 0
mess_mask[mess_mask == 0] <- NA
plot(mess_mask)
writeRaster(mess_mask, paste0(wd.out, mod.name, "MESS_mask_negval.tif"), overwrite = TRUE)


# Mask of all non-negative values
mess.out <- terra::rast(paste0(wd.out, "model_pseudo_allpres_01_rmess_mess.tif"))
mess_positive <- mess.out
mess_positive[mess_positive < 0] <- NA # remove all negative values
# make all remaining values 1 -- just want a mask of non-negative values
m <- matrix(c(0, Inf, 1),
            nrow = 1, byrow = T)
mess_positive <- terra::classify(mess_positive, m, include.lowest = TRUE)
plot(mess_positive)
writeRaster(mess_positive, paste0(wd.out, mod.name, "MESS_non-negative.tif"), overwrite = TRUE)



#### Calculate thresholds ####

## Inspect the sensitivity and specificity
# p must be Raster format (not SpatRaster)
p <- raster::raster(paste0(wd.out, mod.name,"final_preds_thawfix.tif"))
p.pred <- raster::extract(p, occ)
r <- roc(occ$occ, p.pred)
r # 0.7477

df <- data.frame(spec=r$specificities, sens=r$sensitivities, thresh=r$thresholds)
colnames(df) <- c("specificity", "sensitivity", "threshold")
write.csv(df, paste0(wd.out, mod.name, "thresh_sen_spec.csv"), row.names = F)

df[which.min(abs(0.95-df$sens)), ] ## The threshold that yields 95% sensitivity: 0.2349814
df[which.min(abs(df$spec-df$sens)), ] ## The threshold that minimises the difference between sensitivity and specificity: 0.5267021
df[which.max(df$spec+df$sens), ]## Maximized difference between sensitivity and specificity: 0.501773


# >  Reclassify prediction raster according to thresholds ####
# 95% sensitivity threshold
rcl <- matrix(c(0,0.2349814,0, 0.2349814,1,1), nrow=2, byrow=T) 
p95 <- reclassify(p, rcl)
raster::writeRaster(p95, file=paste0(wd.out, mod.name, "final_preds_thawfix_SenThresh95.tiff"), overwrite = TRUE)
t <- terra::rast(paste0(wd.out, mod.name, "final_preds_thawfix_SenThresh95.tiff"))
# Threshold that maximizes difference between sensitivity and specificity
rcl <- matrix(c(0,0.501773,0, 0.501773,1,1), nrow=2, byrow=T) 
maxss <- reclassify(p, rcl)
raster::writeRaster(maxss, file=paste0(wd.out, mod.name, "final_preds_thawfix_MaxSenSpec.tiff"), overwrite = TRUE)
maxss <- rast(paste0(wd.out, mod.name, "final_preds_thawfix_MaxSenSpec.tiff"))
## Threshold that minimizes difference between sensitivity and specificity
rcl <- matrix(c(0,0.5267021,0, 0.5267021,1,1), nrow=2, byrow=T) 
minss <- reclassify(p, rcl)
raster::writeRaster(minss, file=paste0(wd.out, mod.name, "final_preds_thawfix_MinSenSpec.tiff"), overwrite = TRUE)
minss <- rast(paste0(wd.out, mod.name, "final_preds_thawfix_MinSenSpec.tiff"))
# Thresholded value at survey points
# Minimum threshold
minss <- raster::raster(paste0(wd.out, mod.name, "final_preds_thawfix_MinSenSpec.tiff"))
# All Presence points
occ_pres <- occ[occ$occ == 1, ]
occ_pres_ss <- as.data.frame(terra::extract(minss, occ_pres))
occ_pres_ss <- cbind(occ_pres_ss , occ_pres[, c("source", "occ", "geometry")])
colnames(occ_pres_ss)[1] <- "Min_SenSpec_Threshold"
sf::st_write(occ_pres_ss, paste0(wd.out, mod.name, "MinSenSpec_pres_occ.gpkg"), delete_dsn = TRUE)

sum(occ_pres_ss$Min_SenSpec_Threshold, na.rm = TRUE); nrow(occ_pres_ss) # 421 of 620 (0.67) presence points captured by threshold 
length(which(is.na(occ_pres_ss$Min_SenSpec_Threshold))) # 0 NAs

# pseudo-absence points
occ_abs <- occ[occ$occ == 0, ]
occ_abs_ss <- as.data.frame(terra::extract(minss, occ_abs))
occ_abs_ss <- cbind(occ_abs_ss , occ_abs[, c("source", "occ", "geometry")])
colnames(occ_abs_ss)[1] <- "Min_SenSpec_Threshold"
sf::st_write(occ_abs_ss, paste0(wd.out, mod.name, "MinSenSpec_abs_occ.gpkg"), delete_dsn = TRUE)

nrow(occ_abs_ss[occ_abs_ss$Min_SenSpec_Threshold == 0,]); nrow(occ_abs_ss) # 422 of 621 (0.68) absence points captured by threshold  
length(which(is.na(occ_abs_ss$Min_SenSpec_Threshold))) # 2 NAs


# NDMP Presence
occ_ndmp <- occ
occ_ndmp <- occ_ndmp[occ_ndmp$source == "NDMP", ]

occ_ndmp_ss <- as.data.frame(terra::extract(minss, occ_ndmp))
occ_ndmp_ss <- cbind(occ_ndmp_ss, occ_ndmp[, c("source", "occ", "geometry")])
colnames(occ_ndmp_ss)[1] <- "Min_SenSpec_Threshold"
sf::st_write(occ_ndmp_ss, paste0(wd.out, mod.name, "MinSenSpec_NDMP_occ.gpkg"), delete_dsn = TRUE)

sum(occ_ndmp_ss$Min_SenSpec_Threshold, na.rm = TRUE); nrow(occ_ndmp_ss) # 82 of 96 (0.85) 
length(which(is.na(occ_ndmp_ss$Min_SenSpec_Threshold))) # 0 NAs

# LRD+NDD Presence
occ_lrdndd <- occ[occ$source == "LRD+NDD", ]
occ_lrdndd_ss <- as.data.frame(terra::extract(minss, occ_lrdndd))
occ_lrdndd_ss <- cbind(occ_lrdndd_ss , occ_lrdndd[, c("source", "occ", "geometry")])
colnames(occ_lrdndd_ss)[1] <- "Min_SenSpec_Threshold"
sf::st_write(occ_lrdndd_ss, paste0(wd.out, mod.name, "MinSenSpec_LRDNDD_occ.gpkg"), delete_dsn = TRUE)

sum(occ_lrdndd_ss$Min_SenSpec_Threshold, na.rm = TRUE); nrow(occ_lrdndd_ss) # 339 of 524 (0.65) 
length(which(is.na(occ_lrdndd_ss$Min_SenSpec_Threshold))) # 0 NAs




#### Calculate R-squared for all models in best subset ####
dat.c <- read.csv(paste0(wd.out, mod.name, "dat_stnd_cal.csv"))

# best models from final dredge before averaging -> mlqid.comb_best
load(file=paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_best_vars_best_mods"))

rsq.fn <- function(x) {
  (x$null.deviance - x$deviance) / x$null.deviance
}

rsq <- unlist(lapply(mlqid.comb_best, rsq.fn))
mean(rsq) # 0.1824049
sd(rsq) # 0.001799512




