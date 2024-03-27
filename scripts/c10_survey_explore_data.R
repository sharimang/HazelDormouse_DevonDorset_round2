#### EXPLORING DATA FOR SURVEY DESIGN #####
## Exploring model outputs, error, etc. to see where the model could be improved ##
## Exploring relationship between variable values and model outputs ##
## Create ranges of interest for each variable and interaction from which we want to place survey sites ##


## Written by: Shari Mang ## 
## Date: February 2024 ##



#### 1. SET UP ####
pacman::p_load(ggplot2, RColorBrewer, classInt, colorspace, terra, spatialEco, dyplr, here, conflicted)
conflicted::conflict_prefer("here", "here")

wd.dat <- here("data/")
wd.dist <- here("data/distribution/")
wd.out <- here("model/thaw_pseudo_abs/")
wd.plot <- here("model/thaw_pseudo_abs/plots/")
wd.env <- here("data/")
wd.survey <- here("survey/")
wd.clim.100 <- here("data_raw/climate_100m/")

mod.name <- "model_pseudo_allpres_01_"

bng_epsg <- "EPSG:27700"


# > Load Files ####

# All values from environmental data
env_vals_raw <- read.csv(paste0(wd.out, mod.name, "env_raster_values.csv"))  ## Has updated THaW values ##
vars.plus <- colnames(env_vals_raw)[2:ncol(env_vals_raw)] # first column is cell number -- this also has climate means 2022
vars <- vars.plus[!grepl("_2022", vars.plus)] # keep only the variables included in the model
  
# Predicted suitability and Standard Error
pred <- terra::rast(paste0(wd.out, mod.name,"final_preds_thawfix.tif")) 

var_names <- c("Mean conifer cover", "Aspect WE",
               "Mean spring rainfall", "Scrub",
               "Mean spring sun", "Max spring temp", 
               "Treehedge", "Slope",
               "Ancient woodland", "Change broadleaf cover",
               "Aspect NS", "Min winter temp")

# All enviro data as raster, wiith full coverage of climate data (2019 centred 5-year mean)
env_r <- terra::rast(paste0(wd.dat, "enviro_data_combined_full_clim_m01_plus.tif"))

# Study area
devdor_buff <- terra::vect("gis/devondorset_buffer5km.shp")
#crs(devdor_buff) 

# Devon-Dorset county borders
dd_county <- st_read("gis/os_boundary/bdline_gb_DevDor.shp")

# NDMP sites 
ndmp_bounds <- terra::vect("gis/NDMP_boundaries_Oct21_devdor.shp")
#crs(ndmp_bounds)

# All (non-zero) values of environmental variables and values from survey points --> not matched to cell number
vals_true_sample <- read.csv(paste0(wd.out, mod.name, "compare_variable_values.csv"))
vals_true <- vals_true_sample[vals_true_sample$source == "raster", ]

# Occurrence data 
# All presence and pseudo-absence (cleaned - used in model)
occ <- st_read(paste0(wd.dist, "occ100m1990_dd_clean.shp"))

# Survey absences -- full dataset, from previous project 
surv_abs <- st_read(paste0(wd.dist, "surv_ab_24Jan2022.shp"))

# Urban areas --> to be excluded from possible survey sites 
st_layers(paste0(wd.dat, "OS_Built_Up_Areas/OS_Open_Built_Up_Areas.gpkg")) 
urban <- st_read(paste0(wd.dat, "OS_Built_Up_Areas/OS_Open_Built_Up_Areas.gpkg"), layer = "OS_Open_Built_Up_Areas")




#### 2. RELATIONSHIP BETWEEN VARIABLES AND MODEL ####
# Extract prediction and SE values into dataframe
pred_dat <- as.data.frame(terra::extract(pred, 1:ncell(pred))) # has both prediction and standard error
pred_dat$cell_num <- 1:ncell(pred)
head(pred_dat) # has both fit and se.fit

# Combine environmental data with predictions
env_pred <- merge(pred_dat, env_vals_raw, by = "cell_num")
View(env_pred)

mean(env_pred$fit, na.rm = TRUE) # 0.4043283
median(env_pred$fit, na.rm = TRUE)  # 0.3930963
mean(env_pred$se.fit, na.rm = TRUE) # 0.05304083
median(env_pred$se.fit, na.rm = TRUE)  # 0.04759482
max(env_pred$se.fit, na.rm = TRUE) # 1.079294

# Make data subsets for plotting
# SE< 0.2 and Pred. >0.25
ep_sub1 <- env_pred[(env_pred$fit > 0.25 & env_pred$se.fit < 0.2), ]
ep_sub1 <- ep_sub1 %>%
  drop_na()
max(ep_sub1$scrub)

# High error but high prediction -> SE>0.2 and Pred. >0.25
ep_sub2 <- env_pred[(env_pred$fit > 0.25 & env_pred$se.fit > 0.2), ]
ep_sub2 <- ep_sub2 %>%
  drop_na()




#### Data distribution Suitability and Standard Error ####
pred_hist <- pred_dat[,c("fit", "se.fit")]
pred_hist <- pred_dat %>%
  gather(key = "metric", value = "value") %>%
  drop_na()

# Predicted suitability
ggplot(subset(pred_hist, metric %in% "fit"), 
       aes(x = value)) + 
  geom_histogram(position="identity", fill = "grey", colour = "black") + 
  xlab("Predicted Suitability")  +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))
ggsave(last_plot(), file = paste0(wd.plot, mod.name, "hist_pred_suitability.png"),  
       width = 30, height = 21, units = "cm")

# Standard error
ggplot(subset(pred_hist, metric %in% "se.fit"), 
       aes(x = value)) + 
  geom_histogram(position="identity", fill = "grey", colour = "black") + 
  xlab("Standard Error")  +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))
ggsave(last_plot(), file = paste0(wd.plot, mod.name, "hist_standard_error.png"),  
       width = 30, height = 21, units = "cm")


#### Explore data ####
# Classify outliers in SE
# Calculate z score
pred_dat$se_zscore <- (abs(pred_dat$se.fit - mean(pred_dat$se.fit, na.rm = TRUE))/sd(pred_dat$se.fit, na.rm = TRUE))
length(which(pred_dat$se_zscore>3)) # 17990 occurrences where z-score > 3

# remove rows with a zscore greater than 3
pred_dat_subset <- pred_dat[!pred_dat$se_zscore>3, ]
max(pred_dat_subset$se.fit, na.rm = TRUE) # 0.2009323
median(pred_dat_subset$se.fit, na.rm = TRUE) # 0.06051902
# 965230 cells remaining 

# Check what intervals for SE
se_inter <- classInt::classIntervals(pred_dat_subset$se.fit, n = 3, style = "equal", intervalClosure = "left") 

# Make raster for each interval 
pred.se <- pred$se.fit
# Low
se_inter$brks[2] # 0.06698305
pred_se_low <- terra::clamp(pred.se, upper = 0.06698304, values = FALSE)
names(pred_se_low) <- "se.fit_low"
# Medium
se_inter$brks[3] # 0.1339577
pred_se_med <- terra::clamp(pred.se, lower = se_inter$brks[2], upper = 0.1339576, values = FALSE)
names(pred_se_med) <- "se.fit_med"
# High
pred_se_high <- terra::clamp(pred.se, lower = se_inter$brks[3], values = FALSE)
names(pred_se_high) <- "se.fit_high"

terra::writeRaster(pred_se_low, paste0(wd.out, mod.name, "pred_se_low.tif"), overwrite = TRUE)
terra::writeRaster(pred_se_med, paste0(wd.out, mod.name, "pred_se_med.tif"), overwrite = TRUE)
terra::writeRaster(pred_se_high, paste0(wd.out, mod.name, "pred_se_high.tif"), overwrite = TRUE)



#### Plot Variable VS SE --> Subset to SE<0.2 and Pred >0.25 ####

# Hex plot 
i <- 1
for(i in 1:length(vars)) {
  ymin <- min(ep_sub1[, vars[i]], na.rm = TRUE)
  ymax <- max(ep_sub1[, vars[i]], na.rm = TRUE)
  
  p <- ggplot(ep_sub1, aes(x = se.fit, y = .data[[vars[i]]])) + 
    geom_hex() +
    colorspace::scale_fill_continuous_sequential(palette = "Viridis",
                                                 labels = scales::label_comma(),
                                                 trans = "log",
                                                 name = "log(count)"
                                                 ) +
    ylab(var_names[i]) + 
    xlab("Standard Error") +
    ggtitle("Data: SE<0.2 and Pred.>0.25") +
    scale_y_continuous(limits = c(ymin, ymax)) + 
    scale_x_continuous(limits = c(0,0.20)) + 
    geom_smooth(color = "red") + 
    theme_classic() + 
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 18), 
          legend.text = element_text(size = 14),
          legend.title = element_text(size=14),
          plot.title = element_text((size = 20)))
  ggsave(p, file = paste0(wd.plot, mod.name, "SE_subset_vs_", vars[i], ".png"),  
         width = 30, height = 21, units = "cm")
}


# # Point version 
# for(i in 1:length(vars)) {
#   ymin <- min(ep_sub1[, vars[i]], na.rm = TRUE)
#   ymax <- max(ep_sub1[, vars[i]], na.rm = TRUE)
#   
#   p <- ggplot(ep_sub1, aes(x = se.fit, y = .data[[vars[i]]])) + 
#     geom_point(size = 0.75)+
#     ylab(vars[i]) + 
#     xlab("Standard Error") +
#     ggtitle("Data: SE<0.2 and Pred.>0.25") +
#     scale_y_continuous(limits = c(ymin,ymax)) + 
#     scale_x_continuous(limits = c(0,0.20)) + # SE cut off at 0.20
#     geom_smooth() + 
#     theme_classic() + 
#     theme(axis.text = element_text(size = 16),
#           axis.title = element_text(size = 18))
#   ggsave(p, file = paste0(wd.plot, mod.name, "SE_subset_vs_", vars[i], ".png"),  
#          width = 30, height = 21, units = "cm")
# }


#### Plot Variable VS Predicted Suitability --> Subset to SE<0.2 and Pred >0.25 ####
# Point with SE colour
i <- 1
for(i in 1:length(vars)) {
  ymin <- min(ep_sub1[, vars[i]], na.rm = TRUE)
  ymax <- max(ep_sub1[, vars[i]], na.rm = TRUE)
  
  p <- ggplot(ep_sub1, aes(x = fit, y = .data[[vars[i]]])) + 
    geom_point(size = 1.25, aes(colour = se.fit))+
    colorspace::scale_color_continuous_sequential(palette = "Purple-Yellow", rev = FALSE,
                                                  name = "SE") +
    ylab(var_names[i]) + 
    xlab("Predicted Suitability") +
    ggtitle("Data: SE<0.2 and Pred.>0.25") +
    scale_y_continuous(limits = c(ymin, ymax)) + 
    scale_x_continuous(limits = c(0.25,1)) + # Pred above 0.25
    geom_smooth(color = "red") + 
    theme_classic() + 
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 18), 
          legend.text = element_text(size = 14),
          legend.title = element_text(size=14),
          plot.title = element_text((size = 20)))
  
  ggsave(p, file = paste0(wd.plot, mod.name, "pred_sub_pt_vs_", vars[i], ".png"),  
         width = 30, height = 21, units = "cm")
}

# Hex plot
i <- 1
for(i in 1:length(vars)) {
  ymin <- min(ep_sub1[, vars[i]], na.rm = TRUE)
  ymax <- max(ep_sub1[, vars[i]], na.rm = TRUE)
  
  p <- ggplot(ep_sub1, aes(x = fit, y = .data[[vars[i]]])) + 
    geom_hex() +
    colorspace::scale_fill_continuous_sequential(palette = "Viridis",
                                                 labels = scales::label_comma(),
                                                 trans = "log",
                                                 name = "log(count)"
                                                 ) +
    ylab(var_names[i]) + 
    xlab("Predicted Suitability") +
    ggtitle("Data: SE<0.2 and Pred.>0.25") +
    scale_y_continuous(limits = c(ymin, ymax)) + 
    scale_x_continuous(limits = c(0.25,1)) + 
    geom_smooth(color = "red") + 
    theme_classic() + 
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 18), 
          legend.text = element_text(size = 14),
          legend.title = element_text(size=14),
          plot.title = element_text((size = 20)))
  
  ggsave(p, file = paste0(wd.plot, mod.name, "pred_sub_hex_vs_", vars[i], ".png"),  
         width = 30, height = 21, units = "cm")
}



#### Plot Variable VS Predicted Suitability --> High error: Subset to SE>0.2 and Pred >0.25 ####
i <- 1
for(i in 1:length(vars)) {
  ymin <- min(ep_sub2[, vars[i]], na.rm = TRUE)
  ymax <- max(ep_sub2[, vars[i]], na.rm = TRUE)
  
  p <- ggplot(ep_sub2, aes(x = fit, y = .data[[vars[i]]])) + 
    geom_point(size = 0.75, aes(colour = se.fit))+
    colorspace::scale_color_continuous_sequential(palette = "Purple-Yellow", rev = FALSE,
                                                  name = "SE") +
    ylab(var_names[i]) + 
    xlab("Predicted Suitability") +
    ggtitle("Data: SE>0.2 and Pred.>0.25") +
    scale_y_continuous(limits = c(ymin, ymax)) + 
    scale_x_continuous(limits = c(0.25,1)) + # Pred above 0.25
    geom_smooth(color = "red") + 
    theme_classic() + 
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 18), 
          legend.text = element_text(size = 14),
          legend.title = element_text(size=14),
          plot.title = element_text((size = 20)))
  
  ggsave(p, file = paste0(wd.plot, mod.name, "pred_pnt_highSE_vs_", vars[i], ".png"),  
         width = 30, height = 21, units = "cm")
}



#### Plot Variable VS SE --> All data ####
all_plot <- env_pred %>%
  drop_na()
max(all_plot$se.fit, na.rm = TRUE) # max SE = 1.079
# Hex plot 
i <- 1
for(i in 1:length(vars)) {
  ymin <- min(all_plot[, vars[i]], na.rm = TRUE)
  ymax <- max(all_plot[, vars[i]], na.rm = TRUE)
  
  p <- ggplot(all_plot, aes(x = se.fit, y = .data[[vars[i]]])) + 
    geom_hex() +
    colorspace::scale_fill_continuous_sequential(palette = "Viridis",
                                                 labels = scales::label_comma(),
                                                 trans = "log",
                                                 name = "log(count)"
    ) +
    ylab(var_names[i]) + 
    xlab("Standard Error") +
    ggtitle("Data: Full dataset") +
    scale_y_continuous(limits = c(ymin, ymax)) + 
    scale_x_continuous(limits = c(0,1.1)) + 
    geom_smooth(color = "red") + 
    theme_classic() + 
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 18), 
          legend.text = element_text(size = 14),
          legend.title = element_text(size=14),
          plot.title = element_text((size = 20)))
  ggsave(p, file = paste0(wd.plot, mod.name, "SE_all_vs_", vars[i], ".png"),  
         width = 30, height = 21, units = "cm")
}


#### Plot Variable VS Predicted Suitability --> All data ####
# Point with SE colour
i <- 1
for(i in 1:length(vars)) {
  ymin <- min(all_plot[, vars[i]], na.rm = TRUE)
  ymax <- max(all_plot[, vars[i]], na.rm = TRUE)
  
  p <- ggplot(all_plot, aes(x = fit, y = .data[[vars[i]]])) + 
    geom_point(size = 0.75, aes(colour = se.fit))+
    colorspace::scale_color_continuous_sequential(palette = "Purple-Yellow", rev = FALSE,
                                                  name = "SE") +
    ylab(var_names[i]) + 
    xlab("Predicted Suitability") +
    ggtitle("Data: Full dataset") +
    scale_y_continuous(limits = c(ymin, ymax)) + 
    scale_x_continuous(limits = c(0.0,1)) + 
    geom_smooth(color = "red") + 
    theme_classic() + 
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 18), 
          legend.text = element_text(size = 14),
          legend.title = element_text(size=14),
          plot.title = element_text((size = 20)))
  
  ggsave(p, file = paste0(wd.plot, mod.name, "pred_all_pt_vs_", vars[i], ".png"),  
         width = 30, height = 21, units = "cm")
}




#### 3. SURVEY THRESHOLDS FOR VARIABLES ####

#### Create Base Layer ####
# Reference region from which survey sites can be selected
# Using treehedge and scrub cover as initial threshold --> excluding locations without enough vegetation cover/habitat

# > Explore Treehedge ####
# Above 0.02 cover
th_gt.02 <- terra::clamp(env_r$treehedge, lower = 0.02, values = FALSE)
#writeRaster(th_gt.02, paste0(wd.survey, "treehedge_above_02.tif"))
th_gt.02 <- terra::rast(paste0(wd.survey, "treehedge_above_02.tif"))
area_th_gt.02 <- terra::expanse(th_gt.02, unit = "km") # 9646.8 km^2

## Get the difference with threshold and with prediction >0.25 
# Sensitivity threshold 
pred_thresh <- terra::rast(paste0(wd.out, mod.name, "final_preds_thawfix_MaxSenSpec.tiff"))
# area of the threshold that is not included in the treehedge threshold 
diff_thresh_th02 <- terra::mask(pred_thresh, th_gt.02, inverse = TRUE)
# remove zero and calculate the area
diff_thresh_th02 <- ifel(diff_thresh_th02 == 1, 1, NA)
plot(diff_thresh_th02)
writeRaster(diff_thresh_th02, paste0(wd.survey, "diff_threshold_treehedge_02.tif"), overwrite = TRUE)

area_pred_thresh <-terra::expanse(pred_thresh, unit = "km") # 11027.98 km^2
area_diff_thresh_th02 <- terra::expanse(diff_thresh_th02, unit = "km") # 40.5 km^2 (0.004 of total threshold area)



# > Explore Scrub ####
scrub <- env_r$scrub
plot(scrub)
scr_gt_0001 <- terra::clamp(scrub, lower = 0.00009, values = FALSE)
writeRaster(scr_gt_0001, paste0(wd.survey, "scrub_above_0001.tif"), overwrite = TRUE)
# Fills in a lot of the threshold-treehedge gaps but also adds a lot more other area
# still a lot of gaps in Dorset - where most of the gaps are
# much of it is open fields or really small hedge bits in/between fields 
# values of scrub in many gaps is 0 (scanned in Dorset) so reducing scrub threshold won't help

# Isolate high scrub areas 
# Quantile values in var_summary
# Try 20th quantile from previous exploration -> lower value = 0.006117703
scrub_q20 <- terra::clamp(scrub, lower = 0.006117703, values = F) 
a <- terra::expanse(scrub_q20, unit = "km") # 483.1 km^2
# add in some of the scrubby areas that are missed by the treehedge filter
terra::writeRaster(scrub_q20, paste0(wd.survey, "scrub_q20.tif"), overwrite = TRUE)
scrub_q20 <- terra::rast(paste0(wd.survey, "scrub_q20.tif"))

# upper quantile 
int <- classInt::classIntervals(env_vals_raw$scrub, n = 4, style = "quantile", intervalClosure = "left") # na.rm = TRUE
# upper quantile = [0.001699923,0.9995893] 

#Destandardize values 
std <- seq(-1, 4, length.out=6) 
destd <- (std * destdize[2,"scrub"]) + destdize[1,"scrub"]
# -1  0  1  2  3  4
# -0.0003109697  0.0013037511  0.0029184720  0.0045331928 0.0061479137  0.0077626345



# > Combine into reference region ####
# Combine treehedge >0.02 with scrub >0.006117703 (20th quantile of 20)
th02_scr006 <- terra::mosaic(th_gt.02, scrub_q20, fun = "sum")
# writeRaster(th02_scr006, paste0(wd.survey, "combined_treehedge02_scrub006.tif"), overwrite = T)
th02_scr006 <- terra::rast(paste0(wd.survey, "combined_treehedge02_scrub006.tif"))

# Remove existing occurrence records
# Combine occurrence data and survey absences into single shapefile 
colnames(occ); colnames(surv_abs)
# Occurrence data
occ_simple <- occ %>%
  dplyr::select(easting, northing, source, occ, geometry)

# Survey absences 
surv_abs_simple <- sf::st_transform(surv_abs, bng_epsg)
surv_abs_simple <- surv_abs_simple %>%
  dplyr::select(easting, northing, ptes, occ, geometry) %>%
  dplyr::rename(source = ptes)

# Combine 
occ_all <- rbind(occ_simple, surv_abs_simple)
# add unique ID number -> row number
occ_all <- occ_all %>%
  dplyr::mutate(id=row_number()) %>%
  dplyr::relocate(id, .before = easting)
# st_write(occ_all, paste0(wd.survey, "occurrence_all_combined.gpkg"), delete_dsn = TRUE)
occ_all <- sf::st_read(paste0(wd.survey, "occurrence_all_combined.gpkg"))

# Sample these points from a raster template for region 
# Make a raster with values only in cells with an occurrence record
rast_tmp <- env_r$anc_wood
# convert all values to 1
m <- matrix(c(-1,Inf,1), ncol=3, byrow=T) 
rast_tmp <- terra::classify(rast_tmp, m) ## entire region is assigned 1
names(rast_tmp) <- "template_100m"
terra::writeRaster(rast_tmp, paste0(wd.survey, "template_100m.tif"), overwrite = TRUE)

# mask with occurrence records --> keep pixels that overlap with points 
occ_rast <- terra::mask(rast_tmp, occ_all)
names(occ_rast) <- "occ_records"
terra::writeRaster(occ_rast, paste0(wd.survey, "occurrence_points.tif"), overwrite = T)
# raster with value 1 for all places with an occurrence record
occ_rast <- terra::rast(paste0(wd.survey, "occurrence_points.tif"))

# Remove locations with occurrence records 
ref_region_survey <- terra::mask(th02_scr006, occ_rast, inverse = TRUE)
plot(ref_region_survey)

# Remove urban areas 
ref_region_survey <- terra::mask(ref_region_survey, urban, inverse = TRUE)


# Remove buffer around Counties
# Dissolve county borders into single polygon
dd_combine <- spatialEco::sf_dissolve(dd_county)
st_write(dd_combine, "gis/os_boundary/devdor_dissolved.shp")
# crop study area to county border
ref_region_survey <- terra::mask(ref_region_survey, dd_combine)

# Change all non-zero values to 1 within reference region
m <- matrix(c(0,Inf,1), ncol=3, byrow=T)
ref_region_survey <- terra::classify(ref_region_survey, m) ## Make all values 1 (excluding 0)
names(ref_region_survey) <- "reference_region_survey"
# terra::writeRaster(ref_region_survey, paste0(wd.survey, "reference_region_survey.tif"), overwrite = T)
ref_region_survey <- terra::rast(paste0(wd.survey, "reference_region_survey.tif"))



#### Treehedge ####
# Targeting treehedge, use greater than 0.02 
# When in an interaction, lower treehedge  is < 0.40
th_gt.02 <- terra::rast(paste0(wd.survey, "treehedge_above_02.tif"))

# Trim to reference region
treehedge_survey_vals <- terra::mask(th_gt.02, ref_region_survey)
names(treehedge_survey_vals) <- "treehedge_survey_vals"
terra::writeRaster(treehedge_survey_vals, paste0(wd.survey, "treehedge_survey_vals.tif"), overwrite = TRUE)
# Change all non-zero values to 1 and change 0 to NA --> keep values only in sample-able cells
m <- matrix(c(0,Inf,1,
              -1,0,NA), 
            nrow = 2, ncol=3, byrow=T)
treehedge_survey <- terra::classify(treehedge_survey_vals, m)
names(treehedge_survey) <- "treehedge_survey"
terra::writeRaster(treehedge_survey, paste0(wd.survey, "treehedge_survey.tif"), overwrite = TRUE)


# Split into low and high ranges --> for interactions
# Intervals 
th_vals <- env_vals_raw[!env_vals_raw$treehedge < 0.02, "treehedge"]
treehedge_int <- classInt::classIntervals(th_vals, n = 3, style = "equal", intervalClosure = "left")
treehedge_int$brks # 0.02000024 0.34661970 0.67323916 0.99985862

# Apply thirds to data w/ treehedge > 0.02
# Remove values in middle interval from already truncated data
m <- matrix(c(treehedge_int$brks[2], treehedge_int$brks[3], NA),
            ncol = 3, byrow = T)
treehedge_low_up_third <- terra::classify(th_gt.02, m, right = F)
names(treehedge_low_up_third) <- "treehedge_up_low_third"
# Trim to reference region
treehedge_low_up_third <- terra::mask(treehedge_low_up_third, ref_region_survey)
terra::writeRaster(treehedge_low_up_third, paste0(wd.survey, "treehedge_up_low_third.tif"), overwrite = TRUE)
treehedge_low_up_third <- terra::rast(paste0(wd.survey, "treehedge_up_low_third.tif"))

# Split lower and upper regions 
treehedge_low_survey <- terra::clamp(treehedge_low_up_third, upper = treehedge_int$brks[2], values = FALSE)
treehedge_high_survey <- terra::clamp(treehedge_low_up_third, lower = treehedge_int$brks[3], values = FALSE)
# Change all non-zero values to 1 and change 0 to NA --> keep values only in sample-able cells
m <- matrix(c(0,Inf,1,
              -1,0,NA), 
            nrow = 2, ncol=3, byrow=T)
treehedge_low_survey <- terra::classify(treehedge_low_survey, m)
treehedge_high_survey <- terra::classify(treehedge_high_survey, m)
names(treehedge_low_survey) <- "treehedge_low_survey"
names(treehedge_high_survey) <- "treehedge_high_survey"
terra::writeRaster(treehedge_low_survey, paste0(wd.survey, "treehedge_low_survey.tif"), overwrite = TRUE)
terra::writeRaster(treehedge_high_survey, paste0(wd.survey, "treehedge_high_survey.tif"), overwrite = TRUE)



#### Scrub ####
# Target all high scrub area; scrub >0.007
scrub_007 <- terra::clamp(scrub, lower = 0.007, values = F) 
# terra::writeRaster(scrub_007, paste0(wd.survey, "scrub_007.tif"), overwrite = TRUE)
scrub_007 <- terra::rast(paste0(wd.survey, "scrub_007.tif"))

# Trim to reference region
scrub_survey_vals <- terra::mask(scrub_007, ref_region_survey)
names(scrub_survey_vals) <- "scrub_survey_vals"
terra::writeRaster(scrub_survey_vals, paste0(wd.survey, "scrub_survey_vals.tif"), overwrite = TRUE)
# Change all non-zero values to 1 and change 0 to NA --> keep values only in sample-able cells
m <- matrix(c(0,Inf,1,
              -1,0,NA), 
            nrow = 2, ncol=3, byrow=T)
scrub_survey <- terra::classify(scrub_survey_vals, m)
names(scrub_survey) <- "scrub_survey"
terra::writeRaster(scrub_survey, paste0(wd.survey, "scrub_survey.tif"), overwrite = TRUE)



####  Ancient woodland ####
# Survey cells with cover >0 but exclude are in an NDMP site 

# Ancient woodland outside NDMP areas
anc_wood_no_ndmp <- terra::mask(env_r$anc_wood, ndmp_bounds, inverse = TRUE) # inverse preserves the non-overlapping regions
# terra::writeRaster(anc_wood_no_ndmp, paste0(wd.survey, "anc_wood_no_ndmp.tif"), overwrite = TRUE)
anc_wood_no_ndmp <- terra::rast(paste0(wd.survey, "anc_wood_no_ndmp.tif"))

# trim to reference region 
anc_wood_survey_vals <- terra::mask(anc_wood_no_ndmp, ref_region_survey) # keep overlapping regions
names(anc_wood_survey_vals) <- "anc_wood_survey_vals"
terra::writeRaster(anc_wood_survey_vals, paste0(wd.survey, "anc_wood_survey_vals.tif"), overwrite = TRUE)
# anc_wood_survey_vals <- terra::rast(paste0(wd.survey, "anc_wood_survey_vals.tif"))

# Change all non-zero values to 1 and change 0 to NA --> keep values only in sample-able cells
m <- matrix(c(0,Inf,1,
              -1,0,NA), 
            nrow = 2, ncol=3, byrow=T)
anc_wood_survey <- terra::classify(anc_wood_survey_vals, m)
names(anc_wood_survey) <- "anc_wood_survey"
terra::writeRaster(anc_wood_survey, paste0(wd.survey, "anc_wood_survey.tif"), overwrite = TRUE)


# Area of woodland outside NDMP
# remove 0 and calculate area of anc wood without NDMP
aw_no_ndmp_rm0 <- terra::subst(anc_wood_no_ndmp, 0, NA)
area_aw_no_ndmp_rm0 <- terra::expanse(aw_no_ndmp_rm0, unit = "km") 
# 1032.193 km^2 woodland outside NDMP sites 

# Ancient woodland area in NDMP
anc_wood_ndmp <- terra::mask(env_r$anc_wood, ndmp_bounds)
plot(anc_wood_ndmp)
terra::writeRaster(anc_wood_ndmp, paste0(wd.survey, "anc_wood_ndmp.tif"), overwrite = TRUE)

cellSize(anc_wood_ndmp, unit = "km") # each cell is ~0.01km^2 (10,000m^2)
area_aw_ndmp <- terra::expanse(anc_wood_ndmp, unit = "km") # summed area for all non-NA raster cells 
# 14.55471 km^2
# remove zeros, calculate area
aw_ndmp_rm0 <- terra::subst(anc_wood_ndmp, 0, NA)
area_aw_ndmp_rm0 <- terra::expanse(aw_ndmp_rm0, unit = "km") 
# 7.202175 km^2

# Total anc_wood area (non-zero)
anc_wood <- env_r$anc_wood
anc_wood_rm0 <- terra::subst(anc_wood, 0, NA)
area_aw_rm0 <- terra::expanse(anc_wood_rm0, unit = "km")
# 1039.395 km^2

# histogram non-zero values 
min(vals_true$value[vals_true$variable == "anc_wood"]) # 1.990403e-24
# majority are very small values
ggplot(subset(vals_true, variable %in% "anc_wood"), 
       aes(x = value)) + 
  geom_histogram(position="identity") + 
  xlab("ancient woodland")  +
  theme_classic() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))



#### Max spring rainfall ####
# Use most recent 5 year mean --> 2022 layer
summary(env_vals_raw$rainfall_spr_5yrmn, na.rm = T); summary(env_vals_raw$rainfall_spr_5yrmn_2022, na.rm = T)
# 2022 mean has lower min and max values

# Remove values equal to and above 350mm - not being considered 
rain_spr_survey_vals <- terra::clamp(env_r$rainfall_spr_5yrmn_2022, upper = 350, values = FALSE)

# Find upper and lower thirds of range being considered
# Remove values >350 from raster val data frame 
rain_vals <- env_vals_raw[!env_vals_raw$rainfall_spr_5yrmn_2022 > 350, "rainfall_spr_5yrmn_2022"]
summary(rain_vals, na.rm= T)
rain_spr_int <- classInt::classIntervals(rain_vals, n = 3, style = "equal", intervalClosure = "left")
rain_spr_int$brks # 128.6959 202.4637 276.2316 349.9995

# Remove values in middle interval from already truncated data
m <- matrix(c(rain_spr_int$brks[2], rain_spr_int$brks[3], NA),
            ncol = 3, byrow = T)
rain_spr_survey_vals <- terra::classify(rain_spr_survey_vals, m, right = FALSE) # right = FALSE means interval is open on right and closed on left, as it is in the interval we're calling from
rain_spr_survey_vals <- terra::crop(rain_spr_survey_vals, devdor_buff, mask = TRUE)
terra::writeRaster(rain_spr_survey_vals, paste0(wd.survey, "rainfall_spr_5yrmn_low_up_third.tif"), overwrite = TRUE) # raster with lower and upper third (truncated to below 350mm)

# Trim to reference region
rain_spr_survey_vals <- terra::mask(rain_spr_survey_vals, ref_region_survey)
names(rain_spr_survey_vals) <- "rainfall_spr_5yrmn_survey_vals"
terra::writeRaster(rain_spr_survey_vals, paste0(wd.survey, "rainfall_spr_5yrmn_survey_vals.tif"), overwrite = TRUE)
# Convert all values to 1
m <- matrix(c(0,Inf,1), ncol=3, byrow=T)
rain_spr_survey <- terra::classify(rain_spr_survey_vals, m)
names(rain_spr_survey) <- "rainfall_spr_5yrmn_survey"
terra::writeRaster(rain_spr_survey, paste0(wd.survey, "rainfall_spr_5yrmn_survey.tif"), overwrite = TRUE)

# Lower and upper survey regions individually 
rain_spr_low_survey_vals <- terra::clamp(rain_spr_survey_vals, upper = rain_spr_int$brks[2], values = FALSE)
rain_spr_high_survey_vals <- terra::clamp(rain_spr_survey_vals, lower = rain_spr_int$brks[3], values = FALSE)
# Convert all values to 1
m <- matrix(c(0,Inf,1), ncol=3, byrow=T)
rain_spr_low_survey <- terra::classify(rain_spr_low_survey_vals, m)
rain_spr_high_survey <- terra::classify(rain_spr_high_survey_vals, m)
writeRaster(rain_spr_low_survey, paste0(wd.survey, "rainfall_spr_5yrmn_low_survey.tif"), overwrite = TRUE)
writeRaster(rain_spr_high_survey, paste0(wd.survey, "rainfall_spr_5yrmn_high_survey.tif"), overwrite = TRUE)

rain_spr_low_survey <- terra::rast(paste0(wd.survey, "rainfall_spr_5yrmn_low_survey.tif"))
rain_spr_high_survey <- terra::rast(paste0(wd.survey, "rainfall_spr_5yrmn_high_survey.tif"))


# Destandardize values 
std <- seq(-1, 5, length.out=7) 
destd <- (std * destdize[2,"rainfall_spr_5yrmn"]) + destdize[1,"rainfall_spr_5yrmn"] ## 2022 mean isn't in destandardize table so just check with original file
# -1  0  1  2  3  4  5
# 160.7303 206.1920 251.6536 297.1153 342.5770 388.0386 433.5003



#### Min winter temperature ####
# Use most recent 5 year mean --> 2022 layer
summary(env_vals_raw$tasmin_win_5yrmn, na.rm = T); summary(env_vals_raw$tasmin_win_5yrmn_2022, na.rm = T)
# 2022 mean has slightly lower min and max values 

# Find upper and lower third 
min_win_t_int <- classInt::classIntervals(env_vals_raw$tasmin_win_5yrmn_2022, n = 3, style = "equal", intervalClosure = "left")
min_win_t_int$brks # 1.242301 2.961575 4.680849 6.400123

# Remove values in middle interval 
m <- matrix(c(min_win_t_int$brks[2], min_win_t_int$brks[3], NA),
            ncol = 3, byrow = T)
min_win_t_survey_vals <- terra::classify(env_r$tasmin_win_5yrmn_2022, m, right = FALSE) # right = FALSE means interval is open on right and closed on left, as it is in the interval we're calling from
min_win_t_survey_vals <- terra::crop(min_win_t_survey_vals, devdor_buff, mask = TRUE)
terra::writeRaster(min_win_t_survey_vals, paste0(wd.survey, "tasmin_win_5yrmn_low_up_third.tif"), overwrite = TRUE) # raster with lower and upper third

# Trim to region
min_win_t_survey_vals <- terra::mask(min_win_t_survey_vals, ref_region_survey)
names(min_win_t_survey_vals) <- "tasmin_win_5yrmn_survey_vals"
terra::writeRaster(min_win_t_survey_vals, paste0(wd.survey, "tasmin_win_5yrmn_survey_vals.tif"), overwrite = TRUE)
# Convert all values to 1
m <- matrix(c(0,Inf,1), ncol=3, byrow=T)
min_win_t_survey <- terra::classify(min_win_t_survey_vals, m)
names(min_win_t_survey) <- "tasmin_win_5yrmn_survey"
terra::writeRaster(min_win_t_survey, paste0(wd.survey, "tasmin_win_5yrmn_survey.tif"), overwrite = TRUE)

# Lower and upper survey regions individually 
min_win_t_low_survey_vals <- terra::clamp(env_r$tasmin_win_5yrmn_2022, upper = min_win_t_int$brks[2], values = FALSE)
min_win_t_high_survey_vals <- terra::clamp(env_r$tasmin_win_5yrmn_2022, lower = min_win_t_int$brks[3], values = FALSE)
# Convert all values to 1
m <- matrix(c(0,Inf,1), ncol=3, byrow=T)
min_win_t_low_survey <- terra::classify(min_win_t_low_survey_vals, m)
min_win_t_high_survey <- terra::classify(min_win_t_high_survey_vals, m)
names(min_win_t_low_survey) <- "tasmin_win_5yrmn_low_survey"
names(min_win_t_high_survey) <- "tasmin_win_5yrmn_high_survey"
writeRaster(min_win_t_low_survey, paste0(wd.survey, "tasmin_win_5yrmn_low_survey.tif"), overwrite = TRUE)
writeRaster(min_win_t_high_survey, paste0(wd.survey, "tasmin_win_5yrmn_high_survey.tif"), overwrite = TRUE)

min_win_t_low_survey <- terra::rast(paste0(wd.survey, "tasmin_win_5yrmn_low_survey.tif"))
min_win_t_high_survey <- terra::rast(paste0(wd.survey, "tasmin_win_5yrmn_high_survey.tif"))


#### Max spring temperature ####
# Use most recent 5 year mean --> 2022 layer
summary(env_vals_raw$tasmax_spr_5yrmn, na.rm = T); summary(env_vals_raw$tasmax_spr_5yrmn_2022, na.rm = T)
# 2022 mean has slightly higher min and max temperature 

# Target restricted to above 14degC
max_spr_t_survey_vals <- terra::clamp(env_r$tasmax_spr_5yrmn_2022, lower = 14, values = FALSE)
max_spr_t_survey_vals <- terra::crop(max_spr_t_survey_vals, devdor_buff, mask = TRUE)
plot(max_spr_t_survey_vals)
# terra::writeRaster(max_spr_t_survey_vals, paste0(wd.survey, "tasmax_spr_5yrmn_above14.tif"), overwrite = TRUE)

# Trim to reference region
max_spr_t_survey_vals <- terra::mask(max_spr_t_survey_vals, ref_region_survey)
names(max_spr_t_survey_vals) <- "tasmax_spr_5yrmn_survey_vals"
terra::writeRaster(max_spr_t_survey_vals, paste0(wd.survey, "tasmax_spr_5yrmn_survey_vals.tif"), overwrite = TRUE)
# Convert all values to 1
m <- matrix(c(0,Inf,1), ncol=3, byrow=T)
max_spr_t_survey <- terra::classify(max_spr_t_survey_vals, m)
names(max_spr_t_survey) <- "tasmax_spr_5yrmn_survey"
terra::writeRaster(max_spr_t_survey, paste0(wd.survey, "tasmax_spr_5yrmn_survey.tif"), overwrite = TRUE)


# > Explore data ####
# Min winter temperature had different effects last project between NDMP and NDD sites 
dat <- read.csv(here("data/env_occ_clean.csv")) 
# Focus on data source 
dat.source <- dat[, c("source",vars)]
dat.source <- dat.source %>%
  gather(-source, key = "var", value = "value") 

dat.source <- dat.source[dat.source$source!= "pseudo-absence", ]

min(dat.source$value[dat.source$source == "NDMP" & dat.source$var == "tasmin_win_5yrmn"]) #1.441857
summary(dat.source$value[dat.source$source == "NDMP" & dat.source$var == "tasmin_win_5yrmn"]) 
min(dat.source$value[dat.source$source == "LRD+NDD" & dat.source$var == "tasmin_win_5yrmn"]) # 0.5048469
summary(dat.source$value[dat.source$source == "LRD+NDD" & dat.source$var == "tasmin_win_5yrmn"])
## similar range, median, and Qs between both
temp <- dat.source[dat.source$source == "NDMP",]
nrow(temp[temp$var == "tasmin_win_5yrmn" & temp$value > 3.5,]) # 18 NDMP points in high min win temps
temp <- dat.source[dat.source$source == "LRD+NDD",]
nrow(temp[temp$var == "tasmin_win_5yrmn" & temp$value > 3.5,]) # 120 LRD points in high min win temps

min(env_vals_raw$tasmin_win_5yrmn, na.rm = T)
max(env_vals_raw$tasmin_win_5yrmn, na.rm = T)



#### Sun spring ####
# Use most recent 5 year mean --> 2022 layer
summary(env_vals_raw$sun_spr_5yrmn, na.rm = T); summary(env_vals_raw$sun_spr_5yrmn_2022, na.rm = T)
# 2022 mean has high min and max values 

# Find upper and lower thirds 
sun_spr_int <- classInt::classIntervals(env_vals_raw$sun_spr_5yrmn_2022, n = 3, style = "equal", intervalClosure = "left")
sun_spr_int$brks #  156.1325 181.9197 207.7070 233.4943

# Remove values in middle interval
m <- matrix(c(sun_spr_int$brks[2], sun_spr_int$brks[3], NA),
            ncol = 3, byrow = T)
sun_spr_survey_vals <- terra::classify(env_r$sun_spr_5yrmn_2022, m, right = FALSE) # right = FALSE means interval is open on right and closed on left, as it is in the interval we're calling from
sun_spr_survey_vals <- terra::crop(sun_spr_survey_vals, devdor_buff, mask = TRUE)
# terra::writeRaster(sun_spr_survey_vals, paste0(wd.survey, "sun_spr_5yrmn_low_up_third.tif"), overwrite = TRUE)
sun_spr_survey_vals <- terra::rast(paste0(wd.survey, "sun_spr_5yrmn_low_up_third.tif"))

# Trim to reference region
sun_spr_survey_vals <- terra::mask(sun_spr_survey_vals, ref_region_survey)
names(sun_spr_survey_vals) <- "sun_spr_5yrmn_survey_vals"
terra::writeRaster(sun_spr_survey_vals, paste0(wd.survey, "sun_spr_5yrmn_survey_vals.tif"), overwrite = TRUE)
# Convert all values to 1
m <- matrix(c(0,Inf,1), ncol=3, byrow=T)
sun_spr_survey <- terra::classify(sun_spr_survey_vals, m)
names(sun_spr_survey) <- "sun_spr_5yrmn_survey"
terra::writeRaster(sun_spr_survey, paste0(wd.survey, "sun_spr_5yrmn_survey.tif"), overwrite = TRUE)

# Lower and upper survey regions individually 
sun_spr_low_survey_vals <- terra::clamp(env_r$sun_spr_5yrmn_2022, upper = sun_spr_int$brks[2], values = FALSE)
sun_spr_high_survey_vals <- terra::clamp(env_r$sun_spr_5yrmn_2022, lower = sun_spr_int$brks[3], values = FALSE)
# Convert all values to 1
m <- matrix(c(0,Inf,1), ncol=3, byrow=T)
sun_spr_low_survey <- terra::classify(sun_spr_low_survey_vals, m)
names(sun_spr_low_survey) <- "sun_spr_low_survey"
sun_spr_high_survey <- terra::classify(sun_spr_high_survey_vals, m)
names(sun_spr_high_survey) <- "sun_spr_high_survey"

terra::writeRaster(sun_spr_low_survey, paste0(wd.survey, "sun_spr_5yrmn_low_survey.tif"), overwrite = TRUE)
terra::writeRaster(sun_spr_high_survey, paste0(wd.survey, "sun_spr_5yrmn_high_survey.tif"), overwrite = TRUE)



#### Conifer Cover (mean) ####
# All areas with conifer cover >0
# High cover (upper third) used for interactions

# Trim to reference region 
conif_survey_vals <- terra::mask(env_r$combo_conif_m, ref_region_survey) # keep overlapping regions
names(conif_survey_vals) <- "combo_conif_m_survey_vals"
terra::writeRaster(conif_survey_vals, paste0(wd.survey, "combo_conif_m_survey_vals.tif"), overwrite = TRUE) # still contains 0s
#conif_survey_vals <- terra::rast(paste0(wd.survey, "combo_conif_m_survey_vals.tif"))

# Change all non-zero values to 1 and change 0 to NA --> keep values only in sample-able cells
m <- matrix(c(0,Inf,1,
              -1,0,NA), 
            nrow = 2, ncol=3, byrow=T)
conif_survey <- terra::classify(conif_survey_vals, m)
names(conif_survey) <- "combo_conif_m_survey"
terra::writeRaster(conif_survey, paste0(wd.survey, "combo_conif_m_survey.tif"), overwrite = TRUE)

# Find upper and lower thirds 
conif_int <- classInt::classIntervals(env_vals_raw$combo_conif_m, n = 3, style = "equal", intervalClosure = "left")
conif_int$brks # 0.0000 0.3125 0.6250 0.9375

# Upper survey region 
conif_high_survey_vals <- terra::clamp(env_r$combo_conif_m, lower = conif_int$brks[3], values = FALSE)
# convert all non-zero values to 1 
m <- matrix(c(0,Inf,1), ncol=3, byrow=T)
conif_high_survey <- terra::classify(conif_high_survey_vals, m)
names(conif_high_survey) <- "combo_conif_m_high_survey"
terra::writeRaster(conif_high_survey, paste0(wd.survey, "combo_conif_m_high_survey.tif"), overwrite = TRUE)




#### 4. EXPLORE COMBINATIONS OF VARIABLE SUBSETS ####

#### Minimum Winter Temperature pairs ####
# min_win_t_low_survey and min_win_t_high_survey
# > Scrub High ####
# -> the base scrub survey region is only focused on high scrub 
min_win_H_scrub_H <-  terra::mask(min_win_t_high_survey, scrub_survey)
#plot(min_win_H_scrub_H)
min_win_L_scrub_H <-  terra::mask(min_win_t_low_survey, scrub_survey)
#plot(min_win_L_scrub_H)
a <- terra::expanse(min_win_H_scrub_H, unit = "km") # 11.3 km^2
b <- terra::expanse(min_win_L_scrub_H, unit = "km") # 134.1
terra::writeRaster(min_win_H_scrub_H, paste0(wd.survey, "min_win_temp_HIGH_scrub_HIGH.tif"), overwrite = TRUE)
terra::writeRaster(min_win_L_scrub_H, paste0(wd.survey, "min_win_temp_LOW_scrub_HIGH.tif"), overwrite = TRUE)

# > Treehedge Low ####
# treehedge_low_survey
min_win_H_treehedge_L <- terra::mask(min_win_t_high_survey, treehedge_low_survey)
min_win_L_treehedge_L <- terra::mask(min_win_t_low_survey, treehedge_low_survey)
a <- terra::expanse(min_win_H_treehedge_L, unit = "km") # 83.2
b <- terra::expanse(min_win_L_treehedge_L, unit = "km") # 3019.9
terra::writeRaster(min_win_H_treehedge_L, paste0(wd.survey, "min_win_temp_HIGH_treehedge_LOW.tif"), overwrite = T)
terra::writeRaster(min_win_L_treehedge_L, paste0(wd.survey, "min_win_temp_LOW_treehedge_LOW.tif"), overwrite = T)


#### Spring Sunshine ####
# sun_spr_survey -> upper and lower thirds 
# sun_spr_low_survey; sun_spr_high_survey
# > Scrub ####
sun_spr_L_scrub_H <- terra::mask(sun_spr_low_survey, scrub_survey)
sun_spr_H_scrub_H <- terra::mask(sun_spr_high_survey, scrub_survey)
a <- terra::expanse(sun_spr_L_scrub_H, unit = "km") # 164.7
b <- terra::expanse(sun_spr_H_scrub_H, unit = "km") # 19.0
terra::writeRaster(sun_spr_L_scrub_H, paste0(wd.survey, "sun_spr_LOW_scrub_HIGH.tif"), overwrite = T)
terra::writeRaster(sun_spr_H_scrub_H, paste0(wd.survey, "sun_spr_HIGH_scrub_HIGH.tif"), overwrite = T)

# > Spring Temperature ####
# max_spr_t_survey
sun_spr_L_max_spr_t_H <- terra::mask(sun_spr_low_survey, max_spr_t_survey)
sun_spr_H_max_spr_t_H <- terra::mask(sun_spr_high_survey, max_spr_t_survey)
a <- terra::expanse(sun_spr_L_max_spr_t_H, unit = "km") # 993.6
b <- terra::expanse(sun_spr_H_max_spr_t_H, unit = "km") # 90.8
terra::writeRaster(sun_spr_L_max_spr_t_H, paste0(wd.survey, "sun_spr_LOW_max_spr_t_HIGH.tif"), overwrite = T)
terra::writeRaster(sun_spr_H_max_spr_t_H, paste0(wd.survey, "sun_spr_HIGH_max_spr_t_HIGH.tif"), overwrite = T)

# > Conifer Cover ####
# conif_high_survey
sun_spr_L_conif_H <- terra::mask(sun_spr_low_survey, conif_high_survey)
sun_spr_H_conif_H <- terra::mask(sun_spr_high_survey, conif_high_survey)
a <- terra::expanse(sun_spr_L_conif_H, unit = "km") # 91.1
b <- terra::expanse(sun_spr_H_conif_H, unit = "km") # 5.2
terra::writeRaster(sun_spr_L_conif_H, paste0(wd.survey, "sun_spr_LOW_conif_HIGH.tif"), overwrite = T)
terra::writeRaster(sun_spr_H_conif_H, paste0(wd.survey, "sun_spr_HIGH_conif_HIGH.tif"), overwrite = T)



#### Spring Rainfall ####
# rain_spr_low_survey; rain_spr_high_survey
# rain_spr_survey

# > Treehedge ####
rain_spr_H_treehedge_L <- terra::mask(rain_spr_high_survey, treehedge_low_survey)
rain_spr_H_treehedge_H <- terra::mask(rain_spr_high_survey, treehedge_high_survey)
a <- terra::expanse(rain_spr_H_treehedge_L, unit = "km") # 318.3
b <- terra::expanse(rain_spr_H_treehedge_H, unit = "km") # 31.0
terra::writeRaster(rain_spr_H_treehedge_L, paste0(wd.survey, "rain_spr_HIGH_treehedge_LOW.tif"), overwrite = T)
terra::writeRaster(rain_spr_H_treehedge_H, paste0(wd.survey, "rain_spr_HIGH_treehedge_HIGH.tif"), overwrite = T)

# > Conifer cover ####
rain_spr_L_conif_H <- terra::mask(rain_spr_low_survey, conif_high_survey)
rain_spr_H_conif_H <- terra::mask(rain_spr_high_survey, conif_high_survey)
a <- terra::expanse(rain_spr_L_conif_H, unit = "km") # 94.1
b <- terra::expanse(rain_spr_H_conif_H, unit = "km") # 9.8
terra::writeRaster(rain_spr_L_conif_H, paste0(wd.survey, "rain_spr_LOW_conif_HIGH.tif"), overwrite = T)
terra::writeRaster(rain_spr_H_conif_H, paste0(wd.survey, "rain_spr_HIGH_conif_HIGH.tif"), overwrite = T)

