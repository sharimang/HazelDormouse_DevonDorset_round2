#### Evaluate Model ####

## Plot variable response curves ##


## Written by: Shari Mang ##
## Date: December 2023 ##
## Based on code by: Regan Early, written in 2021 ##



#### Set up ####
library(ggplot2) 
library(interactions) ## Needed to make interact_plot run and needs to be version for R 3.3.0 or later
library(effects)
library(terra)


wd.out <- here("model/thaw_pseudo_abs/")
wd.plot <- here("model/thaw_pseudo_abs/plots/")

mod.name <- "model_pseudo_allpres_01_"


#### Data ####
destdize <- read.csv(paste0(wd.out, "model_pseudo_allpres_01_dat_means_sd_outlier_rm.csv"))
# row 1 = mean, row 2 = SD
# Load standardized data
dat_stnd <- read.csv(paste0(wd.out, mod.name, "dat_stnd_scrub_outliers_removed.csv")) # equal to Regan's dat.m
dat.m <- dat_stnd
# load final averaged model
load(file=paste0(wd.out, mod.name,"final_averaged_thawfix")) # final_m01
final_m01_coefs <- read.csv(paste0(wd.out, mod.name,"final_averaged_thawfix.csv"))
# Data used in model
dat.c <- read.csv(paste0(wd.out, mod.name, "dat_stnd_cal.csv"))

vars <- as.vector(colnames(final_m01$coefficients))
vars <- vars[!grepl("Intercept|I|:", vars)] 

var_names <- c("Mean conifer cover", "Aspect WE",
               "Mean spring rainfall", "Scrub",
               "Mean spring sun", "Max spring temp", 
               "Treehedge", "Slope",
               "Ancient woodland", "Change broadleaf cover",
               "Aspect NS", "Min winter temp")

mns <- apply(dat.m[,vars], 2, mean) 


#### Graph response curves ####

# Run through all variables.
# Use standardized data to use for predictions
# can't extract family from averaged model -- pull from glm that went into final dredge (or any with same family)
m <- glm(as.formula(paste0("occ ~ ", vars)), family=binomial, dat.c, na.action=na.fail) 
ilink <- family(m)$linkinv # get the link function from the model

# y.limits for graphs
ylims <- rep(1.0, length(vars))
names(ylims) <- vars

png(paste0(wd.plot, mod.name, "final_response_curves.png"), width=700, height=500)
par(mfrow=c(3,4))
par(mar=c(4,4.5,2,1)) # c(bottom, left, top, right)

i <- vars[1]
for(i in vars) { 
  newdat <- data.frame(matrix(data=mns, 
                              byrow=T, nrow=10000, ncol=length(vars), dimnames=(list(NULL,vars)))) ## Make a dataframe with mean (standardised) values of each variable
  newdat[,i] <- seq(min(dat.m[,i]), max(dat.m[,i]), length.out=10000) #
  var_df <- as.data.frame(subset(newdat, select = i)) # isolate variable of interest to bind to prediction output
  
  pdat <- predict(final_m01, newdat, type = "link", se.fit = TRUE) # Predict on the link scale
  pdat <- cbind(var_df, predict(final_m01, newdat, type = "link", se.fit = TRUE)[1:2])
  pdat <- transform(pdat, Fitted = ilink(fit), Upper = ilink(fit + (2 * se.fit)), ## Backtransform to the probability scale using the inverse of the link function
                    Lower = ilink(fit - (2 * se.fit)))
  
  ### Calculate the raw (destandardised) values of the explanatory variable
  ## Multiply by standard deviation (row 2 of destdize) and add mean (row 1 of destdize)
  x.std <- seq(min(newdat[,i]), max(newdat[,i]), length.out=6) ## Identify six values for plotting
  x.destd <- round((x.std * destdize[2,i]) + destdize[1,i], 2) ## Multiply standardised value by sd of original data and add the mean
  
  # if(i %in% "scrub") { ## If variable is scrub, multiply by 100 as values are so small axis is all 0
  #   x.destd <- round(((x.std * destdize[2,i]) + destdize[1,i])*100, 1) ## scrub value is multiplied by 100
  # } else {
  #   x.destd <- round((x.std * destdize[2,i]) + destdize[1,i], 1) ## Multiply standardised value by sd of original data and add the mean
  # }  
  
  plot(newdat[,i], pdat$Fitted, type="n", xlab=var_names[vars==i],
       ylab="Prob. occurrence", xaxt="n", yaxt="n",
       ylim=c(0,ylims[i]),
       cex.lab=1.5) 
  axis(side=1, at=x.std, labels=x.destd, cex.axis=1.2) ## Add x axis on original scale
  axis(side=2, at=seq(0,1.5,0.3), cex.axis=1.2) ##, labels=c(0,format(ylims[i], scientific = FALSE)))
  l.ci <- pdat$Lower
  u.ci <- pdat$Upper
  polygon(c(newdat[,i], rev(newdat[,i])), c(c(l.ci, rev(u.ci))), col = 'grey80', border = NA)
  lines(newdat[,i], pdat$Fitted)
}
dev.off()
## Note: using type = "response" is similar to transforming fit when predicted with type = "link"
## Difference between "response" fit and "link" Fitted column ~ 0.00035




##### Interaction plots ##### 
# Can't enter a model average object into interact_plot so make a dummy model
# Make a dummy model to alter coefficients
var.inter <- as.character(final_m01_coefs$Variables, as.factor=F)
var.inter <-subset(var.inter, !grepl("Intercept", var.inter))
var.inter <- paste0(var.inter, collapse=" + ")
dummy <- glm(as.formula(paste0("occ ~ ", var.inter)), family=binomial, dat.m, na.action=na.fail) ## Make model with calibration data

## Update dummy model with coefficients from model average
coef <- names(dummy$coefficients)
coef[26] <- final_m01_coefs$Variables[26] ## for some reason the variable position in interaction is swapped
n <- 1
for(n in 1:length(coef)) {
  dummy$coefficients[coef[n]] <- final_m01_coefs[final_m01_coefs$Variables==coef[n],"MeanCoefficients"]
}


## Identify the interactors
var.inter <- as.character(final_m01_coefs$Variables, as.factor=F)
interactors <- var.inter[grepl(":", var.inter)]
interactors <- strsplit(interactors, ":") ## 9 interactors

## Remember can't trust SD because taken from original dummy model
## run manually, as interactors[[1]][1] (etc) doesn't seem to work

## Start with interactions in all of the best models:
jpeg(paste0(wd.plot, mod.name, "final_response_int_ConifM-RainSpr.png"), width=600, height=400)
interact_plot(dummy, pred=combo_conif_m, modx=rainfall_spr_5yrmn, interval=T, rug=T)
dev.off()

jpeg(paste0(wd.plot, mod.name, "final_response_int_ConifM-SunSpr.png"), width=600, height=400)
interact_plot(dummy, pred=combo_conif_m, modx=sun_spr_5yrmn, interval=T, rug=T)
dev.off()

jpeg(paste0(wd.plot, mod.name, "final_response_int_AspectWE-RainSpr.png"), width=600, height=400)
interact_plot(dummy, pred=OS_Terrain_100_WE, modx=rainfall_spr_5yrmn, interval=T, rug=T)
dev.off()

jpeg(paste0(wd.plot, mod.name, "final_response_int_RainSpr-Treehedge.png"), width=600, height=400)
interact_plot(dummy, pred=rainfall_spr_5yrmn, modx=treehedge, interval=T, rug=T)
dev.off()

jpeg(paste0(wd.plot, mod.name, "final_response_int_Scrub-SunSpr.png"), width=600, height=400)
interact_plot(dummy, pred=scrub, modx=sun_spr_5yrmn, interval=T, rug=T) #  modx.values = "terciles", 
dev.off()

jpeg(paste0(wd.plot, mod.name, "final_response_int_SunSpr-MaxSpr.png"), width=600, height=400)
interact_plot(dummy, pred=sun_spr_5yrmn, modx=tasmax_spr_5yrmn,  interval=T, rug=T)
dev.off()

# less important
jpeg(paste0(wd.plot, mod.name, "final_response_int_SunSpr-Treehedge.png"), width=600, height=400)
interact_plot(dummy, pred=sun_spr_5yrmn, modx=treehedge, interval=T, rug=T)
dev.off()

jpeg(paste0(wd.plot, mod.name, "final_response_int_AspectNS-RainSpr.png"), width=600, height=400)
interact_plot(dummy, pred=OS_Terrain_100_NS, modx=rainfall_spr_5yrmn, interval=T, rug=T)
dev.off()

jpeg(paste0(wd.plot, mod.name, "final_response_int_MaxSpr-Treehedge.png"), width=600, height=400)
interact_plot(dummy, pred=tasmax_spr_5yrmn, modx=treehedge, interval=T, rug=T)
dev.off()






#### Explore Standard Error ####
# Final model prediction 
pred <- terra::rast(paste0(wd.out, mod.name,"final_preds_thawfix.tif"))
pred_dat <- as.data.frame(terra::extract(pred, 1:ncell(pred))) # has both prediction and standard error

ci_range <- readRDS(paste0(wd.out, mod.name, "variable_predict_CI_range.rds"))

pred.se <- terra::rast(paste0(wd.out, mod.name, "final_preds_thawfix_se.tiff"))
pred_se_dat <- as.data.frame(terra::extract(pred.se, 1:ncell(pred.se)))
colnames(pred_se_dat)[1] <- "se_predict"

pred_dat$se_predict <- pred_se_dat$se_predict # not clear that there's any difference between se.fit and se from the prediction
pred_dat$check <- pred_dat$se.fit - pred_dat$se_predict
any(pred_dat$check != 0, na.rm = T) # all values are zero, se sources are same

View(pred_dat)
tmp <- pred_dat %>%
  drop_na()
tmp <- tmp[tmp$se.fit > 0.50, ]



# Histogram of SE and prediction values 
pred_hist <- pred_dat %>%
  gather(key = "metric", value = "value")
unique(pred_hist$metric)

ggplot(subset(pred_hist, metric %in% "fit"),
       aes(x = value)) + 
  geom_histogram(position="identity", colour = "black", fill = "lightgrey") + 
  xlab("Model prediction")  +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18)) #+
  # facet_wrap(~ metric) 
min(pred_dat$fit, na.rm = T); max(pred_dat$fit, na.rm = T) # 2.582535e-06, 0.9944764

ggplot(subset(pred_hist, metric %in% "se.fit"),
       aes(x = value)) + 
  geom_histogram(position="identity", colour = "black", fill = "lightgrey") + 
  xlab("Model standard error")  +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18)) #+
# facet_wrap(~ metric) 

min(pred_dat$se.fit, na.rm = T); max(pred_dat$se.fit, na.rm = T) # 8.407495e-06, 0.7913738

# Prediction vs SE
colnames(pred_dat)
ggplot(pred_dat, aes(x = fit, y = se.fit)) + 
  geom_point(size = 2)+
  ylab("Standard Error") + 
  xlab("Model fit") +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) + 
  scale_x_continuous(limits = c(0,1)) + 
  theme_classic() + 
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))



# Compare SE to the CI range calculated from individual variable predictions 
pred_dat <- cbind(pred_dat, ci_range$range_ci_overall)
colnames(pred_dat) [3] <- "range_ci_calculated"
median(pred_dat$se.fit, na.rm = T)
ggplot(pred_dat, aes(x = range_ci_calculated, y = se.fit)) + 
  geom_point(size = 2)+
  ylab("Standard Error") + 
  xlab("Calculated CI Range") +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) + 
  scale_x_continuous() + 
  theme_classic() + 
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))
min(pred_dat$se.fit, na.rm = T); max(pred_dat$se.fit, na.rm = T) # 8.407495e-06, 0.7913738
min(pred_dat$range_ci_calculated, na.rm = T); max(pred_dat$range_ci_calculated, na.rm = T) # 0.1964527; 1.664959



# Compare CI to SE from model with only linear variables from final version of model 
p_se_lin <- terra::rast(paste0(wd.out, mod.name, "finalv4_linear_only_se.tif"))
lin_se_dat <- as.data.frame(terra::extract(p_se_lin, 1:ncell(p_se_lin)))

lin_se_dat <- cbind(lin_se_dat, ci_range$range_ci_overall)
colnames(lin_se_dat)[2] <- "range_ci_calculated"

ggplot(lin_se_dat, aes(x = range_ci_calculated, y = se_linear_only_v4)) + 
  geom_point(size = 2)+
  ylab("Standard Error: linear vars only") + 
  xlab("Calculated CI Range") +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) + 
  scale_x_continuous() + 
  theme_classic() + 
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))

min(lin_se_dat$se_linear_only_v4, na.rm = T); max(lin_se_dat$se_linear_only_v4, na.rm = T) #2.503394e-15; 0.3705535


# as histogram
lin_se_hist <- lin_se_dat %>%
  gather(key = "metric", value = "value")
unique(lin_se_hist$metric)

ggplot(lin_se_hist,
       aes(x = value, colour = metric)) + 
  geom_histogram(position="identity", fill = "white") + 
 # xlab("Model prediction")  +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18), 
        legend.text = element_text(size = 14),
        legend.title = element_text(size=14),
        plot.title = element_text(size = 18), 
        strip.text = element_text(size = 14)) +
  facet_wrap(~ metric) 



# Compare SE linear only to SE full model 
lin_se_dat <- cbind(lin_se_dat, pred_dat$se.fit)
colnames(lin_se_dat)[3] <- "se_full_modv4"


ggplot(lin_se_dat, aes(x = se_full_modv4, y = se_linear_only_v4)) + 
  geom_point(size = 2)+
  ylab("Standard Error: linear vars only") + 
  xlab("Standard Error: Full model") +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) + 
  scale_x_continuous(limits = c(0,1)) + 
  theme_classic() + 
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))



