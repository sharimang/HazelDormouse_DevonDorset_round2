#### Species Distribution Model ####
## Model 01 --> All presence data and pseudo-absences, incorporating THaW data ##
## Developing GLM for the SDM ###

## Written by: Shari Mang ##
## Date: November 2023 ##
## Based on code by: Regan Early, written in 2021 ##


## ~~ Overview ~~ ##
## Include THaW data as treehedge and scrub
## Absences are pseudo-absence 
## All occurrence records used (i.e., NDMP, LRD, and NDD)
## Starting point for variable selection are those that were included in models from previous project




#### SET UP ####
#pacman::p_install(car, snow)
pacman::p_load(dplyr, MuMIn, car, snow, lme4, glmmTMB, DHARMa,here, effects, tidyverse, conflicted)
conflicted::conflict_prefer("here", "here")

bng_epsg <- "EPSG:27700"
bng <- "+init=epsg:27700"

wd.dat <- here("data/")
wd.orig <- here("model/thaw_pseudo_abs/")
wd.out <- here("model/thaw_pseudo_abs/")

# File outputs for this model named: model_pseudo_allpres_01_
mod.name <- "model_pseudo_allpres_01_"


#### Load Data ####
dat <- read.csv(here("data/env_occ_clean.csv")) 
colnames(dat)
any(is.na(dat[, c(8:ncol(dat))])) # model won't run with NAs

vars <- c('combo_broadl_m', 'combo_broadl_c', 'combo_conif_m', 'combo_conif_c', 'anc_wood', 
          'treehedge', 'scrub',
          'OS_Terrain_100_NS', 'OS_Terrain_100_WE','OS_Terrain_100_slope_pct', 
          "tasrng_win_5yrmn", "tasmax_spr_5yrmn", "tasmin_win_5yrmn", 
          "sun_spr_5yrmn", "rainfall_spr_5yrmn")


#### Format data for model ####

# > Remove scrub outliers ####
# Data heavily skewed so removing outliers based on z-score 
# Calculate z score
dat$scrub_zscore <- (abs(dat$scrub-mean(dat$scrub))/sd(dat$scrub))
length(which(dat$scrub_zscore>3)) # 31 occurrences where z-score > 3

# remove rows with a zscore greater than 3
dat <- dat[!dat$scrub_zscore>3, ]
# write.csv(dat, here("data/env_occ_clean_scrub_outliers_rm.csv"))
dat <- read.csv(here("data/env_occ_clean_scrub_outliers_rm.csv"))

dat.plot <- dat[, c("source","scrub")]
dat.plot <- dat.plot %>%
  gather(-source, key = "var", value = "value") 

ggplot(dat.plot, aes(x=as.factor(source), y=value)) +
  geom_boxplot() +
  facet_wrap(~ var, scales = "free") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) 


# > Standardize data ####

# First save the values used to standardize, so as to graph and map later
# For data with scrub outlier removed.
means <- apply(na.omit(dat[,vars]), 2, mean)
sds <- apply(na.omit(dat[,vars]), 2, sd)
destdize <- rbind(means, sds)
write.csv(destdize, paste0(wd.out, "model_pseudo_allpres_01_dat_means_sd_outlier_rm.csv"), row.names=F) 

# Remove NAs
# Shouldn't be any but run just to check
dat_stnd <- dat
any(is.na(dat[, c("occ", vars, "source")])) # No NAs in data
# dat_stnd <- na.omit(dat_stnd[,c("occ", vars, "source")]) 

# Standardize variables
dat_stnd[,vars] <- MuMIn::stdize(dat_stnd[,vars]) 
dat_stnd <- as.data.frame(cbind(dat_stnd, id=c(1:nrow(dat_stnd))))
# write.csv(dat_stnd, paste0(wd.out, mod.name, "dat_stnd_scrub_outliers_removed.csv"), row.names = FALSE)
dat_stnd <- read.csv(paste0(wd.out, mod.name, "dat_stnd_scrub_outliers_removed.csv"))
## ~~~ If including survey absences need to trim data here to have same number of samples ~~~ ##

# > Calibration and validation data #### --> scrub outliers removed
dat.c <- sample_frac(cbind(dat_stnd), size=0.7)
#write.csv(dat.c, paste0(wd.out, mod.name, "dat_stnd_cal.csv"), row.names=F)
dat.c <- read.csv(paste0(wd.out, mod.name, "dat_stnd_cal.csv")) ## has been updated with scrub outliers removed

dat.v <- dat_stnd[!(dat_stnd$id %in% dat.c$id),] 
#write.csv(dat.v, paste0(wd.out, mod.name,  "dat_stnd_val.csv"), row.names=F)
#dat.v <- read.csv(paste0(wd.out, mod.name, "dat_stnd_val.csv")) ## has been updated with scrub outliers removed



#### Set up the cluster ####
# Number of cores limited by computing resources - best on server
clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
clust <- try(makeCluster(getOption("cl.cores", 30), type = clusterType))
clusterExport(clust, "dat.c") ## Send the data to each node of the cluster
clusterCall(clust, function() library(MuMIn)) ## Call the function to load the library on each node of the cluster
# stopCluster(clust)


#### PART 1 ####
#### Model with linear terms only ####
vars.m <- paste(vars, collapse=" + ")
# Make model with calibration data.
ml <- glm(as.formula(paste0("occ ~ ", vars.m)), family=binomial, dat.c, na.action=na.fail)
## A model to be dredged must not have na.action as na.omit or na.exclude, or dredge will fail.

## Check for variance inflation - an effect of collinear explanatory variables
car::vif(ml, singular.ok = TRUE) ## combo_broadl_m (4.32) and treehedge (5.94) are borderline high (6-7)

# Explanatory deviance
dev_e <- 1-(ml$deviance/ml$null.deviance) # only 0.129

## Dredge model with linear terms. 
ml_drg <- dredge(ml, beta="sd", evaluate=T, trace=2, cluster=clust) # remove cluster when not using
save(ml_drg, file=paste0(wd.out, mod.name, "linear_dredged_all"))
load(paste0(wd.out, mod.name, "linear_dredged_all"))

# Extract model runs with dAIC<2 -- indistinguishable from each other re: model performance
ml_drg_df <- ml_drg[ml_drg$delta<=2,]
ml_drg_df <- as.data.frame(ml_drg_df)
write.csv(ml_drg_df, paste0(wd.out, mod.name, "linear_dredged_best_df"))
# 35 best models
# rarely included -> "combo_conif_c", "OS_Terrain_100_slope_pct", "tasrng_win_5yrmn"

# # Save the best model runs
ml_best <- get.models(ml_drg, subset=delta<2)
save(ml_best, file=paste0(wd.out, mod.name, "linear_dredged_best"))



#### PART 2 ####
load(file=paste0(wd.out, mod.name, "linear_dredged_best"))

#### Model with linear and quadratic terms ####
# Get the linear terms that were selected in the top models and make into quadratic terms
vars.mld <- unique(unlist(lapply(ml_best, function(x) {rownames(summary(x)$coefficients)})))
vars.mld <- vars.mld[vars.mld!="(Intercept)"] # Kept all 15 variables
# remove combo_broadl_m as it's correlated with treehedge --> more interested in the latter
vars.mld <- vars.mld[vars.mld != "combo_broadl_m"]

# Remove rarely included variables - don't try as quadratic.
lin_only <- c("combo_conif_c", "OS_Terrain_100_slope_pct", "tasrng_win_5yrmn")
i <- 1
for(i in 1:length(lin_only)) {
  vars.mld <- vars.mld[vars.mld!=lin_only[i]]
}
# Quadratic terms
vars.mq <- paste0("I(", vars.mld, "^2)")
# all linear and quadratic terms.
vars <- vars[vars != "combo_broadl_m"] # take out broadleaf mean
vars.mlq <- paste(c(vars, vars.mq), collapse=" + ")

# Dredge the model with all linear terms retained in the best model subset, and all of their quadratic terms
mlq <- glm(as.formula(paste0("occ ~ ", vars.mlq)), family=binomial, dat.c, na.action=na.fail)
# check variance inflation
# car::vif(mlq)
# OS_Terrain_100_NS = 8.821634e+04
# OS_Terrain_100_WE = 9.066272e+03
# I(combo_broadl_m^2) = 8.905545e+00
# I(OS_Terrain_100_NS^2) = 4.296280e+07
# I(OS_Terrain_100_WE^2) = 4.280311e+07
# I(treehedge^2) = 4.132180e+00

# Terrains may be because of inflation with their own quadratics
# Treehedge and broadl are correlated  ---> broadl taken out, now lower 

# Retain quadratic term only if linear term is added
msubset <- expression(dc("anc_wood", "I(anc_wood^2)") &
                        dc("combo_broadl_c", "I(combo_broadl_c^2)" ) &
                        dc("combo_conif_m", "I(combo_conif_m^2)") & ##
                        dc("OS_Terrain_100_NS", "I(OS_Terrain_100_NS^2)") &
                        dc("OS_Terrain_100_WE", "I(OS_Terrain_100_WE^2)") &
                        dc("rainfall_spr_5yrmn", "I(rainfall_spr_5yrmn^2)") &
                        dc("scrub", "I(scrub^2)") &
                        dc("sun_spr_5yrmn", "I(sun_spr_5yrmn^2)") & ##
                        dc("tasmax_spr_5yrmn", "I(tasmax_spr_5yrmn^2)") & ##
                        dc("tasmin_win_5yrmn", "I(tasmin_win_5yrmn^2)") &
                        dc("treehedge", "I(treehedge^2)"))


# Dredge subset of variables
mlqd <- dredge(mlq, beta="sd", evaluate=T, trace=2, subset=msubset, cluster=clust)
#save(mlqd, file=paste0(wd.out, mod.name, "linear_quadratic_dredged_subset"))
load(paste0(wd.out, mod.name, "linear_quadratic_dredged_subset"))

# Isolate best models
mlqd_df <- mlqd[mlqd$delta<=2,]
View(mlqd_df)
mlqd_df <- as.data.frame(mlqd_df)
#write.csv(mlqd_df, paste0(wd.out, mod.name, "linear_quadratic_dredged_subset_best_df"))

mlqd_best <- get.models(mlqd, subset=delta<2) # 80 best models
#save(mlqd_best, file=paste0(wd.out, mod.name, "linear_quadratic_dredged_subset_best"))
load(paste0(wd.out, mod.name, "linear_quadratic_dredged_subset_best"))

# Summarize the number of "best" models that each variable was included in
df_name <- c("variable", "num_models_included")
df <- data.frame(matrix(nrow = 0, ncol = length(df_name)))
names(df) <- df_name

i <- 1
for(i in 2:ncol(mlqd_df)) { # i = 1 is intercept column; don't want
  n <- nrow(mlqd_df[!is.na(mlqd_df[[i]]),]) # if column has all NAs then get warning but just assigned n as zero so it's correct
  r <- cbind(names(mlqd_df)[i], n)
  df[c(i-1),] <- r
}
View(mlqd_df_summary)
mlqd_df_summary <- df
mlqd_df_summary <- mlqd_df_summary %>%
  mutate(num_models_included = as.numeric(num_models_included)) %>%
  mutate(prop = round(num_models_included/nrow(mlqd_df),3)) %>%
  dplyr::filter(!(variable %in% c("df", "logLik", "AICc", "delta", "weight", "(Intercept)")))

#write.csv(mlqd_df_summary, paste0(wd.out, mod.name, "linear_quadratic_dredged_best_summary.csv"), row.names = F)

unique(unlist(lapply(mlqd_best, function(x) {rownames(summary(x)$coefficients)}))) ## 21 variables.

## Check for variance inflation within the model runs
for(i in 1:length(mlqd_best)) {print(vif(mlqd_best[[i]]), singular.ok = TRUE)}
# Inspect models (ask whether fixed effect relationships meaningful)
for(i in 1:length(mlqd_best)) {print(summary(mlqd_best[[i]]))}



#### PART 3 ####

#### Model with interactions between linear terms ####
## Create list of interactions to test -> removes those not in any "best" models
load(paste0(wd.out, mod.name, "linear_quadratic_dredged_subset_best")) # mlqd_best
vars.mlqd <- unique(unlist(lapply(mlqd_best, function(x) {rownames(summary(x)$coefficients)})))
vars.mlqd <- vars.mlqd[vars.mlqd!="(Intercept)"]

# Make interactions with only the linear terms - just used those commonly in the models (at least 25%) 
lin_var <- vars.mlqd[!grepl("I", vars.mlqd)] 
exclude <- c("combo_broadl_c", "tasrng_win_5yrmn", "OS_Terrain_100_slope_pct", "anc_wood")
lin_var <- lin_var[!lin_var %in% exclude] # 9 variables remain

vars.mi <- combn(lin_var, 2) # 36 combinations
vars.mi <- apply(vars.mi, MARGIN=2, FUN=function(x) {paste(x, collapse="*")}) # Turns into interactions for formula


## Group interactions into subsets, otherwise there are too many fixed terms for model selection to work
# only add 4 interactions -> max 25 variables for sake of run speed
  ## Each subset was run separately as background job as they each take hours to days to run ##
vars.mi.1 <- paste(c(vars.mlqd, vars.mi[1:4]), collapse= " + ")
vars.mi.2 <- paste(c(vars.mlqd, vars.mi[5:8]), collapse= " + ")
vars.mi.3 <- paste(c(vars.mlqd, vars.mi[9:12]), collapse= " + ") 
vars.mi.4 <- paste(c(vars.mlqd, vars.mi[13:16]), collapse= " + ") 
vars.mi.5 <- paste(c(vars.mlqd, vars.mi[17:20]), collapse= " + ") 
vars.mi.6 <- paste(c(vars.mlqd, vars.mi[21:24]), collapse= " + ") 
vars.mi.7 <- paste(c(vars.mlqd, vars.mi[25:28]), collapse= " + ") 
vars.mi.8 <- paste(c(vars.mlqd, vars.mi[29:32]), collapse= " + ")
vars.mi.9 <- paste(c(vars.mlqd, vars.mi[33:36]), collapse= " + ") 

mlqi.full1 <- glm(as.formula(paste0("occ ~ ", vars.mi.1)), family=binomial, dat.c, na.action=na.fail) 
mlqi.full2 <- glm(as.formula(paste0("occ ~ ", vars.mi.2)), family=binomial, dat.c, na.action=na.fail) 
mlqi.full3 <- glm(as.formula(paste0("occ ~ ", vars.mi.3)), family=binomial, dat.c, na.action=na.fail)
mlqi.full4 <- glm(as.formula(paste0("occ ~ ", vars.mi.4)), family=binomial, dat.c, na.action=na.fail)
mlqi.full5 <- glm(as.formula(paste0("occ ~ ", vars.mi.5)), family=binomial, dat.c, na.action=na.fail)
mlqi.full6 <- glm(as.formula(paste0("occ ~ ", vars.mi.6)), family=binomial, dat.c, na.action=na.fail)
mlqi.full7 <- glm(as.formula(paste0("occ ~ ", vars.mi.7)), family=binomial, dat.c, na.action=na.fail) 
mlqi.full8 <- glm(as.formula(paste0("occ ~ ", vars.mi.8)), family=binomial, dat.c, na.action=na.fail) 
mlqi.full9 <- glm(as.formula(paste0("occ ~ ", vars.mi.9)), family=binomial, dat.c, na.action=na.fail)


# Retain quadratic term only if linear term is added
msubset <- expression(dc("combo_conif_m", "I(combo_conif_m^2)") & 
                        dc("OS_Terrain_100_NS", "I(OS_Terrain_100_NS^2)") &
                        dc("OS_Terrain_100_WE", "I(OS_Terrain_100_WE^2)") &
                        dc("rainfall_spr_5yrmn", "I(rainfall_spr_5yrmn^2)") &
                        dc("scrub", "I(scrub^2)") & 
                        dc("tasmax_spr_5yrmn", "I(tasmax_spr_5yrmn^2)") & 
                        dc("tasmin_win_5yrmn", "I(tasmin_win_5yrmn^2)") &
                        dc("treehedge", "I(treehedge^2)"))



# Group 1 -> Run in background script "sdm_update_interaction_grp1.R" --> FINISHED
mlqid.1 <- dredge(mlqi.full1, beta="sd", evaluate=T, trace=2, subset = msubset, cluster=clust)
#save(mlqid.1, file=paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp1")) #### name changed so it doesn't get overwritten.
load(paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp1"))

mlqid.1_df <- as.data.frame(mlqid.1[mlqid.1$delta<=2,])
mlqid.1_best <- get.models(mlqid.1, subset=delta<2)
#write.csv(mlqid.1_df , paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp1_best_df.csv"), row.names = F)
#save(mlqid.1_best, file=paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp1_best"))
mlqid.1_df <- read.csv(paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp1_best_df.csv"))


# Group 2 -> Run in background script "sdm_update_interaction_grp2.R" --> FINISHED
mlqid.2 <- dredge(mlqi.full2, beta="sd", evaluate=T, trace=2, subset = msubset, cluster=clust)
#save(mlqid.2, file=paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp2")) 
load(paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp2"))

mlqid.2_df <- as.data.frame(mlqid.2[mlqid.2$delta<=2,])
mlqid.2_best <- get.models(mlqid.2, subset=delta<2)
#write.csv(mlqid.2_df , paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp2_best_df.csv"), row.names = F)
#save(mlqid.2_best, file=paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp2_best"))
mlqid.2_df <- read.csv(paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp2_best_df.csv"))


# Group 3 -> Run in background script "sdm_update_interaction_grp3.R" -> FINISHED
mlqid.3 <- dredge(mlqi.full3, beta="sd", evaluate=T, trace=2, subset = msubset, cluster=clust)
#save(mlqid.3, file=paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp3")) 
load(paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp3"))

mlqid.3_df <- as.data.frame(mlqid.3[mlqid.3$delta<=2,])
mlqid.3_best <- get.models(mlqid.3, subset=delta<2)
#write.csv(mlqid.3_df, paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp3_best_df.csv"), row.names = F)
#save(mlqid.3_best, file=paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp3_best"))
mlqid.3_df <- read.csv(paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp3_best_df.csv"))


# Group 4 -> Run in background script "sdm_update_interaction_grp4.R" -> FINISHED
mlqid.4 <- dredge(mlqi.full4, beta="sd", evaluate=T, trace=2, subset = msubset, cluster=clust)
#save(mlqid.4, file=paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp4")) 
load(paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp4"))

mlqid.4_df <- as.data.frame(mlqid.4[mlqid.4$delta<=2,])
mlqid.4_best <- get.models(mlqid.4, subset=delta<2)
#write.csv(mlqid.4_df, paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp4_best_df.csv"), row.names = F)
#save(mlqid.4_best, file=paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp4_best"))
mlqid.4_df <- read.csv(paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp4_best_df.csv"))


# Group 5 -> Run in background script "sdm_update_interaction_grp5.R" -> FINISHED
mlqid.5 <- dredge(mlqi.full5, beta="sd", evaluate=T, trace=2, subset = msubset, cluster=clust)
#save(mlqid.5, file=paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp5"))
load(paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp5"))

mlqid.5_df <- as.data.frame(mlqid.5[mlqid.5$delta<=2,])
mlqid.5_best <- get.models(mlqid.5, subset=delta<2)
#write.csv(mlqid.5_df, paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp5_best_df.csv"), row.names = F)
#save(mlqid.5_best, file=paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp5_best"))
mlqid.5_df <- read.csv(paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp5_best_df.csv"))


# Group 6 -> Run in background script "sdm_update_interaction_grp6.R" -> FINISHED
mlqid.6 <- dredge(mlqi.full6, beta="sd", evaluate=T, trace=2, subset = msubset, cluster=clust)
#save(mlqid.6, file=paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp6")) 
load(paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp6"))

mlqid.6_df <- as.data.frame(mlqid.6[mlqid.6$delta<=2,])
mlqid.6_best <- get.models(mlqid.6, subset=delta<2)
#write.csv(mlqid.6_df, paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp6_best_df.csv"), row.names = F)
#save(mlqid.6_best, file=paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp6_best"))
mlqid.6_df <- read.csv(paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp6_best_df.csv"))


# Group 7 -> Run in background script "sdm_update_interaction_grp7.R" -> FINISHED
mlqid.7 <- dredge(mlqi.full7, beta="sd", evaluate=T, trace=2, subset = msubset, cluster=clust)
#save(mlqid.7, file=paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp7")) 
load(paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp7"))

mlqid.7_df <- as.data.frame(mlqid.7[mlqid.7$delta<=2,])
mlqid.7_best <- get.models(mlqid.7, subset=delta<2)
#write.csv(mlqid.7_df, paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp7_best_df.csv"), row.names = F)
#save(mlqid.7_best, file=paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp7_best"))
mlqid.7_df <- read.csv(paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp7_best_df.csv"))


# Group 8 -> Run in background script "sdm_update_interaction_grp8.R" -> FINISHED
mlqid.8 <- dredge(mlqi.full8, beta="sd", evaluate=T, trace=2, subset = msubset, cluster=clust)
#save(mlqid.8, file=paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp8")) 
load(paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp8"))

mlqid.8_df <- as.data.frame(mlqid.8[mlqid.8$delta<=2,])
mlqid.8_best <- get.models(mlqid.8, subset=delta<2)
#write.csv(mlqid.8_df, paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp8_best_df.csv"), row.names = F)
#save(mlqid.8_best, file=paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp8_best"))
mlqid.8_df <- read.csv(paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp8_best_df.csv"))


# Group 9 -> Run in background script "sdm_update_interaction_grp9.R -> FINISHED
mlqid.9 <- dredge(mlqi.full9, beta="sd", evaluate=T, trace=2, subset = msubset, cluster=clust)
#save(mlqid.9, file=paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp9")) 
load(paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp9"))

mlqid.9_df <- as.data.frame(mlqid.9[mlqid.9$delta<=2,])
mlqid.9_best <- get.models(mlqid.9, subset=delta<2)
#write.csv(mlqid.9_df, paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp9_best_df.csv"), row.names = F)
#save(mlqid.9_best, file=paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp9_best"))
mlqid.9_df <- read.csv(paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_grp9_best_df.csv"))



# Summarize the number of "best" models that each variable was included in
# Combine all data into a list and remove unnecessary columns
grp_df <- list(mlqid.1_df, mlqid.2_df, mlqid.3_df, mlqid.4_df, mlqid.5_df, mlqid.6_df, mlqid.7_df, mlqid.8_df, mlqid.9_df)
grp_df <- lapply(grp_df, function(x) subset(x, select = -c(X.Intercept., df, logLik, AICc, delta, weight)))
# Intercept named this way (with X.) if read in from csv -- adjust as needed

empty_name <- c("subset", "variable", "num_models_included") # Template for output dataframe to save out
empty <- data.frame(matrix(nrow = 0, ncol = length(empty_name)))
names(empty) <- empty_name

# to be populated with the output dataframes
summary_list <- list()

i <- 1
j <- 1
for(j in 1:length(grp_df)) { # for each item in list
  tmp <- as.data.frame(grp_df[j]) # isolate and convert to dataframe
  grp <- paste0("subset_", j) # create a name for this subset
  df <- empty # copy the template so it can begin empty each time
  for(i in 1:ncol(tmp)) { 
    n <- nrow(tmp[!is.na(tmp[[i]]),]) # if column has all NAs then get warning but just assigned n as zero so it's correct
    r <- cbind(grp, names(tmp)[i], n)
    df[i,] <- r 
  }
  df <- df %>%
    mutate(num_models_included = as.numeric(num_models_included)) %>%
    mutate(prop = round(num_models_included/nrow(tmp),3))
  
  summary_list[[j]] <- df
  write.csv(df, paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_summary_", grp, ".csv"), row.names = F)
}


# Combine all groups into single object and subset to best models.
mlqid <- merge(mlqid.1, mlqid.2)
mlqid <- merge(mlqid, mlqid.3)
mlqid <- merge(mlqid, mlqid.4)
mlqid <- merge(mlqid, mlqid.5)
mlqid <- merge(mlqid, mlqid.6)
mlqid <- merge(mlqid, mlqid.7)
mlqid <- merge(mlqid, mlqid.8)
mlqid <- merge(mlqid, mlqid.9)
mlqid_df <- as.data.frame(mlqid[mlqid$delta<=2,])
mlqid_best <- get.models(mlqid, subset=delta<2)
View(mlqid_df)
write.csv(mlqid_df, paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_combine_best_df.csv"), row.names = F)
save(mlqid_best, file=paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_combine_best"))

# Summarize the number of "best" models that each variable was included in across all subsets
temp <- subset(mlqid_df, select = -c(df, logLik, AICc, delta, weight))

df_name <- c("variable", "num_models_included")
df <- data.frame(matrix(nrow = 0, ncol = length(df_name)))
names(df) <- df_name

i <- 1
for(i in 2:ncol(temp)) { # i = 1 is intercept column; don't want
  n <- nrow(temp[!is.na(temp[[i]]),]) # if column has all NAs then get warning but just assigned n as zero so it's correct
  r <- cbind(names(temp)[i], n)
  df[c(i-1),] <- r 
}

mlqid_df_summary <- df
mlqid_df_summary <- mlqid_df_summary %>%
  mutate(num_models_included = as.numeric(num_models_included)) %>%
  mutate(prop = round(num_models_included/nrow(mlqid_df),3))

write.csv(mlqid_df_summary, paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_best_combine_summary.csv"), row.names = F) 
mlqid_df_summary <- read.csv(paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_best_combine_summary.csv"))
# View(mlqid_df_summary)


#### PART 4 #### 
# Reduce variables to those deemed best based on subsets 
# Pull variable names from the mlqid_best which summarized the best from all subsets, but remove and add vars as indicated in notes 

# Variables in best models when dredges combined 
vars.mlqid.all <- mlqid_df_summary$variable[mlqid_df_summary$num_models_included !=0] 
# Manually remove those deemed to have low importance
exclude <- c("rainfall_spr_5yrmn:tasmin_win_5yrmn", "scrub:tasmax_spr_5yrmn")
# remove variables from variable list
vars.mlqid.all <- vars.mlqid.all[!(vars.mlqid.all %in% exclude)] # 21 variables remain

# Interactions that were in the majority of the local best model subset 
extra.inter <- c("combo_conif_m:rainfall_spr_5yrmn", "combo_conif_m:sun_spr_5yrmn", 
                 "OS_Terrain_100_NS:rainfall_spr_5yrmn", "OS_Terrain_100_WE:rainfall_spr_5yrmn", 
                 "sun_spr_5yrmn:tasmax_spr_5yrmn", "sun_spr_5yrmn:treehedge",
                 "tasmax_spr_5yrmn:treehedge")

# Combine the variables from best model when dredges combined and the extra interactions of interest 
vars.max <- paste(c(vars.mlqid.all, extra.inter), collapse = " + ") # 28 variables


#### Run model with this variable subset ####
mlqi.comb <- glm(as.formula(paste0("occ ~ ", vars.max)), family=binomial, dat.c, na.action=na.fail) 
# vif(mlqi.comb) # conifer mean is high with it's quadratic - not an issue

# Retain quadratic term only if linear term is added
msubset <- expression(dc("combo_conif_m", "I(combo_conif_m^2)") &  
                        dc("OS_Terrain_100_WE", "I(OS_Terrain_100_WE^2)") &
                        dc("rainfall_spr_5yrmn", "I(rainfall_spr_5yrmn^2)") &
                        dc("scrub", "I(scrub^2)") & 
                        dc("tasmax_spr_5yrmn", "I(tasmax_spr_5yrmn^2)") & 
                        dc("tasmin_win_5yrmn", "I(tasmin_win_5yrmn^2)") &
                        dc("treehedge", "I(treehedge^2)"))


# Dredge the max variable model
mlqid.comb <- dredge(mlqi.comb, beta="sd", evaluate=T, trace=2, subset=msubset, cluster=clust) 
save(mlqid.comb, file=paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_best_vars"))

# Isolate best models
mlqid.comb_df <- mlqid.comb[mlqid.comb$delta<=2,]  
mlqid.comb_df <- as.data.frame(mlqid.comb_df)
# write.csv(mlqid.comb_df, paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_best_vars_best_mods.csv"), row.names = F)
mlqid.comb_df <- read.csv(paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_best_vars_best_mods.csv"))
View(mlqid.comb_df)

mlqid.comb_best <- get.models(mlqid.comb, subset=delta<2)
save(mlqid.comb_best, file=paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_best_vars_best_mods"))

# Summarize the number of "best" models that each variable was included in across all subsets
temp <- subset(mlqid.comb_df, select = -c(df, logLik, AICc, delta, weight))

df_name <- c("variable", "num_models_included")
df <- data.frame(matrix(nrow = 0, ncol = length(df_name)))
names(df) <- df_name

i <- 1
for(i in 2:ncol(temp)) { # i = 1 is intercept column; don't want
  n <- nrow(temp[!is.na(temp[[i]]),]) # if column has all NAs then get warning but just assigned n as zero so it's correct
  r <- cbind(names(temp)[i], n)
  df[c(i-1),] <- r 
}

mlqid.comb_df_summary <- df
mlqid.comb_df_summary <- mlqid.comb_df_summary %>%
  mutate(num_models_included = as.numeric(num_models_included)) %>%
  mutate(prop = round(num_models_included/nrow(mlqid.comb_df),3))

write.csv(mlqid.comb_df_summary, paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_best_vars_best_mods_summary.csv"), row.names = F) 



#### PART 5 ####
#### Produce & save final model ####
# Load the final dredged model output
load(paste0(wd.out, mod.name, "linear_quadratic_interactions_dredged_best_vars")) # mlqid.comb

final_m01 <- model.avg(mlqid.comb, subset = delta<=2, fit=T) ## Average the models with delta AIC < value and ensure the models are actually fit
s <- summary(final_m01)

# Output the coefficients in a table
final_m01_coefs <- final_m01$coefficients["full",] ## Average parameter estimates over all models (value will be 0 in models wehre the parameter wasn't included)
vars <- names(final_m01_coefs)
coefs <- as.vector(final_m01_coefs)
final_m01_coefs <- data.frame(Variables=vars, MeanCoefficients=coefs)
rownames(final_m01_coefs) <- vars

# Output summary of model with standard error
final_m01_summary_full <- as.data.frame(s$coefmat.full)
View(final_m01_summary_full)
final_m01_summary_full$variable <- rownames(final_m01_summary_full)
final_m01_summary_full <- final_m01_summary_full %>%
  relocate(variable, .before = Estimate) %>%
  rename(mean_coefficient = Estimate)

# Write outputs
save(final_m01, file=paste0(wd.out, mod.name,"final_averaged_thawfix"))
write.csv(final_m01_coefs, file=paste0(wd.out, mod.name,"final_averaged_thawfix.csv"), row.names=F)
write.csv(final_m01_summary_full, paste0(wd.out, mod.name, "final_averaged_ModSummary_thawfix.csv"), row.names = F)

View(final_m01_coefs)
