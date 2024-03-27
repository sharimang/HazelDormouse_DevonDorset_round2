#### Data Exploration ####

## Visualizing relationships between occurrence and environmental variables ##
## Generating plots and data summaries ##

## Written by: Shari Mang ##
## Date: November 2023 ##
## Based on code by: Regan Early, written in 2021 ##

# Plotting info - might be helpful 
# http://www.sthda.com/english/wiki/be-awesome-in-ggplot2-a-practical-guide-to-be-highly-effective-r-software-and-data-visualization


#### Set up ####

# pacman::p_install()
pacman::p_load(
  sf,
  tidyverse, 
  here,
  ggplot2,
  tidyr,
  dplyr,
  scales,
  ggpubr,
  MuMIn, # stdize()
  corrplot,
  Hmisc,  # rcorr()
  conflicted
)
conflicted::conflict_prefer("here", "here")

wd.fig <- here("model/explore_data/")

dat <- read.csv(here("data/env_occ_clean.csv")) 


### Plotting ####
vars <- c('combo_broadl_m', 'combo_broadl_c', 'combo_conif_m', 'combo_conif_c', 'anc_wood', 
          'OS_Terrain_100_NS', 'OS_Terrain_100_WE', 'OS_Terrain_100_slope_pct', 'treehedge', 'scrub', 
          "tasrng_win_5yrmn", "tasmax_spr_5yrmn", "tasmin_win_5yrmn", "sun_spr_5yrmn", "rainfall_spr_5yrmn")

# > Histograms for enviro variables 
dat.1 <- dat[,vars[1:8]]
dat.2 <- dat[,vars[9:15]]

dat.1 %>%
  gather() %>%
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram() + 
  theme_bw()
ggsave(plot = last_plot(), file = paste0(wd.fig, "hist_variables_1.png"),  
       width = 30, height = 25, units = "cm")

dat.2 %>%
  gather() %>%
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram()+ 
  theme_bw()
ggsave(plot = last_plot(), file = paste0(wd.fig, "hist_variables_2.png"),  
       width = 30, height = 25, units = "cm")


#### Observations across environmental variables ####
# > Distribution of presence and absence across environmental variables 

dat.hist <- dat[,c("occ",vars)]
dat.hist <- dat.hist %>%
  dplyr::mutate(pres_abs = if_else(dat.hist$occ == 1, "Present", "Absent")) %>%
  dplyr::select(-occ)
dat.hist <- dat.hist %>%
  gather(key = "env_var", value = "value", - pres_abs)

i <- 1
for(i in 1:length(vars)) {
  p <- ggplot(subset(dat.hist, env_var %in% vars[i]), 
         aes(x = value, color = pres_abs)) + 
    geom_histogram(position="identity", fill = "white") + 
    xlab(vars[i])  +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 18), 
          legend.text = element_text(size = 14),
          legend.title = element_text(size=14),
          plot.title = element_text(size = 18), 
          strip.text = element_text(size = 14)) +
    facet_wrap(~ pres_abs) 
  ggsave(p, file = paste0(wd.fig, "hist_occ_", vars[i], ".png"),  
         width = 30, height = 21, units = "cm")
}


# # If I want the histograms to overlap use this
# ggplot(dat, aes(x = combo_broadl_m, color = occ)) + 
#   geom_histogram(data=subset(dat, occ == '0'),color = "green", fill = "white", alpha = 0.5) + 
#   geom_histogram(data=subset(dat, occ == '1'),color = "blue", fill = "white", alpha = 0.5)



# > Variation within enviro variables across data sources ####
dat.plot <- dat[, c("source",vars)]
dat.plot <- dat.plot %>%
  gather(-source, key = "var", value = "value") 

ggplot(dat.plot, aes(x=as.factor(source), y=value)) +
  geom_boxplot() +
  facet_wrap(~ var, scales = "free") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) 
ggsave(last_plot(), file = paste0(wd.fig, "vars_boxplot_source.png"),  
       width = 25, height = 30, units = "cm")

# Histograms of environmental variables split by data source
dat.plot.h <- dat.plot[dat.plot$source!= "pseudo-absence", ]

i <- 1
for(i in 1:length(vars)) {
  p <- ggplot(subset(dat.plot.h, var %in% vars[i]), 
              aes(x = value, color = source)) + 
    geom_histogram(position="identity", fill = "white") + 
    xlab(vars[i])  +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 18), 
          legend.text = element_text(size = 14),
          legend.title = element_text(size=14),
          plot.title = element_text(size = 18), 
          strip.text = element_text(size = 14)) +
    facet_wrap(~ source) 
  ggsave(p, file = paste0(wd.fig, "hist_source_", vars[i], ".png"),  
         width = 30, height = 21, units = "cm")
}


#### Correlation ####
corrs1 <- cor(dat[,vars], method = "spearman") #, use="complete.obs"
# plot correlation
jpeg(paste0(wd.fig, "corrplot_spearman.jpg"), width = 700, height = 700)
#corrplot::corrplot(corrs1, method = "number")
corrplot::corrplot.mixed(corrs1, tl.pos='lt', tl.cex=0.9, number.cex=0.9, addCoefasPercent=T)
dev.off()
cordf <- as.data.frame(round(corrs1, digits = 3))

# > Extract correlations > 0.7 ####
corrs_sp <- Hmisc::rcorr(as.matrix(dat[,vars]), type = "spearman") 
corrs_sp$r[lower.tri(corrs_sp$r)] <- NA # stops duplicates
a <- as.data.frame(which(abs(corrs_sp$r) >= 0.7 & abs(corrs_sp$r) != 1, arr.ind = TRUE))
r <- rownames(corrs_sp$r)[a$row]
c <- colnames(corrs_sp$r)[a$col]

corrs_sp_out <- data.frame(matrix(nrow=nrow(a), ncol=2, dimnames=list(c(), c("vars", "corrs"))), stringsAsFactors=F)
corrs_sp_out$vars <- paste0(r, " & ", c)
corrs_sp_out$corrs <- apply(a, MARGIN=1, FUN=function(x) { corrs_sp$r[x[1], x[2]] })
corrs_sp_out$method <- "spearman"






