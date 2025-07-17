# Fit models of verticality using sdmTMB (SECOND ROUND OF MODELS USING FHD AND REVISED CANOPY HEIGHT)
# Step 1: compare random effects structure using cross validation, reml = T, include all fixed effects
# Step 2: compare fixed effects structure for different hypotheses regarding environmental predictors using cross validation, reml = F
# Step 3: refit the best model using all data
# Step 4: predict the best model to future env data

library(sdmTMB)
library(lme4)
library(tidyverse)
library(terra)
library(car)
library(data.table)
library(DHARMa)
source("scripts/00_functions/manuscript_functions.R")
source("scripts/00_functions/00_plot_functions.R")

# mammals forest only 50km resolution
raw = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/mammals_comdat.csv")
v = vect("data/original/rnaturalearth_world.shp")
v = project(v, "+proj=cea +datum=WGS84")


# remove cells with fewer than five species
raw = raw %>% 
  # subset to richness >= 5
  filter(rich >= 5) %>% 
  filter(!is.na(vert.mean.ses)) %>% 
  drop_na(!fhd) %>% 
  dplyr::select(!canopy_height)

# plot relationship between SES verticality and env predictors
raw %>% 
  dplyr::select(vert.mean.ses, biome:tmin_cold, elev, veg_den, veg_complexity,
                canopy_height2, fhd, precip_warm) %>%
  pivot_longer(cols = 3:17, names_to = "var", values_to = "val") %>% 
  ggplot(aes(x = val, y = vert.mean.ses, color = biome)) +
  geom_point(pch = ".") +
  facet_wrap(~var, scales = "free") +
  theme_classic()

dat.t = raw %>% 
  mutate(log_clim_velocity = log10(clim_velocity),
         log_tmin_cold = log10(tmin_cold + 50),
         log_precip_dry = log10(precip_dry + 1),
         log_precip_wet = log10(precip_wet))

dat.t %>% 
  dplyr::select(vert.mean.ses, biome:tmin_cold, elev, veg_den, veg_complexity,
                canopy_height2, fhd, precip_warm, log_clim_velocity, log_precip_dry, 
                log_tmin_cold, log_precip_wet) %>%
  pivot_longer(cols = 3:21, names_to = "var", values_to = "val") %>% 
  ggplot(aes(x = val, y = vert.mean.ses, color = biome)) +
  geom_point(pch = ".") +
  facet_wrap(~var, scales = "free") +
  theme_classic()

# scale the data
dat = dat.t %>%
  mutate_at(.vars = vars(mat:tmin_cold, elev, veg_den, veg_complexity, 
                         canopy_height2, fhd, log_precip_dry, precip_warm), .funs = scale)


# test correlation between predictor variables
dat %>% 
  dplyr::select(canopy_height2, veg_den, fhd, veg_complexity,
                tmax_warm, tmin_cold, temp_sea,
                precip_wet, log_precip_dry, precip_sea, precip_warm) %>% 
  slice_sample(n = 1000) %>% 
  pairs()

cortestdat = dat %>% 
  dplyr::select(canopy_height2, veg_den, fhd, veg_complexity,
                tmax_warm, tmin_cold, temp_sea, temp_diu,
                precip_wet, log_precip_dry, precip_sea, precip_warm)

cortest = cor(cortestdat)
cor.test(dat$fhd, dat$canopy_height2)

# fhd is very strongly correlated with canopy height
# tmin_cold is very strongly correlated with temp seasonality
# precip_dry is very strongly correlated with precip seasonality
# veg_complexity strongly correlated with canopy height, veg den, and fhd


vif(lm(vert.mean.ses ~ canopy_height2 + veg_den + tmax_warm + I(tmax_warm^2)  + temp_diu +
         tmin_cold + precip_warm + log_precip_dry, data = dat))

# VIF all under 5

dat$x = dat$x/1e5
dat$y = dat$y/1e5
dat$ecoregion = factor(dat$ecoregion)
dat$biome = factor(dat$biome)
dat$biorealm = factor(paste(dat$biome, dat$realm, sep = "_"))
dat$realm = factor(dat$realm)

# SES VERT ----------------------------------------------------------------


# DIURNAL TEMPERATURE RANGE MODELS ----------------------------------------

f1 = formula(vert.mean.ses ~ tmax_warm + I(tmax_warm^2) + tmin_cold + temp_diu +
               precip_warm + log_precip_dry + 
               canopy_height2 + veg_den +
               precip_warm:canopy_height2 + tmin_cold:canopy_height2 + 
               tmax_warm:canopy_height2 + I(tmax_warm^2):canopy_height2 +
               log_precip_dry:canopy_height2 + temp_diu:canopy_height2)

f2 = formula(vert.mean.ses ~ tmax_warm + I(tmax_warm^2) + tmin_cold + temp_diu +
               precip_warm + log_precip_dry + 
               canopy_height2 + veg_den)

f3 = formula(vert.mean.ses ~ 1)

forms = list(f1, f2, f3)


samp = dat %>% sample_n(1000)
samp.cor = ncf::correlog(x = samp$x, y = samp$y, z = samp$vert.mean.ses, increment = 50000/1e5, resamp = 99)
ncf:::plot.correlog(samp.cor)
ncf:::plot.correlog(samp.cor, xlim = c(0,100))

# set initial range as 30 and max.edge as range/5
# https://haakonbakkagit.github.io/btopic104.html
# Bakka, H., J. Vanhatalo, J. Illian, D. Simpson, and H. Rue. 2016. “Accounting for Physical Barriers in Species Distribution Modeling with Non-Stationary Spatial Random Effects.” arXiv preprint arXiv:1608.03787. Norwegian University of Science; Technology, Trondheim, Norway. 

# removed SVC for canopy height
fitmesh = fit_mesh(f1, dat, range = samp.cor$x.intercept, v = v, family = gaussian())

mesh = fitmesh$meshes[[length(fitmesh$meshes)]]
saveRDS(list(mesh = mesh), "results/sdmTMB_models2.3/mammals_sesvert_tempdiu.rds")
out = readRDS("results/sdmTMB_models2.3/mammals_sesvert_tempdiu.rds")
mesh = out$mesh

taxon = "Mammals"
response_var = "SES verticality"
fname_end = "mammals_sesvert"

mods.spatial = compareMods_AIC(f = forms, dat, mesh, taxon = taxon,
                               response_var = response_var, family = gaussian(), reml = F, spatial = "on")
mods.nonspatial = compareMods_AIC(f = forms, dat, mesh, taxon = taxon,
                                  response_var = response_var, family = gaussian(), reml = F, spatial = "off")

saveRDS(list(mesh = mesh, mods.spatial = mods.spatial, mods.nonspatial = mods.nonspatial), "results/sdmTMB_models2.3/mammals_sesvert_tempdiu.rds")



plot_resids(mod = fullmod, response_var = "vert.mean.ses", 
            fpath = paste0("figures/residual_checks/sdmTMB2.3/mammals_tempdiu/",fname_end), integer_response = F)


