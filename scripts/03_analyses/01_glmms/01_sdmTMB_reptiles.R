# Fit models of verticality using sdmTMB
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
library(fmesher)
source("scripts/00_functions/00_plot_functions.R")
source("scripts/00_functions/manuscript_functions.R")


# amphibian forest only 50km resolution
raw = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/reptiles_comdat.csv")

# remove cells with fewer than five species
raw = raw %>% 
  # subset to richness >= 5
  filter(rich >= 5) %>% 
  filter(!is.na(vert.mean.ses)) %>% 
  drop_na()

# future env data
env.f = rast(list.files(path = "data/derivative_data/resampled_env_rasters_50km/chelsa_future/2071_2100/ensemble/",
                        pattern = ".tif", full.names = T))
env.f = env.f[[2:7]]

dat.f = raw %>% 
  dplyr::select(x,y,canopy_height, veg_den, clim_velocity, realm)

env.f = terra::extract(env.f, dat.f[,c('x', "y")], ID = F)

dat.f = cbind(dat.f, env.f)
dat.f = dat.f %>% rename(tmax_warm = max_temp_warm, tmin_cold = min_temp_cold)

# remove NA values from future and present data
raw = raw[which(!is.na(dat.f$precip_dry)),]
dat.f = dat.f[which(!is.na(dat.f$precip_dry)),]



# plot relationship between SES verticality and env predictors
# raw %>% 
#   dplyr::select(vert.mean.ses, biome:clim_velocity, elev, veg_den, veg_complexity) %>%
#   pivot_longer(cols = 3:16, names_to = "var", values_to = "val") %>% 
#   ggplot(aes(x = val, y = vert.mean.ses, color = biome)) +
#   geom_point(pch = ".") +
#   facet_wrap(~var, scales = "free") +
#   theme_classic()

dat.t = raw %>% 
  mutate(log_clim_velocity = log10(clim_velocity),
         log_tmin_cold = log10(tmin_cold + 50),
         log_precip_dry = log10(precip_dry),
         log_precip_wet = log10(precip_wet))

# dat.t %>% 
#   dplyr::select(vert.mean.ses, biome:clim_velocity, elev, veg_den, veg_complexity, log_clim_velocity, log_precip_dry, log_tmin_cold, log_precip_wet) %>%
#   pivot_longer(cols = 3:20, names_to = "var", values_to = "val") %>% 
#   ggplot(aes(x = val, y = vert.mean.ses, color = biome)) +
#   geom_point(pch = ".") +
#   facet_wrap(~var, scales = "free") +
#   theme_classic()

# scale the data
dat = dat.t %>%
  mutate_at(.vars = vars(canopy_height:tmin_cold, elev, veg_den, veg_complexity, log_clim_velocity, log_precip_dry, precip_warm), 
            .funs = function(x) scale(x)[,1]) %>% 
  dplyr::select(!precip_dry)


# scale future data
dat.f = dat.f %>% 
  # take log of clim_velocity and precip_dry
  mutate(log_clim_velocity = log10(clim_velocity),
         log_precip_dry = log10(precip_dry)) %>% 
  mutate(tmax_warm = (tmax_warm - mean(dat.t$tmax_warm))/sd(dat.t$tmax_warm),
         tmin_cold = (tmin_cold-mean(dat.t$tmin_cold))/sd(dat.t$tmin_cold),
         log_precip_dry = (log10(precip_dry) - mean(dat.t$precip_dry))/sd(dat.t$precip_dry),
         precip_wet = (precip_wet - mean(dat.t$precip_wet))/sd(dat.t$precip_wet),
         precip_sea = (precip_sea - mean(dat.t$precip_sea))/sd(dat.t$precip_sea),
         precip_warm = (precip_warm - mean(dat.t$precip_warm))/sd(dat.t$precip_warm),
         canopy_height = (canopy_height - mean(dat.t$canopy_height))/sd(dat.t$canopy_height),
         veg_den = (veg_den - mean(dat.t$veg_den))/sd(dat.t$veg_den),
         log_clim_velocity = (log10(clim_velocity) - mean(dat.t$clim_velocity))/sd(dat.t$clim_velocity))       

# test correlation between predictor variables
dat %>% 
  dplyr::select(canopy_height, veg_den, veg_complexity,
                tmax_warm, tmin_cold, temp_sea,
                precip_wet, log_precip_dry, precip_sea,
                log_clim_velocity) %>% 
  pairs()

plot((dat$log_precip_dry - dat.f$log_precip_dry), (dat$precip_sea - dat.f$precip_sea))

# veg_complexity strongly correlated with canopy height and veg den
# temp sea strongly correlated with tmin_cold
vif(lm(vert.mean.ses ~ canopy_height + veg_den + veg_complexity + tmax_warm + tmin_cold + temp_sea + precip_sea + precip_warm + log_precip_dry + log_clim_velocity, data = dat))

# remove temp_sea
vif(lm(vert.mean.ses ~ canopy_height + veg_den + veg_complexity + tmax_warm + tmin_cold + precip_sea + log_precip_dry + precip_warm + log_clim_velocity, data = dat))

# remove veg_complexity
vif(lm(vert.mean.ses ~ canopy_height + veg_den + tmax_warm + tmin_cold  + precip_sea + log_precip_dry + precip_warm + log_clim_velocity, data = dat))

# remove precip_sea
vif(lm(vert.mean.ses ~ canopy_height + veg_den + tmax_warm + tmin_cold + log_precip_dry + precip_warm + log_clim_velocity, data = dat))


# VIF all under 4

#f1 = formula(vert.mean.ses ~ canopy_height + veg_den + I(tmax_warm^2) + tmax_warm + tmin_cold + precip_wet + log_precip_dry + precip_warm + log_clim_velocity)
f1 = formula(vert.mean.ses ~ tmax_warm + I(tmax_warm^2) + tmin_cold +
               precip_warm + log_precip_dry +
               canopy_height + veg_den + 
               log_clim_velocity +
               precip_warm:canopy_height + tmin_cold:canopy_height + 
               tmax_warm:canopy_height + I(tmax_warm^2):canopy_height + log_precip_dry:canopy_height)

dat$x = dat$x/1e5
dat$y = dat$y/1e5
dat$ecoregion = factor(dat$ecoregion)
dat$biome = factor(dat$biome)
dat$biorealm = factor(paste(dat$biome, dat$realm, sep = "_"))
dat$realm = as.factor(dat$realm)

dat.f$x = dat.f$x/1e5
dat.f$y = dat.f$y/1e5
dat.f$realm = factor(dat.f$realm)



# SES VERT ----------------------------------------------------------------

# * - Set up mesh for analysis ------------------------------------------------

# set up spatial mesh

taxon = "Reptiles"
response_var = "SES verticality"
fname_end = "reptiles_sesvert"

# first estimate range of spatial autocorrelation
samp = dat %>% sample_n(1000)
samp.cor = ncf::correlog(x = samp$x, y = samp$y, z = samp$vert.mean.ses, increment = 50000/1e5, resamp = 99)
ncf:::plot.correlog(samp.cor)
ncf:::plot.correlog(samp.cor, xlim = c(0,100))

# set initial range as 30 and max.edge as range/5
# https://haakonbakkagit.github.io/btopic104.html
# Bakka, H., J. Vanhatalo, J. Illian, D. Simpson, and H. Rue. 2016. “Accounting for Physical Barriers in Species Distribution Modeling with Non-Stationary Spatial Random Effects.” arXiv preprint arXiv:1608.03787. Norwegian University of Science; Technology, Trondheim, Norway. 

fitmesh = fit_mesh(f1, dat, range = samp.cor$x.intercept, v = v, family = gaussian())
mesh = fitmesh$meshes[[length(fitmesh$meshes)]]
saveRDS(list(mesh = mesh), "results/sdmTMB_models2/reptiles_sesvert.rds")
# out = readRDS("results/sdmTMB_models2/mammals_sesvert.rds")
# mesh = out$mesh

# save(fitmesh, file = paste0("results/sdmTMB_models/model_selection/",fname_end,".RData"))
# load(paste0("results/sdmTMB_models/model_selection/",fname_end,".RData"))


# * - Compare models with AIC -------------------------------------------------


# Compare models using all points with AIC
# Compare models using all points with AIC
forms = list(f1,
             update(f1, . ~ .-precip_warm:canopy_height),
             update(f1, . ~ .-log_precip_dry:canopy_height),
             update(f1, . ~ .-tmin_cold:canopy_height),
             update(f1, . ~ .-tmax_warm:canopy_height - I(tmax_warm^2):canopy_height),
             update(f1, . ~ .-precip_warm:canopy_height - log_precip_dry:canopy_height),
             update(f1, . ~ .-precip_warm:canopy_height - tmin_cold:canopy_height),
             update(f1, . ~ .-precip_warm:canopy_height - tmax_warm:canopy_height - I(tmax_warm^2):canopy_height),
             update(f1, . ~ .-log_precip_dry:canopy_height - tmin_cold:canopy_height),
             update(f1, . ~ .-log_precip_dry:canopy_height - tmax_warm:canopy_height - I(tmax_warm^2):canopy_height),
             update(f1, . ~ .-tmin_cold:canopy_height - tmax_warm:canopy_height - I(tmax_warm^2):canopy_height),
             update(f1, . ~ .-precip_warm:canopy_height - log_precip_dry:canopy_height - tmin_cold:canopy_height),
             update(f1, . ~ .-precip_warm:canopy_height - log_precip_dry:canopy_height - tmax_warm:canopy_height - I(tmax_warm^2):canopy_height),
             update(f1, . ~ .-log_precip_dry:canopy_height - tmin_cold:canopy_height - tmax_warm:canopy_height - I(tmax_warm^2):canopy_height),
             update(f1, . ~ .-precip_warm:canopy_height - log_precip_dry:canopy_height - tmin_cold:canopy_height - tmax_warm:canopy_height - I(tmax_warm^2):canopy_height))


compMods_aic = compareMods_AIC(f = forms, dat, mesh, taxon = taxon, response_var = response_var, family = gaussian(), reml = F)
saveRDS(list(mesh = mesh, compMods_aic = compMods_aic),  "results/sdmTMB_models2/reptiles_sesvert.rds")

compMods_aic$modsel
sanity(compMods_aic$modlist[[3]])
summary(compMods_aic$modlist[[3]])

# refit model using REML
bestmod = compareMods_AIC(f = list(forms[[3]]), dat, mesh, taxon = taxon, response_var = response_var, family = gaussian(), reml = T)
bestmod = bestmod$modlist[[1]]
saveRDS(list(mesh = mesh, compMods_aic = compMods_aic, bestmod = bestmod), "results/sdmTMB_models2/reptiles_sesvert.rds")

# 
# compMods_aic = compareMods_AIC(f1, dat, mesh, taxon = taxon, response_var = response_var, family = gaussian())
# sanity(compMods_aic$mods[[1]])
# sanity(compMods_aic$mods[[2]])
# sanity(compMods_aic$mods[[3]])
# sanity(compMods_aic$mods[[4]])
# 
# compMods_aic[[2]]
# save(compMods_aic, fitmesh, file = paste0("results/sdmTMB_models/model_selection/",fname_end,".RData"))
# load(paste0("results/sdmTMB_models/model_selection/",fname_end,".RData"))
# 
# compMods_aic[[2]]

# * - Cross validation --------------------------------------------------------

# compare random effects structure using cross validation
# random effects include random intercept, spatial random field, and spatially varying coefficient
# spatial random field and spatially varying coefficient both depend on mesh
# preliminary assessment indicated that residuals displayed strong spatial autocorrelation when
# spatial effects were not accounted for in any way and when biome or biorealm were included as random intercepts in the model

# set up five random folds for cross validation
# set.seed(2345)
# folds = sample(1:5, size = nrow(dat), replace = T)
# 
# compMods_cv = compare_cv(f1, dat, mesh, folds, parallel = F, taxon = taxon, response_var = response_var)
# 
# save(compMods_aic, compMods_cv, fitmesh, file = paste0("results/sdmTMB_models/model_selection/",fname_end,".RData"))
# compMods_cv$compMods_cv

# * - residual check ----------------------------------------------------------

#load(paste0("results/sdmTMB_models/model_selection/",fname_end,".RData"))
out = readRDS("results/sdmTMB_models2/reptiles_sesvert.rds")
bestmod = out$bestmod
plot_resids(mod = bestmod, response_var = "vert.mean.ses", 
            fpath = paste0("figures/residual_checks/sdmTMB2/reptiles/",fname_end), integer_response = F)

# * - plot model coefs for comparison models ---------------------------

#plot_compMods_coefs(mods = compMods_aic$mods, fname = paste0("figures/model_selection/",fname_end,".png"))

# * - predict svc + realm model to the future ---------------------------------------

predict_future(mod = bestmod, newdata = dat.f, type = "response",
               fpath = paste0("results/sdmTMB_models2/predictions/",fname_end,".rds"))



# MEAN VERTICALITY --------------------------------------------------------

# ls = ls()
# a = which(ls == "dat" | ls == "dat.f")
# rm(list = ls[-a])
# source("scripts/00_functions/manuscript_functions.R")
# source("scripts/00_functions/00_plot_functions.R")
# 
# v = vect("data/original/rnaturalearth_world.shp")
# v = project(v, "+proj=cea +datum=WGS84")
# f1 = formula(vert.mean ~ canopy_height + veg_den + I(tmax_warm^2) + tmax_warm + tmin_cold + precip_wet + log_precip_dry + log_clim_velocity)
# 
# # plot relationship between mean verticality and env predictors
# dat %>% 
#   dplyr::select(vert.mean, biome:clim_velocity, elev, veg_den, veg_complexity, log_precip_dry, log_clim_velocity) %>%
#   pivot_longer(cols = 3:17, names_to = "var", values_to = "val") %>% 
#   ggplot(aes(x = val, y = vert.mean, color = biome)) +
#   geom_point(pch = ".") +
#   facet_wrap(~var, scales = "free") +
#   theme_classic()
# 
# 
# 
# taxon = "Reptiles"
# response_var = "Mean verticality"
# fname_end = "reptiles_meanvert"
# 
# # Set up mesh for analysis ------------------------------------------------
# 
# # set up spatial mesh
# 
# # first estimate range of spatial autocorrelation
# samp = dat %>% sample_n(1000)
# samp.cor = ncf::correlog(x = samp$x, y = samp$y, z = samp$vert.mean, increment = 50000/1e5, resamp = 99)
# ncf:::plot.correlog(samp.cor)
# ncf:::plot.correlog(samp.cor, xlim = c(0,100))
# 
# # set initial range as 30 and max.edge as range/5
# # https://haakonbakkagit.github.io/btopic104.html
# # Bakka, H., J. Vanhatalo, J. Illian, D. Simpson, and H. Rue. 2016. “Accounting for Physical Barriers in Species Distribution Modeling with Non-Stationary Spatial Random Effects.” arXiv preprint arXiv:1608.03787. Norwegian University of Science; Technology, Trondheim, Norway. 
# 
# # NOTE: including realm as random effect led to not ok model (sigma_G smaller than 0.01 - consider omitting this part of the model)
# f1 = formula(vert.mean ~ canopy_height + veg_den + I(tmax_warm^2) + tmax_warm + tmin_cold + precip_wet + log_precip_dry + log_clim_velocity)
# 
# fitmesh = fit_mesh(f = f1, dat, range = samp.cor$x.intercept, v = v, family = Beta())
# 
# mesh = fitmesh$meshes[[length(fitmesh$meshes)]]
# 
# save(fitmesh, file = paste0("results/sdmTMB_models/model_selection/", fname_end, ".RData"))
# load(paste0("results/sdmTMB_models/model_selection/", fname_end, ".RData"))
# 
# 
# # * - Compare models with AIC -------------------------------------------------
# 
# # Compare models using all points with AIC
# 
# compMods_aic = compareMods_AIC(f = f1, dat, mesh, taxon = taxon, response_var = response_var, family = Beta(link = "logit"))
# sanity(compMods_aic$mods[[1]])
# sanity(compMods_aic$mods[[2]])
# sanity(compMods_aic$mods[[3]])
# sanity(compMods_aic$mods[[4]])
# 
# compMods_aic[[2]]
# save(compMods_aic, fitmesh, file = paste0( "results/sdmTMB_models/model_selection/", fname_end, ".RData"))
# #load(paste0("results/sdmTMB_models/model_selection/", fname_end, ".RData"))
# 
# compMods_aic$compare
# 
# 
# # realm does not contribute to the model, but maybe best to keep it for consistency
# 
# # * - Cross validation --------------------------------------------------------
# 
# # compare random effects structure using cross validation
# # random effects include random intercept, spatial random field, and spatially varying coefficient
# # spatial random field and spatially varying coefficient both depend on mesh
# # preliminary assessment indicated that residuals displayed strong spatial autocorrelation when
# # spatial effects were not accounted for in any way and when biome or biorealm were included as random intercepts in the model
# 
# # set up five random folds for cross validation
# set.seed(2345)
# folds = sample(1:5, size = nrow(dat), replace = T)
# 
# compMods_cv = compare_cv_beta(f1, dat, mesh, folds, parallel = F, taxon = taxon, response_var = response_var)
# 
# load("tempfiles/amph_mods.RData")
# load("tempfiles/amph_mods2.RData")
# load("tempfiles/amph_mods3.RData")
# 
# lapply(compMods_cv$mods$mod.realm.cv$models, sanity)
# 
# save(compMods_aic, compMods_cv, fitmesh, file = paste0("results/sdmTMB_models/model_selection/", fname_end, ".RData"))
# 
# compMods_cv$compMods_cv
# 
# # * - residual check ----------------------------------------------------------
# 
# load(paste0("results/sdmTMB_models/model_selection/", fname_end, ".RData"))
# 
# plot_resids(mod = compMods_aic$mods$mod.realm.svc, response_var = "vert.mean", 
#             fpath = paste0("figures/residual_checks/", fname_end))
# 
# # * - plot model coefs for comparison models ---------------------------
# 
# plot_compMods_coefs(mods = compMods_aic$mods, fname = paste0("figures/model_selection/", fname_end, ".png"))
# 
# # * - predict svc + realm model to the future ---------------------------------------
# 
# predict_future(mod = compMods_aic$mods$mod.realm.svc, newdata = dat.f, type = "response",
#                fpath = paste0("results/sdmTMB_models/", fname_end,".RData"))
# 
# 
# 
# 
# # PROPORTION ARBOREAL --------------------------------------------------------
# 
# ls = ls()
# a = which(ls == "dat" | ls == "dat.f")
# rm(list = ls[-a])
# 
# source("scripts/00_functions/manuscript_functions.R")
# source("scripts/00_functions/00_plot_functions.R")
# 
# v = vect("data/original/rnaturalearth_world.shp")
# v = project(v, "+proj=cea +datum=WGS84")
# f1 = formula(p.arb ~ canopy_height + veg_den + I(tmax_warm^2) + tmax_warm + tmin_cold + precip_wet + log_precip_dry + log_clim_velocity)
# 
# # plot relationship between proportion arboreality and env predictors
# dat %>% 
#   dplyr::select(p.arb, biome:clim_velocity, elev, veg_den, veg_complexity, log_precip_dry, log_clim_velocity) %>%
#   pivot_longer(cols = 3:17, names_to = "var", values_to = "val") %>% 
#   ggplot(aes(x = val, y = p.arb, color = biome)) +
#   geom_point(pch = ".") +
#   facet_wrap(~var, scales = "free") +
#   theme_classic()
# 
# 
# 
# taxon = "Reptiles"
# response_var = "Proportion Arboreal"
# fname_end = "reptiles_parb"
# wts = dat$rich
# 
# # Set up mesh for analysis ------------------------------------------------
# 
# # set up spatial mesh
# 
# # first estimate range of spatial autocorrelation
# samp = dat %>% sample_n(1000)
# samp.cor = ncf::correlog(x = samp$x, y = samp$y, z = samp$vert.mean, increment = 50000/1e5, resamp = 99)
# ncf:::plot.correlog(samp.cor)
# ncf:::plot.correlog(samp.cor, xlim = c(0,100), ylim = c(-1,1))
# 
# # set initial range as 30 and max.edge as range/5
# # https://haakonbakkagit.github.io/btopic104.html
# # Bakka, H., J. Vanhatalo, J. Illian, D. Simpson, and H. Rue. 2016. “Accounting for Physical Barriers in Species Distribution Modeling with Non-Stationary Spatial Random Effects.” arXiv preprint arXiv:1608.03787. Norwegian University of Science; Technology, Trondheim, Norway. 
# 
# dat$narb = dat$rich*dat$p.arb
# dat$notarb = dat$rich*(1-dat$p.arb)
# f1 = formula(cbind(narb,notarb) ~ canopy_height + veg_den + I(tmax_warm^2) + tmax_warm + tmin_cold + precip_wet + log_precip_dry + log_clim_velocity, weights = rich)
# 
# 
# fitmesh = fit_mesh(f = f1, dat, range = samp.cor$x.intercept, v = v, family = binomial())
# 
# mesh = fitmesh$meshes[[length(fitmesh$meshes)]]
# 
# save(fitmesh, file = paste0("results/sdmTMB_models/model_selection/", fname_end, ".RData"))
# load(paste0("results/sdmTMB_models/model_selection/", fname_end, ".RData"))
# 
# 
# # * - Compare models with AIC -------------------------------------------------
# 
# # Compare models using all points with AIC
# 
# compMods_aic = compareMods_AIC(f = f1, dat, mesh, taxon = taxon, response_var = response_var, family = binomial())
# sanity(compMods_aic$mods[[1]])
# sanity(compMods_aic$mods[[2]])
# sanity(compMods_aic$mods[[3]])
# sanity(compMods_aic$mods[[4]])
# 
# compMods_aic[[2]]
# save(compMods_aic, fitmesh, file = paste0( "results/sdmTMB_models/model_selection/", fname_end, ".RData"))
# load(paste0("results/sdmTMB_models/model_selection/", fname_end, ".RData"))
# 
# compMods_aic$compare
# 
# 
# # realm does not contribute to the model, but maybe best to keep it for consistency
# 
# # * - Cross validation --------------------------------------------------------
# 
# # compare random effects structure using cross validation
# # random effects include random intercept, spatial random field, and spatially varying coefficient
# # spatial random field and spatially varying coefficient both depend on mesh
# # preliminary assessment indicated that residuals displayed strong spatial autocorrelation when
# # spatial effects were not accounted for in any way and when biome or biorealm were included as random intercepts in the model
# 
# # set up five random folds for cross validation
# set.seed(2345)
# folds = sample(1:5, size = nrow(dat), replace = T)
# 
# compMods_cv = compare_cv_binom(f1, dat, mesh, folds, parallel = F, taxon = taxon, response_var = response_var)
# 
# lapply(compMods_cv$mods$mod.realm.cv$models, sanity)
# 
# save(compMods_aic, compMods_cv, fitmesh, file = paste0("results/sdmTMB_models/model_selection/", fname_end, ".RData"))
# 
# compMods_cv$compMods_cv
# 
# # * - residual check ----------------------------------------------------------
# # issues with this for parb
# load(paste0("results/sdmTMB_models/model_selection/", fname_end, ".RData"))
# 
# plot_resids(mod = compMods_aic$mods$mod.realm.svc, response_var = "p.arb", 
#             fpath = paste0("figures/residual_checks/", fname_end))
# 
# # * - plot model coefs for comparison models ---------------------------
# 
# plot_compMods_coefs(mods = compMods_aic$mods, fname = paste0("figures/model_selection/", fname_end, ".png"))
# 
# # * - predict svc + realm model to the future ---------------------------------------
# 
# predict_future(mod = compMods_aic$mods$mod.realm.svc, newdata = dat.f, type = "response",
#                fpath = paste0("results/sdmTMB_models/", fname_end,".RData"))




# # Cross validation --------------------------------------------------------
# 
# # compare random effects structure using cross validation
# # preliminary assessment indicated that residuals displayed strong spatial autocorrelation when
# # spatial effects were not accounted for in any way and when biome or biorealm were included as random intercepts in the model
# 
# # set up spatial mesh
# 
# # use 300km mesh - relatively fine with reasonable fitting speed
# mesh300km = make_mesh(data = dat, xy_cols = c("x", "y"), cutoff = 3)
# plot(mesh300km)
# 
# # set up five random folds for cross validation
# folds = sample(1:5, size = nrow(dat), replace = T)
# 
# # 1: fit random intercept of ecoregion without spatial random field
# library(future)
# plan(multisession)
# 
# mod.ecoregion.cv = sdmTMB_cv(update(f1, ~ . + (1|ecoregion)), 
#                              data = dat,
#                              mesh = mesh300km,
#                              spatial = "off",
#                              reml = T, 
#                              fold_ids = folds,
#                              k_folds = length(unique(folds)))
# mod.ecoregion.cv$converged
# lapply(mod.ecoregion.cv$models, sanity)
# # gradients not okay for folds 2 and 5
# 
# mod.mesh.cv = sdmTMB_cv(f1, 
#                         data = dat,
#                         mesh = mesh300km,
#                         spatial = "on",
#                         reml = T, 
#                         fold_ids = folds,
#                         k_folds = length(unique(folds)))
# 
# mod.mesh.cv$converged
# lapply(mod.mesh.cv$models, sanity)
# # all okay
# 
# mod.mesh.biome.cv = sdmTMB_cv(update(f1, ~ . +(1|biome)), 
#                               data = dat,
#                               mesh = mesh300km,
#                               spatial = "on",
#                               reml = T, 
#                               fold_ids = folds,
#                               k_folds = length(unique(folds)))
# 
# mod.mesh.biome.cv$converged
# lapply(mod.mesh.biome.cv$models, sanity)
# # first model - gradients not okay, other models - all okay
# 
# mod.mesh.biorealm.cv = sdmTMB_cv(update(f1, ~ . +(1|biorealm)), 
#                                  data = dat,
#                                  mesh = mesh300km,
#                                  spatial = "on",
#                                  reml = T, 
#                                  fold_ids = folds,
#                                  k_folds = length(unique(folds)))
# 
# mod.mesh.biorealm.cv$converged
# lapply(mod.mesh.biorealm.cv$models, sanity)
# # all okay
# 
# mod.mesh.ecoregion.cv = sdmTMB_cv(update(f1, ~ . +(1|ecoregion)), 
#                                   data = dat,
#                                   mesh = mesh300km,
#                                   spatial = "on",
#                                   reml = T, 
#                                   fold_ids = folds,
#                                   k_folds = length(unique(folds)))
# 
# mod.mesh.ecoregion.cv$converged
# lapply(mod.mesh.ecoregion.cv$models, sanity)
# # 2 gradient not okay
# 
# # compare sum_logLik
# mod.ecoregion.cv$sum_loglik
# mod.mesh.cv$sum_loglik
# mod.mesh.biome.cv$sum_loglik
# mod.mesh.biorealm.cv$sum_loglik
# mod.mesh.ecoregion.cv$sum_loglik # best model for random effects structure
# 
# mod.mesh.spatialvar.cv = sdmTMB_cv(update(f1, ~ . +(1|ecoregion)), 
#                                    data = dat,
#                                    mesh = mesh300km,
#                                    spatial = "on",
#                                    reml = T, 
#                                    fold_ids = folds,
#                                    spatial_varying = ~ 0 + canopy_height,
#                                    k_folds = length(unique(folds)))
# mod.mesh.spatialvar.cv$converged
# lapply(mod.mesh.spatialvar.cv$models, sanity)
# 
# mod.mesh.ecoregion.cv$sum_loglik
# mod.mesh.spatialvar.cv$sum_loglik # best
# 
# 
# summ_coefs = function(x, nm){
#   mods = x$models
#   test = lapply(mods, sanity)
#   test = sapply(test, "[[", 9) # get T/F for all_ok
#   mods = mods[which(test)] # subset mods to okay models
#   mods = lapply(mods, tidy)
#   for (i in 1:length(mods)) {mods[[i]]$mod = as.character(i)}
#   mods = rbindlist(mods)
#   mods = mods %>% 
#     group_by(term) %>% 
#     summarise(est = mean(estimate), se = mean(std.error)) %>% 
#     mutate(modname = nm)
#   return(mods)
# }
# 
# mod.ecoregion.summ = summ_coefs(mod.ecoregion.cv, "mod.ecoregion")
# mod.mesh.summ = summ_coefs(mod.mesh.cv, "mod.mesh")
# mod.mesh.biome.summ = summ_coefs(mod.mesh.biome.cv, "mod.mesh.biome")
# mod.mesh.biorealm.summ = summ_coefs(mod.mesh.biorealm.cv, "mod.mesh.biorealm")
# mod.mesh.ecoregion.summ = summ_coefs(mod.mesh.ecoregion.cv, "mod.mesh.ecoregion")
# mod.mesh.spatialvar.summ = summ_coefs(mod.mesh.spatialvar.cv, "mod.mesh.spatialvar")
# 
# mods = rbind(mod.ecoregion.summ, mod.mesh.summ, mod.mesh.biome.summ, mod.mesh.biorealm.summ, mod.mesh.ecoregion.summ, mod.mesh.spatialvar.summ)
# ggplot(mods) +
#   geom_pointrange(aes(x = est, xmax = est+se, xmin = est-se, y = term, color = modname),
#                   position = position_dodge2(width=0.5))
# 
# 
# 
# # Compare models with different fixed effects -----------------------------
# 
# m.full = sdmTMB_cv(formula = update(f1, ~ . +(1|ecoregion)),
#                    data = dat,
#                    mesh = mesh300km,
#                    spatial = "on",
#                    reml = F, 
#                    spatial_varying = ~0 + canopy_height,
#                    fold_ids = folds,
#                    k_folds = length(unique(folds)))
# m.full$converged
# lapply(m.full$models, sanity)
# # all okay
# 
# 
# m.nopoly = sdmTMB_cv(formula = update(f1, ~ . -I(tmax_warm^2) +(1|ecoregion)),
#                      data = dat,
#                      mesh = mesh300km,
#                      spatial = "on",
#                      reml = F,
#                      spatial_varying = ~0 + canopy_height,
#                      fold_ids = folds,
#                      k_folds = length(unique(folds)))
# m.nopoly$converged
# lapply(m.nopoly$models, sanity)
# 
# m.noclimvel = sdmTMB_cv(formula = update(f1, ~ . -log_clim_velocity +(1|ecoregion)),
#                         data = dat,
#                         mesh = mesh300km,
#                         spatial = "on",
#                         reml = F,
#                         spatial_varying = ~0 + canopy_height,
#                         fold_ids = folds,
#                         k_folds = length(unique(folds)))
# m.noclimvel$converged
# lapply(m.noclimvel$models, sanity)
# 
# m.noclimvel_nopoly = sdmTMB_cv(formula = update(f1, ~ . -I(tmax_warm^2) -log_clim_velocity +(1|ecoregion)),
#                         data = dat,
#                         mesh = mesh300km,
#                         spatial = "on",
#                         reml = F,
#                         spatial_varying = ~0 + canopy_height,
#                         fold_ids = folds,
#                         k_folds = length(unique(folds)))
# 
# m.noclimvel_nopoly$converged
# lapply(m.noclimvel_nopoly$models, sanity)
# 
# 
# m.full$sum_loglik # performs best
# m.nopoly$sum_loglik
# m.noclimvel$sum_loglik
# m.noclimvel_nopoly$sum_loglik
# 
# 
# 
# # Refit best model using all the data ------------------------------------------
# control = sdmTMBcontrol(nlminb_loops = 2)
# m.final = sdmTMB(formula = update(f1, . ~ . + (1|ecoregion)),
#                  data = dat,
#                  mesh = mesh300km,
#                  spatial = "on",
#                  spatial_varying = ~ 0 + canopy_height,
#                  reml = T,
#                  control = control)
# 
# m.final
# sanity(m.final)
# 
# pred <- predict(m.final)
# s_m.final = simulate(m.final, nsim = 500, seed = 12345, type = "mle-mvn")
# r_m.final <- DHARMa::createDHARMa(
#   simulatedResponse = s_m.final,
#   observedResponse = dat$vert.mean.ses,
#   fittedPredictedResponse = pred$est_non_rf
# )
# 
# r_m.final <- DHARMa::createDHARMa(
#   simulatedResponse = s_m.final,
#   observedResponse = dat$vert.mean.ses,
#   fittedPredictedResponse = pred$est
# )
# 
# plot(r_m.final)
# testQuantiles(r_m.final)
# testResiduals(r_m.final)
# plotResiduals(r_m.final, rank = F, form = dat$vert.mean.ses)
# plot(dat$vert.mean.ses, pred$est_non_rf)
# plot(dat$vert.mean.ses, pred$est)
# #testSpatialAutocorrelation(r_m.final, x = dat$x, y = dat$y)
# 
# # need to subset predictions to calculate spatial autocorrelation - otherwise there are too many points
# sub = sample(1:nrow(s_m.final), 10000)
# r_m.final.sub = DHARMa::createDHARMa(
#   simulatedResponse = s_m.final[sub,],
#   observedResponse = dat$vert.mean.ses[sub],
#   fittedPredictedResponse = pred$est_non_rf[sub]
# )
# testSpatialAutocorrelation(r_m.final.sub, x = dat$x[sub], y = dat$y[sub])
# 
# 
# 
# # predict models to future env data ---------------------------------------
# dat.f$ecoregion = factor(dat.f$ecoregion)
# dat.f$x = dat.f$x/1e5
# dat.f$y = dat.f$y/1e5
# 
# pred.f = predict(m.final, newdata = dat.f, type = "response")
# head(pred.f)
# 
# pred = predict(m.final, type = "response")
# 
# pred.f$est.dif = pred.f$est - pred$est
# 
# 
# ggplot() +
#   geom_tile(data = pred.f, aes(x*1e5, y*1e5, fill = est.dif)) +
#   scale_fill_continuous_divergingx("spectral") +
#   coord_sf(crs = "+proj=cea +datum=WGS84")
# 
# save(m.final, pred, pred.f, file = "results/sdmTMB_models/reptiles_sesvert.RData")
# 
# 
# # MEAN VERTICALITY --------------------------------------------------------
# 
# head(dat)
# ls = ls()
# rm(list = ls[-c(2,3)])
# 
# # plot relationship between mean verticality and env predictors
# dat %>% 
#   dplyr::select(vert.mean, biome:clim_velocity, elev, veg_den, veg_complexity, log_precip_dry, log_clim_velocity) %>%
#   pivot_longer(cols = 3:18, names_to = "var", values_to = "val") %>% 
#   ggplot(aes(x = val, y = vert.mean, color = biome)) +
#   geom_point(pch = ".") +
#   facet_wrap(~var, scales = "free") +
#   theme_classic()
# 
# # dat$ecoregion = factor(dat$ecoregion)
# # dat$realm = factor(dat$realm)
# # dat$biome = factor(dat$biome)
# # dat$biorealm = factor(paste(dat$biome, dat$realm, spe = "_"))
# # dat$x = dat$x/1e5
# # dat$y = dat$y/1e5
# 
# mesh300km = make_mesh(data = dat, xy_cols = c("x", "y"), cutoff = 3)
# plot(mesh300km)
# 
# # set up five random folds for cross validation
# folds = sample(1:5, size = nrow(dat), replace = T)
# 
# f1 = formula(vert.mean ~ canopy_height + veg_den + I(tmax_warm^2) +tmax_warm + tmin_cold + precip_wet + log_precip_dry + log_clim_velocity)
# 
# # 1: fit random intercept of ecoregion without spatial random field
# library(future)
# plan(multisession)
# 
# mod.ecoregion.cv = sdmTMB_cv(update(f1, ~ . + (1|ecoregion)), 
#                              data = dat,
#                              mesh = mesh300km,
#                              spatial = "off",
#                              family = Beta(),
#                              reml = T, 
#                              fold_ids = folds,
#                              k_folds = length(unique(folds)))
# mod.ecoregion.cv$converged
# lapply(mod.ecoregion.cv$models, sanity)
# 
# 
# 
# mod.mesh.cv = sdmTMB_cv(f1, 
#                         data = dat,
#                         mesh = mesh300km,
#                         spatial = "on",
#                         family = Beta(),
#                         reml = T, 
#                         fold_ids = folds,
#                         k_folds = length(unique(folds)))
# 
# mod.mesh.cv$converged #true
# lapply(mod.mesh.cv$models, sanity)
# # gradient issues with 1 and 4
# # all okay
# 
# mod.mesh.biome.cv = sdmTMB_cv(update(f1, ~ . +(1|biome)), 
#                               data = dat,
#                               mesh = mesh300km,
#                               spatial = "on",
#                               family = Beta(),
#                               reml = T, 
#                               fold_ids = folds,
#                               k_folds = length(unique(folds)))
# 
# mod.mesh.biome.cv$converged
# lapply(mod.mesh.biome.cv$models, sanity)
# # all models have issues
# 
# mod.mesh.biorealm.cv = sdmTMB_cv(update(f1, ~ . +(1|biorealm)), 
#                                  data = dat,
#                                  mesh = mesh300km,
#                                  spatial = "on",
#                                  family = Beta(),
#                                  reml = T, 
#                                  fold_ids = folds,
#                                  k_folds = length(unique(folds)))
# 
# mod.mesh.biorealm.cv$converged 
# lapply(mod.mesh.biorealm.cv$models, sanity)
# # issue with gradient in 3 - all others okay
# 
# mod.mesh.ecoregion.cv = sdmTMB_cv(update(f1, ~ . +(1|ecoregion)), 
#                                   data = dat,
#                                   mesh = mesh300km,
#                                   spatial = "on",
#                                   family = Beta(),
#                                   reml = T, 
#                                   fold_ids = folds,
#                                   k_folds = length(unique(folds)))
# 
# mod.mesh.ecoregion.cv$converged #false
# lapply(mod.mesh.ecoregion.cv$models, sanity)
# # all okay
# 
# # compare sum_logLik
# mod.ecoregion.cv$sum_loglik
# mod.mesh.cv$sum_loglik
# mod.mesh.biome.cv$sum_loglik
# mod.mesh.biorealm.cv$sum_loglik
# mod.mesh.ecoregion.cv$sum_loglik #best
# 
# mod.mesh.spatialvar.cv = sdmTMB_cv(update(f1, ~ . +(1|ecoregion)), 
#                                   data = dat,
#                                   mesh = mesh300km,
#                                   spatial = "on",
#                                   family = Beta(),
#                                   reml = T, 
#                                   fold_ids = folds,
#                                   spatial_varying = ~ 0 + canopy_height,
#                                   k_folds = length(unique(folds)))
# 
# mod.mesh.spatialvar.cv$converged
# lapply(mod.mesh.spatialvar.cv$models, sanity)
# 
# mod.mesh.ecoregion.cv$sum_loglik
# mod.mesh.spatialvar.cv$sum_loglik #best
# 
# 
# summ_coefs = function(x, nm){
#   mods = x$models
#   test = lapply(mods, sanity)
#   test = sapply(test, "[[", 9) # get T/F for all_ok
#   mods = mods[which(test)] # subset mods to okay models
#   mods = lapply(mods, tidy)
#   for (i in 1:length(mods)) {mods[[i]]$mod = as.character(i)}
#   mods = rbindlist(mods)
#   mods = mods %>% 
#     group_by(term) %>% 
#     summarise(est = mean(estimate), se = mean(std.error)) %>% 
#     mutate(modname = nm)
#   return(mods)
# }
# 
# mod.ecoregion.summ = summ_coefs(mod.ecoregion.cv, "mod.ecoregion")
# mod.mesh.summ = summ_coefs(mod.mesh.cv, "mod.mesh")
# #mod.mesh.biome.summ = summ_coefs(mod.mesh.biome.cv, "mod.mesh.biome")
# mod.mesh.biorealm.summ = summ_coefs(mod.mesh.biorealm.cv, "mod.mesh.biorealm")
# mod.mesh.ecoregion.summ = summ_coefs(mod.mesh.ecoregion.cv, "mod.mesh.ecoregion")
# mod.mesh.spatialvar.summ = summ_coefs(mod.mesh.spatialvar.cv, "mod.mesh.spatialvar")
# 
# mods = rbind(mod.ecoregion.summ, mod.mesh.summ, mod.mesh.biorealm.summ, mod.mesh.ecoregion.summ, mod.mesh.spatialvar.summ)
# ggplot(mods %>% filter(term != "(Intercept)")) +
#   geom_pointrange(aes(x = est, xmax = est+se*1.96, xmin = est-se*1.96, y = term, color = modname),
#                   position = position_dodge2(width=0.5)) +
#   geom_vline(xintercept = 0, linetype = "dashed")
# 
# 
# 
# # Compare models with different fixed effects -----------------------------
# m.full = sdmTMB_cv(formula = update(f1, ~ . +(1|ecoregion)),
#                    data = dat,
#                    mesh = mesh300km,
#                    spatial = "on",
#                    family = Beta(),
#                    reml = F, 
#                    spatial_varying = ~ 0 + canopy_height,
#                    fold_ids = folds,
#                    k_folds = length(unique(folds)))
# m.full$converged
# lapply(m.full$models, sanity)
# # all okay
# 
# 
# m.nopoly = sdmTMB_cv(formula = update(f1, ~ . -I(tmax_warm^2) +(1|ecoregion)),
#                      data = dat,
#                      mesh = mesh300km,
#                      spatial = "on",
#                      family = Beta(),
#                      reml = F, 
#                      spatial_varying = ~ 0 + canopy_height,
#                      fold_ids = folds,
#                      k_folds = length(unique(folds)))
# m.nopoly$converged
# lapply(m.nopoly$models, sanity)
# # all okay
# 
# 
# 
# m.noclimvel = sdmTMB_cv(formula = update(f1, ~ . -log_clim_velocity +(1|ecoregion)),
#                         data = dat,
#                         mesh = mesh300km,
#                         spatial = "on",
#                         family = Beta(),
#                         reml = F, 
#                         spatial_varying = ~ 0 + canopy_height,
#                         fold_ids = folds,
#                         k_folds = length(unique(folds)))
# m.noclimvel$converged
# lapply(m.noclimvel$models, sanity)
# 
# m.noclimvel_nopoly = sdmTMB_cv(formula = update(f1, ~ . -I(tmax_warm^2) -log_clim_velocity +(1|ecoregion)),
#                         data = dat,
#                         mesh = mesh300km,
#                         spatial = "on",
#                         family = Beta(),
#                         reml = F, 
#                         spatial_varying = ~ 0 + canopy_height,
#                         fold_ids = folds,
#                         k_folds = length(unique(folds)))
# m.noclimvel_nopoly$converged
# lapply(m.noclimvel_nopoly$models, sanity)
# 
# m.full$sum_loglik  # performs best
# m.nopoly$sum_loglik
# m.noclimvel$sum_loglik
# m.noclimvel_nopoly$sum_loglik
# 
# 
# # Refit best model using all the data ------------------------------------------
# control = sdmTMBcontrol(nlminb_loops = 2)
# m.final = sdmTMB(formula = update(f1, . ~ . + (1|ecoregion)),
#                  data = dat,
#                  mesh = mesh300km,
#                  family = Beta(),
#                  spatial = "on",
#                  spatial_varying = ~ 0 + canopy_height,
#                  reml = T,
#                  control = sdmTMBcontrol(nlminb_loops = 2, iter.max = 2000, eval.max = 4000, start = list(ln_phi = 7.5)))
# 
# 
# m.final
# sanity(m.final)
# 
# sanity(m.full$models[[1]])
# sanity(m.full$models[[2]])
# sanity(m.full$models[[3]])
# sanity(m.full$models[[4]])
# sanity(m.full$models[[5]])
# test1 = m.full$models[[1]]$tmb_obj$env$parList()
# test2 = m.full$models[[2]]$tmb_obj$env$parList()
# test3 = m.full$models[[3]]$tmb_obj$env$parList()
# test4 = m.full$models[[4]]$tmb_obj$env$parList()
# test5 = m.full$models[[5]]$tmb_obj$env$parList()
# test1$ln_phi
# test2$ln_phi
# test3$ln_phi
# test4$ln_phi
# test5$ln_phi
# 
# pred <- predict(m.final)
# s_m.final = simulate(m.final, nsim = 500, seed = 12345, type = "mle-mvn")
# r_m.final <- DHARMa::createDHARMa(
#   simulatedResponse = s_m.final,
#   observedResponse = dat$vert.mean,
#   fittedPredictedResponse = pred$est_non_rf
# )
# 
# r_m.final <- DHARMa::createDHARMa(
#   simulatedResponse = s_m.final,
#   observedResponse = dat$vert.mean,
#   fittedPredictedResponse = pred$est
# )
# 
# plot(r_m.final)
# testQuantiles(r_m.final)
# testResiduals(r_m.final)
# plotResiduals(r_m.final, rank = F, form = dat$vert.mean)
# plot(dat$vert.mean, m.final$family$linkinv(pred$est_non_rf))
# abline(a = 0, b = 1, col = "red")
# plot(dat$vert.mean, m.final$family$linkinv(pred$est))
# abline(a = 0, b = 1, col = "red")
# #testSpatialAutocorrelation(r_m.final, x = dat$x, y = dat$y)
# 
# # need to subset predictions to calculate spatial autocorrelation - otherwise there are too many points
# sub = sample(1:nrow(s_m.final), 10000)
# r_m.final.sub = DHARMa::createDHARMa(
#   simulatedResponse = s_m.final[sub,],
#   observedResponse = dat$vert.mean[sub],
#   fittedPredictedResponse = pred$est_non_rf[sub]
# )
# testSpatialAutocorrelation(r_m.final.sub, x = dat$x[sub], y = dat$y[sub])
# 
# 
# pred.f = predict(m.final, newdata = dat.f, type = "response")
# head(pred.f)
# 
# # predict in response scale
# pred = predict(m.final, type = "response")
# 
# pred.f$est.dif = pred.f$est - pred$est
# 
# ggplot() +
#   geom_tile(data = pred.f, aes(x*1e5, y*1e5, fill = est.dif)) +
#   scale_fill_continuous_divergingx("spectral") +
#   coord_sf(crs = "+proj=cea +datum=WGS84")
# 
# save(m.final, pred, pred.f, file = "results/sdmTMB_models/reptiles_meanvert.RData")
