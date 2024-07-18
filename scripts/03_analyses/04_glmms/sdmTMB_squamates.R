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

# squamate forest only 50km resolution
raw = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/squamates_comdat.csv")

# remove cells with fewer than five species
raw = raw %>% 
  # subset to richness >= 5
  filter(rich >= 5) %>% 
  filter(!is.na(vert.mean.ses)) %>% 
  drop_na()

# future env data
env.f = rast(list.files(path = "data/derivative_data/resampled_env_rasters_50km/chelsa_future/2071_2100/ensemble/",
                        pattern = ".tif", full.names = T))
env.f = env.f[[2:6]]

dat.f = raw %>% 
  dplyr::select(x,y,canopy_height, veg_den, clim_velocity, ecoregion)

env.f = terra::extract(env.f, dat.f[,c('x', "y")], ID = F)

dat.f = cbind(dat.f, env.f)
dat.f = dat.f %>% rename(tmax_warm = max_temp_warm, tmin_cold = min_temp_cold)

# remove NA values from future and present data
raw = raw[which(!is.na(dat.f$precip_dry)),]
dat.f = dat.f[which(!is.na(dat.f$precip_dry)),]


# scale the data
dat = raw %>%
  mutate_at(.vars = vars(canopy_height:clim_velocity, elev, veg_den, veg_complexity), .funs = scale)

# plot relationship between SES verticality and env predictors
dat %>% 
  dplyr::select(vert.mean.ses, biome:clim_velocity, elev, veg_den, veg_complexity) %>%
  pivot_longer(cols = 3:16, names_to = "var", values_to = "val") %>% 
  ggplot(aes(x = val, y = vert.mean.ses, color = biome)) +
  geom_point(pch = ".") +
  facet_wrap(~var, scales = "free") +
  theme_classic()

# scale future data
dat.f = dat.f %>% 
  mutate(tmax_warm = (tmax_warm - mean(raw$tmax_warm))/sd(raw$tmax_warm),
         tmin_cold = (tmin_cold-mean(raw$tmin_cold))/sd(raw$tmin_cold),
         precip_dry = (precip_dry - mean(raw$precip_dry))/sd(raw$precip_dry),
         precip_wet = (precip_wet - mean(raw$precip_wet))/sd(raw$precip_wet),
         precip_sea = (precip_sea - mean(raw$precip_sea))/sd(raw$precip_sea),
         canopy_height = (canopy_height - mean(raw$canopy_height))/sd(raw$canopy_height),
         veg_den = (veg_den - mean(raw$veg_den))/sd(raw$veg_den),
         clim_velocity = (clim_velocity - mean(raw$clim_velocity))/sd(raw$clim_velocity))       

# test correlation between predictor variables
dat %>% 
  dplyr::select(canopy_height, veg_den, veg_complexity,
                tmax_warm, tmin_cold, temp_sea,
                precip_wet, precip_dry, precip_sea,
                clim_velocity) %>% 
  pairs()

# veg_complexity strongly correlated with canopy height and veg den
# temp sea strongly correlated with tmin_cold

vif(lm(vert.mean.ses ~ canopy_height + veg_den + veg_complexity + tmax_warm + tmin_cold + temp_sea + precip_sea + precip_wet + precip_dry + clim_velocity, data = dat))

# remove veg_complexity
vif(lm(vert.mean.ses ~ canopy_height + veg_den + tmax_warm + tmin_cold + temp_sea + precip_sea + precip_wet + precip_dry + clim_velocity, data = dat))

# remove temp_sea
vif(lm(vert.mean.ses ~ canopy_height + veg_den + tmax_warm + tmin_cold + precip_sea + precip_wet + precip_dry + clim_velocity, data = dat))

# VIF all under 4

f1 = formula(vert.mean.ses ~ canopy_height + veg_den + tmax_warm + tmin_cold + precip_sea + precip_wet + precip_dry + clim_velocity)


# Cross validation --------------------------------------------------------

# compare random effects structure using cross validation
# preliminary assessment indicated that residuals displayed strong spatial autocorrelation when
# spatial effects were not accounted for in any way and when biome or biorealm were included as random intercepts in the model

# set up spatial mesh
dat$x = dat$x/1e5
dat$y = dat$y/1e5
dat$ecoregion = factor(dat$ecoregion)
dat$biome = factor(dat$biome)
dat$biorealm = factor(paste(dat$biome, dat$realm, sep = "_"))


# use 300km mesh - relatively fine with reasonable fitting speed
mesh300km = make_mesh(data = dat, xy_cols = c("x", "y"), cutoff = 3)
plot(mesh300km)

# set up five random folds for cross validation
folds = sample(1:5, size = nrow(dat), replace = T)

# 1: fit random intercept of ecoregion without spatial random field
library(future)
plan(multisession)

mod.ecoregion.cv = sdmTMB_cv(update(f1, ~ . + (1|ecoregion)), 
                             data = dat,
                             mesh = mesh300km,
                             spatial = "off",
                             reml = T, 
                             fold_ids = folds,
                             k_folds = length(unique(folds)))
mod.ecoregion.cv$converged
lapply(mod.ecoregion.cv$models, sanity)
# all okay

mod.mesh.cv = sdmTMB_cv(f1, 
                        data = dat,
                        mesh = mesh300km,
                        spatial = "on",
                        reml = T, 
                        fold_ids = folds,
                        k_folds = length(unique(folds)))

mod.mesh.cv$converged
lapply(mod.mesh.cv$models, sanity)
# all okay

mod.mesh.biome.cv = sdmTMB_cv(update(f1, ~ . +(1|biome)), 
                              data = dat,
                              mesh = mesh300km,
                              spatial = "on",
                              reml = T, 
                              fold_ids = folds,
                              k_folds = length(unique(folds)),
                              parallel = F)

mod.mesh.biome.cv$converged
lapply(mod.mesh.biome.cv$models, sanity)
# all okay except #2 = gradients not okay

mod.mesh.biorealm.cv = sdmTMB_cv(update(f1, ~ . +(1|biorealm)), 
                                 data = dat,
                                 mesh = mesh300km,
                                 spatial = "on",
                                 reml = T, 
                                 fold_ids = folds,
                                 k_folds = length(unique(folds)),
                                 parallel = F)

mod.mesh.biorealm.cv$converged
lapply(mod.mesh.biorealm.cv$models, sanity)
# all okay

mod.mesh.ecoregion.cv = sdmTMB_cv(update(f1, ~ . +(1|ecoregion)), 
                                  data = dat,
                                  mesh = mesh300km,
                                  spatial = "on",
                                  reml = T, 
                                  fold_ids = folds,
                                  k_folds = length(unique(folds)))

mod.mesh.ecoregion.cv$converged
lapply(mod.mesh.ecoregion.cv$models, sanity)
# all okay


# compare sum_logLik
mod.ecoregion.cv$sum_loglik
mod.mesh.cv$sum_loglik
mod.mesh.biome.cv$sum_loglik
mod.mesh.biorealm.cv$sum_loglik
mod.mesh.ecoregion.cv$sum_loglik # performed best

summ_coefs = function(x, nm){
  mods = x$models
  test = lapply(mods, sanity)
  test = sapply(test, "[[", 9) # get T/F for all_ok
  mods = mods[which(test)] # subset mods to okay models
  mods = lapply(mods, tidy)
  for (i in 1:length(mods)) {mods[[i]]$mod = as.character(i)}
  mods = rbindlist(mods)
  mods = mods %>% 
    group_by(term) %>% 
    summarise(est = mean(estimate), se = mean(std.error)) %>% 
    mutate(modname = nm)
  return(mods)
}

mod.ecoregion.summ = summ_coefs(mod.ecoregion.cv, "mod.ecoregion")
mod.mesh.summ = summ_coefs(mod.mesh.cv, "mod.mesh")
mod.mesh.biome.summ = summ_coefs(mod.mesh.biome.cv, "mod.mesh.biome")
mod.mesh.biorealm.summ = summ_coefs(mod.mesh.biorealm.cv, "mod.mesh.biorealm")
mod.mesh.ecoregion.summ = summ_coefs(mod.mesh.ecoregion.cv, "mod.mesh.ecoregion")

mods = rbind(mod.ecoregion.summ, mod.mesh.summ, mod.mesh.biome.summ, mod.mesh.biorealm.summ, mod.mesh.ecoregion.summ)
ggplot(mods %>% filter(term != "(Intercept)")) +
  geom_pointrange(aes(x = est, xmax = est+se, xmin = est-se, y = term, color = modname),
                  position = position_dodge2(width=0.5)) +
  geom_vline(xintercept = 0, linetype = "dashed")



# Compare models with different fixed effects -----------------------------

# H1: vegetation
m.veg = sdmTMB_cv(formula = vert.mean.ses ~ canopy_height + veg_den + (1|ecoregion),
                  data = dat,
                  mesh = mesh300km,
                  spatial = "on",
                  reml = F, 
                  fold_ids = folds,
                  k_folds = length(unique(folds)))

m.clim =  sdmTMB_cv(formula = vert.mean.ses ~ tmax_warm + precip_dry + tmin_cold + precip_sea + precip_wet + (1|ecoregion),
                    data = dat,
                    mesh = mesh300km,
                    spatial = "on",
                    reml = F, 
                    fold_ids = folds,
                    k_folds = length(unique(folds)))

m.climvel = sdmTMB_cv(formula = vert.mean.ses ~ clim_velocity + (1|ecoregion),
                      data = dat,
                      mesh = mesh300km,
                      spatial = "on",
                      reml = F, 
                      fold_ids = folds,
                      k_folds = length(unique(folds)))

m.vegclim = sdmTMB_cv(formula = vert.mean.ses ~ canopy_height + veg_den + tmax_warm + precip_dry + tmin_cold + precip_sea + precip_wet + (1|ecoregion),
                      data = dat,
                      mesh = mesh300km,
                      spatial = "on",
                      reml = F, 
                      fold_ids = folds,
                      k_folds = length(unique(folds)))

m.vegclimvel = sdmTMB_cv(formula = vert.mean.ses ~ canopy_height + veg_den + clim_velocity + (1|ecoregion),
                         data = dat,
                         mesh = mesh300km,
                         spatial = "on",
                         reml = F, 
                         fold_ids = folds,
                         k_folds = length(unique(folds)))

m.climclimvel = sdmTMB_cv(formula = vert.mean.ses ~ tmax_warm + precip_dry + tmin_cold + precip_sea + precip_wet + clim_velocity + (1|ecoregion),
                          data = dat,
                          mesh = mesh300km,
                          spatial = "on",
                          reml = F, 
                          fold_ids = folds,
                          k_folds = length(unique(folds)))

m.all = sdmTMB_cv(formula = vert.mean.ses ~ canopy_height + veg_den + tmax_warm + precip_dry + tmin_cold + precip_sea + precip_wet + clim_velocity + (1|ecoregion),
                  data = dat,
                  mesh = mesh300km,
                  spatial = "on",
                  reml = F, 
                  fold_ids = folds,
                  k_folds = length(unique(folds)))

m.veg$sum_loglik
m.clim$sum_loglik
m.climvel$sum_loglik
m.vegclim$sum_loglik # vegclim performed the best by a small margin
m.vegclimvel$sum_loglik
m.climclimvel$sum_loglik
m.all$sum_loglik # best by a tiny margin

summ.veg = summ_coefs(m.veg, "m.veg")
summ.clim = summ_coefs(m.clim, "m.clim")
summ.climvel = summ_coefs(m.climvel, "m.climvel")
summ.vegclim = summ_coefs(m.vegclim, "m.vegclim")
summ.vegclimvel = summ_coefs(m.vegclimvel, "m.vegclimvel")
summ.climclimvel = summ_coefs(m.climclimvel, "m.climclimvel")
summ.all = summ_coefs(m.all, "m.all")

summ.m = rbind(summ.veg, summ.clim, summ.climvel,
               summ.vegclim, summ.vegclimvel, summ.climclimvel, summ.all)

ggplot(summ.m %>% filter(term != "(Intercept)"), aes(est, xmin = est-se*1.96, xmax = est+se*1.96, y = term, color = modname)) +
  geom_pointrange(position = position_dodge2(width = 0.5)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_brewer(palette = "Dark2")


# Refit best model using all the data ------------------------------------------
control = sdmTMBcontrol(nlminb_loops = 2)
m.final = sdmTMB(formula = vert.mean.ses ~ canopy_height + veg_den + tmax_warm + precip_dry + tmin_cold + precip_sea + precip_wet + (1|ecoregion),
                 data = dat,
                 mesh = mesh300km,
                 spatial = "on",
                 reml = T,
                 control = control)

m.final
sanity(m.final)

pred <- predict(m.final)
s_m.final = simulate(m.final, nsim = 500, seed = 12345, type = "mle-mvn")
r_m.final <- DHARMa::createDHARMa(
  simulatedResponse = s_m.final,
  observedResponse = dat$vert.mean.ses,
  fittedPredictedResponse = pred$est_non_rf
)

r_m.final <- DHARMa::createDHARMa(
  simulatedResponse = s_m.final,
  observedResponse = dat$vert.mean.ses,
  fittedPredictedResponse = pred$est
)

plot(r_m.final)
testQuantiles(r_m.final)
testResiduals(r_m.final)
plotResiduals(r_m.final, rank = F, form = dat$vert.mean.ses)
plot(dat$vert.mean.ses, pred$est_non_rf)
abline(a=0, b=1, col = "red")
plot(dat$vert.mean.ses, pred$est)
abline(a=0, b=1, col = "red")
testSpatialAutocorrelation(r_m.final, x = dat$x, y = dat$y)



# predict models to future env data ---------------------------------------
dat.f$ecoregion = factor(dat.f$ecoregion)
dat.f$x = dat.f$x/1e5
dat.f$y = dat.f$y/1e5

pred.f = predict(m.final, newdata = dat.f, type = "response")
head(pred.f)

pred.f$est.dif = pred.f$est - pred$est

ggplot() +
  geom_tile(data = pred.f, aes(x*1e5, y*1e5, fill = est.dif)) +
  scale_fill_continuous_divergingx("spectral") +
  coord_sf(crs = "+proj=cea +datum=WGS84")

save(m.final, pred, pred.f, file = "results/sdmTMB_models/squamates_sesvert.RData")


# MEAN VERTICALITY --------------------------------------------------------

head(dat)
ls = ls()
rm(list = ls[-c(2,3)])

# plot relationship between mean verticality and env predictors
dat %>% 
  dplyr::select(vert.mean, biome:clim_velocity, elev, veg_den, veg_complexity) %>%
  pivot_longer(cols = 3:16, names_to = "var", values_to = "val") %>% 
  ggplot(aes(x = val, y = vert.mean, color = biome)) +
  geom_point(pch = ".") +
  facet_wrap(~var, scales = "free") +
  theme_classic()


dat$ecoregion = factor(dat$ecoregion)
dat$realm = factor(dat$realm)
dat$biome = factor(dat$biome)
dat$biorealm = factor(paste(dat$biome, dat$realm, spe = "_"))
dat$x = dat$x/1e5
dat$y = dat$y/1e5

mesh300km = make_mesh(data = dat, xy_cols = c("x", "y"), cutoff = 3)
plot(mesh300km)

# set up five random folds for cross validation
folds = sample(1:5, size = nrow(dat), replace = T)

f1 = formula(vert.mean ~ canopy_height + veg_den + tmax_warm + tmin_cold + I(tmin_cold^2) + precip_sea + precip_wet + precip_dry + clim_velocity)

# 1: fit random intercept of ecoregion without spatial random field
library(future)
plan(multisession)

mod.ecoregion.cv = sdmTMB_cv(update(f1, ~ . + (1|ecoregion)), 
                             data = dat,
                             mesh = mesh300km,
                             spatial = "off",
                             family = Beta(),
                             reml = T, 
                             fold_ids = folds,
                             k_folds = length(unique(folds)))
mod.ecoregion.cv$converged
lapply(mod.ecoregion.cv$models, sanity)
# all okay

pred <- predict(mod.ecoregion.cv$models[[1]])
s_mod.ecoregion = simulate(mod.ecoregion.cv$models[[1]], nsim = 500, seed = 12345, type = "mle-mvn")
r_mod.ecoregion <- DHARMa::createDHARMa(
  simulatedResponse = s_mod.ecoregion,
  observedResponse = dat$vert.mean,
  fittedPredictedResponse = pred$est
)

plot(r_mod.ecoregion)
testQuantiles(r_mod.ecoregion)
testResiduals(r_mod.ecoregion)
plotResiduals(r_mod.ecoregion, rank = F, form = dat$vert.mean)
plot(dat$vert.mean, mod.ecoregion.cv$models[[1]]$family$linkinv(pred$est))
abline(a = 0, b = 1, col = "red")
testSpatialAutocorrelation(r_mod.ecoregion, x = dat$x, y = dat$y)


mod.mesh.cv = sdmTMB_cv(f1, 
                        data = dat,
                        mesh = mesh300km,
                        spatial = "on",
                        family = Beta(),
                        reml = T, 
                        fold_ids = folds,
                        k_folds = length(unique(folds)))

mod.mesh.cv$converged #true
lapply(mod.mesh.cv$models, sanity)
# gradient issues with 4, rest are okay


mod.mesh.biome.cv = sdmTMB_cv(update(f1, ~ . +(1|biome)), 
                              data = dat,
                              mesh = mesh300km,
                              spatial = "on",
                              family = Beta(),
                              reml = T, 
                              fold_ids = folds,
                              k_folds = length(unique(folds)))

mod.mesh.biome.cv$converged
lapply(mod.mesh.biome.cv$models, sanity)
# no models ok
# sigmas not okay in some models and gradients not okay in some models

mod.mesh.biorealm.cv = sdmTMB_cv(update(f1, ~ . +(1|biorealm)), 
                                 data = dat,
                                 mesh = mesh300km,
                                 spatial = "on",
                                 family = Beta(),
                                 reml = T, 
                                 fold_ids = folds,
                                 k_folds = length(unique(folds)))

mod.mesh.biorealm.cv$converged 
lapply(mod.mesh.biorealm.cv$models, sanity)
# issue with gradient in 1 - all others okay

mod.mesh.ecoregion.cv = sdmTMB_cv(update(f1, ~ . +(1|ecoregion)), 
                                  data = dat,
                                  mesh = mesh300km,
                                  spatial = "on",
                                  family = Beta(),
                                  reml = T, 
                                  fold_ids = folds,
                                  k_folds = length(unique(folds)))

mod.mesh.ecoregion.cv$converged #false
lapply(mod.mesh.ecoregion.cv$models, sanity)
# not okay models for 3,4,5

# compare sum_logLik
mod.ecoregion.cv$sum_loglik
mod.mesh.cv$sum_loglik
#mod.mesh.biome.cv$sum_loglik
mod.mesh.biorealm.cv$sum_loglik
mod.mesh.ecoregion.cv$sum_loglik #best

summ_coefs = function(x, nm){
  mods = x$models
  test = lapply(mods, sanity)
  test = sapply(test, "[[", 9) # get T/F for all_ok
  mods = mods[which(test)] # subset mods to okay models
  mods = lapply(mods, tidy)
  for (i in 1:length(mods)) {mods[[i]]$mod = as.character(i)}
  mods = rbindlist(mods)
  mods = mods %>% 
    group_by(term) %>% 
    summarise(est = mean(estimate), se = mean(std.error)) %>% 
    mutate(modname = nm)
  return(mods)
}

mod.ecoregion.summ = summ_coefs(mod.ecoregion.cv, "mod.ecoregion")
mod.mesh.summ = summ_coefs(mod.mesh.cv, "mod.mesh")
#mod.mesh.biome.summ = summ_coefs(mod.mesh.biome.cv, "mod.mesh.biome")
mod.mesh.biorealm.summ = summ_coefs(mod.mesh.biorealm.cv, "mod.mesh.biorealm")
mod.mesh.ecoregion.summ = summ_coefs(mod.mesh.ecoregion.cv, "mod.mesh.ecoregion")

mods = rbind(mod.ecoregion.summ, mod.mesh.summ, mod.mesh.biorealm.summ, mod.mesh.ecoregion.summ)
ggplot(mods %>% filter(term != "(Intercept)")) +
  geom_pointrange(aes(x = est, xmax = est+se*1.96, xmin = est-se*1.96, y = term, color = modname),
                  position = position_dodge2(width=0.5)) +
  geom_vline(xintercept = 0, linetype = "dashed")



# Compare models with different fixed effects -----------------------------

# H1: vegetation
m.veg = sdmTMB_cv(formula = vert.mean ~ canopy_height + veg_den + (1|ecoregion),
                  data = dat,
                  mesh = mesh300km,
                  spatial = "on",
                  family = Beta(),
                  reml = F, 
                  fold_ids = folds,
                  k_folds = length(unique(folds)))

m.clim =  sdmTMB_cv(formula = vert.mean ~ tmax_warm + precip_dry + I(tmin_cold^2) + tmin_cold + precip_sea + precip_wet + (1|ecoregion),
                    data = dat,
                    mesh = mesh300km,
                    spatial = "on",
                    family = Beta(),
                    reml = F, 
                    fold_ids = folds,
                    k_folds = length(unique(folds)))
lapply(m.clim$models, sanity)

m.clim2 =  sdmTMB_cv(formula = vert.mean ~ tmax_warm + precip_dry + tmin_cold + precip_sea + precip_wet + (1|ecoregion),
                    data = dat,
                    mesh = mesh300km,
                    spatial = "on",
                    family = Beta(),
                    reml = F, 
                    fold_ids = folds,
                    k_folds = length(unique(folds)))
lapply(m.clim2$models, sanity)

m.clim$sum_loglik
m.clim2$sum_loglik # better without quadratic term

m.climvel = sdmTMB_cv(formula = vert.mean ~ clim_velocity + (1|ecoregion),
                      data = dat,
                      mesh = mesh300km,
                      spatial = "on",
                      family = Beta(),
                      reml = F, 
                      fold_ids = folds,
                      k_folds = length(unique(folds)))
lapply(m.climvel$models, sanity)

m.vegclim = sdmTMB_cv(formula = vert.mean ~ canopy_height + veg_den + tmax_warm + precip_dry + tmin_cold + precip_sea + precip_wet + (1|ecoregion),
                      data = dat,
                      mesh = mesh300km,
                      spatial = "on",
                      family = Beta(),
                      reml = F, 
                      fold_ids = folds,
                      k_folds = length(unique(folds)))
lapply(m.vegclim$models, sanity)

m.vegclimvel = sdmTMB_cv(formula = vert.mean ~ canopy_height + veg_den + clim_velocity + (1|ecoregion),
                         data = dat,
                         mesh = mesh300km,
                         spatial = "on",
                         family = Beta(),
                         reml = F, 
                         fold_ids = folds,
                         k_folds = length(unique(folds)))

m.climclimvel = sdmTMB_cv(formula = vert.mean ~ tmax_warm + precip_dry + tmin_cold + precip_sea + precip_wet + clim_velocity + (1|ecoregion),
                          data = dat,
                          mesh = mesh300km,
                          spatial = "on",
                          family = Beta(),
                          reml = F, 
                          fold_ids = folds,
                          k_folds = length(unique(folds)))

m.all = sdmTMB_cv(formula = vert.mean ~ canopy_height + veg_den + tmax_warm + precip_dry + tmin_cold + precip_sea + precip_wet + clim_velocity + (1|ecoregion),
                  data = dat,
                  mesh = mesh300km,
                  spatial = "on",
                  family = Beta(),
                  reml = F, 
                  fold_ids = folds,
                  k_folds = length(unique(folds)))

m.veg$sum_loglik
m.clim$sum_loglik
m.clim2$sum_loglik
m.climvel$sum_loglik
m.vegclim$sum_loglik 
m.vegclimvel$sum_loglik
m.climclimvel$sum_loglik
m.all$sum_loglik # best model

summ.veg = summ_coefs(m.veg, "m.veg")
summ.clim = summ_coefs(m.clim, "m.clim")
summ.clim2 = summ_coefs(m.clim2, "m.clim2")
summ.climvel = summ_coefs(m.climvel, "m.climvel")
summ.vegclim = summ_coefs(m.vegclim, "m.vegclim")
summ.vegclimvel = summ_coefs(m.vegclimvel, "m.vegclimvel")
summ.climclimvel = summ_coefs(m.climclimvel, "m.climclimvel")
summ.all = summ_coefs(m.all, "m.all")

summ.m = rbind(summ.veg, summ.clim2, summ.climvel,
               summ.vegclim, summ.vegclimvel, summ.climclimvel, summ.all)

ggplot(summ.m %>% filter(term != "(Intercept)"), aes(est, xmin = est-se*1.96, xmax = est+se*1.96, y = term, color = modname)) +
  geom_pointrange(position = position_dodge2(width = 0.5)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_brewer(palette = "Dark2")


# Refit best model using all the data ------------------------------------------
control = sdmTMBcontrol(nlminb_loops = 2)
m.final = sdmTMB(formula = vert.mean ~ canopy_height + veg_den + tmax_warm + precip_dry + tmin_cold + precip_sea + precip_wet + clim_velocity + (1|ecoregion),
                 data = dat,
                 mesh = mesh300km,
                 family = Beta(),
                 spatial = "on",
                 reml = T,
                 control = control)

m.final
sanity(m.final)

pred <- predict(m.final)
s_m.final = simulate(m.final, nsim = 500, seed = 12345, type = "mle-mvn")
r_m.final <- DHARMa::createDHARMa(
  simulatedResponse = s_m.final,
  observedResponse = dat$vert.mean,
  fittedPredictedResponse = pred$est_non_rf
)

r_m.final <- DHARMa::createDHARMa(
  simulatedResponse = s_m.final,
  observedResponse = dat$vert.mean,
  fittedPredictedResponse = pred$est
)

plot(r_m.final)
testQuantiles(r_m.final)
testResiduals(r_m.final)
plotResiduals(r_m.final, rank = F, form = dat$vert.mean)
plot(dat$vert.mean, m.final$family$linkinv(pred$est_non_rf))
abline(a = 0, b = 1, col = "red")
plot(dat$vert.mean, m.final$family$linkinv(pred$est))
abline(a = 0, b = 1, col = "red")
testSpatialAutocorrelation(r_m.final, x = dat$x, y = dat$y)


dat.f$ecoregion = factor(dat.f$ecoregion)
dat.f$x = dat.f$x/1e5
dat.f$y = dat.f$y/1e5

pred.f = predict(m.final, newdata = dat.f, type = "response")
head(pred.f)

# predict in response scale
pred = predict(m.final, type = "response")

pred.f$est.dif = pred.f$est - pred$est

ggplot() +
  geom_tile(data = pred.f, aes(x*1e5, y*1e5, fill = est.dif)) +
  scale_fill_continuous_divergingx("spectral") +
  coord_sf(crs = "+proj=cea +datum=WGS84")

save(m.final, pred, pred.f, file = "results/sdmTMB_models/squamates_meanvert.RData")

