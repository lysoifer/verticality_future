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

# mammals forest only 50km resolution
raw = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/mammals_comdat.csv")

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


# plot relationship between SES verticality and env predictors
raw %>% 
  dplyr::select(vert.mean.ses, biome:clim_velocity, elev, veg_den, veg_complexity) %>%
  pivot_longer(cols = 3:16, names_to = "var", values_to = "val") %>% 
  ggplot(aes(x = val, y = vert.mean.ses, color = biome)) +
  geom_point(pch = ".") +
  facet_wrap(~var, scales = "free") +
  theme_classic()

dat.t = raw %>% 
  mutate(log_clim_velocity = log10(clim_velocity),
         log_tmin_cold = log10(tmin_cold + 50),
         log_precip_dry = log10(precip_dry),
         log_precip_wet = log10(precip_wet))

dat.t %>% 
  dplyr::select(vert.mean.ses, biome:clim_velocity, elev, veg_den, veg_complexity, log_clim_velocity, log_precip_dry, log_tmin_cold, log_precip_wet) %>%
  pivot_longer(cols = 3:20, names_to = "var", values_to = "val") %>% 
  ggplot(aes(x = val, y = vert.mean.ses, color = biome)) +
  geom_point(pch = ".") +
  facet_wrap(~var, scales = "free") +
  theme_classic()

# scale the data
dat = dat.t %>%
  mutate_at(.vars = vars(canopy_height:tmin_cold, elev, veg_den, veg_complexity, log_clim_velocity, log_precip_dry), 
            .funs = function(x) scale(x)[,1]) %>% 
  dplyr::select(!precip_dry)

save(dat, file = "scripts/03_analyses/00_testing_ideas/test_mesh_dat_mammals.RData")


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

vif(lm(vert.mean.ses ~ canopy_height + veg_den + veg_complexity + tmax_warm + tmin_cold + temp_sea + precip_sea + precip_wet + log_precip_dry + log_clim_velocity, data = dat))

# remove veg_complexity
vif(lm(vert.mean.ses ~ canopy_height + veg_den + tmax_warm + tmin_cold + temp_sea + precip_sea + precip_wet + log_precip_dry + log_clim_velocity, data = dat))

# remove temp_sea
vif(lm(vert.mean.ses ~ canopy_height + veg_den + tmax_warm + tmin_cold  + precip_sea + precip_wet + log_precip_dry + log_clim_velocity, data = dat))

# remove temp_sea
vif(lm(vert.mean.ses ~ canopy_height + veg_den + tmax_warm + tmin_cold  + precip_wet + log_precip_dry + log_clim_velocity, data = dat))


# VIF all under 4

f1 = formula(vert.mean.ses ~ canopy_height + veg_den + I(tmax_warm^2) + tmax_warm + tmin_cold + precip_wet + log_precip_dry + log_clim_velocity)


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

mesh400km = make_mesh(data = dat, xy_cols = c("x", "y"), cutoff = 4)
plot(mesh400km)

mesh100km = make_mesh(data = dat, xy_cols = c("x", "y"), cutoff = 1)
plot(mesh100km)


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
                              parallel = T)

mod.mesh.biome.cv$converged
lapply(mod.mesh.biome.cv$models, sanity)
# issues with gradient in 2,3

mod.mesh.biorealm.cv = sdmTMB_cv(update(f1, ~ . +(1|biorealm)), 
                                 data = dat,
                                 mesh = mesh300km,
                                 spatial = "on",
                                 reml = T, 
                                 fold_ids = folds,
                                 k_folds = length(unique(folds)),
                                 parallel = T)

mod.mesh.biorealm.cv$converged
lapply(mod.mesh.biorealm.cv$models, sanity)
# all okay

mod.mesh.ecoregion.cv = sdmTMB_cv(update(f1, ~ . +(1|ecoregion)), 
                                  data = dat,
                                  mesh = mesh300km,
                                  spatial = "on",
                                  reml = T, 
                                  fold_ids = folds,
                                  k_folds = length(unique(folds)),
                                  control = sdmTMBcontrol(nlminb_loops = 2))

mod.mesh.ecoregion.cv$converged
lapply(mod.mesh.ecoregion.cv$models, sanity)
# issues with gradient in  1,3,4

test1 = mod.mesh.ecoregion.cv$models[[2]]$tmb_obj$env$parList()
test2 = mod.mesh.ecoregion.cv$models[[5]]$tmb_obj$env$parList()
test1$ln_tau_O
test2$ln_tau_O
test1$ln_phi
test2$ln_phi


# compare sum_logLik
mod.ecoregion.cv$sum_loglik
mod.mesh.cv$sum_loglik
mod.mesh.biome.cv$sum_loglik
mod.mesh.biorealm.cv$sum_loglik
mod.mesh.ecoregion.cv3$sum_loglik # performed best

mod.mesh.spatialvar.cv = sdmTMB_cv(update(f1, ~ . +(1|ecoregion)), 
                                  data = dat,
                                  mesh = mesh300km,
                                  spatial = "on",
                                  reml = T, 
                                  fold_ids = folds,
                                  k_folds = length(unique(folds)),
                                  spatial_varying = ~ 0 + canopy_height,
                                  control = sdmTMBcontrol(nlminb_loops = 2))

mod.mesh.spatialvar.cv$converged
lapply(mod.mesh.spatialvar.cv$models, sanity)

mod.mesh.ecoregion.cv$sum_loglik
mod.mesh.spatialvar.cv$sum_loglik # best


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
mod.mesh.ecoregion.summ = summ_coefs(mod.mesh.ecoregion.cv3, "mod.mesh.ecoregion")
mod.mesh.spatialvar.summ = summ_coefs(mod.mesh.spatialvar, "mod.mesh.spatialvar")

mods = rbind(mod.ecoregion.summ, mod.mesh.summ, mod.mesh.biome.summ, mod.mesh.biorealm.summ, mod.mesh.ecoregion.summ)
ggplot(mods %>% filter(term != "(Intercept)")) +
  geom_pointrange(aes(x = est, xmax = est+se, xmin = est-se, y = term, color = modname),
                  position = position_dodge2(width=0.5)) +
  geom_vline(xintercept = 0, linetype = "dashed")



# Compare models with different fixed effects -----------------------------

m.full = sdmTMB_cv(formula = update(f1, ~ . +(1|ecoregion)),
                   data = dat,
                   mesh = mesh300km,
                   spatial = "on",
                   reml = F, 
                   fold_ids = folds,
                   spatial_varying = ~ 0 + canopy_height,
                   k_folds = length(unique(folds)))
m.full$converged
lapply(m.full$models, sanity)
# all okay


m.nopoly = sdmTMB_cv(formula = update(f1, ~ . -I(tmax_warm^2) +(1|ecoregion)),
                     data = dat,
                     mesh = mesh300km,
                     spatial = "on",
                     reml = F, 
                     fold_ids = folds,
                     spatial_varying = ~ 0 + canopy_height,
                     k_folds = length(unique(folds)))
m.nopoly$converged
lapply(m.nopoly$models, sanity)

m.noclimvel = sdmTMB_cv(formula = update(f1, ~ . -log_clim_velocity +(1|ecoregion)),
                        data = dat,
                        mesh = mesh300km,
                        spatial = "on",
                        reml = F, 
                        fold_ids = folds,
                        spatial_varying = ~ 0 + canopy_height,
                        k_folds = length(unique(folds)))
m.noclimvel$converged
lapply(m.noclimvel$models, sanity)

m.noclimvel_nopoly = sdmTMB_cv(formula = update(f1, ~ . -I(tmax_warm^2) -log_clim_velocity +(1|ecoregion)),
                        data = dat,
                        mesh = mesh300km,
                        spatial = "on",
                        reml = F, 
                        fold_ids = folds,
                        spatial_varying = ~ 0 + canopy_height,
                        k_folds = length(unique(folds)))
m.noclimvel_nopoly$converged
lapply(m.noclimvel_nopoly$models, sanity)

m.full$sum_loglik 
m.nopoly$sum_loglik
m.noclimvel$sum_loglik 
m.noclimvel_nopoly$sum_loglik # performs best


# Refit best model using all the data ------------------------------------------
control = sdmTMBcontrol(nlminb_loops = 2)
m.final = sdmTMB(formula = update(f1, . ~ . -I(tmax_warm) - log_clim_velocity + (1|ecoregion)),
                 data = dat,
                 mesh = mesh300km,
                 spatial = "on",
                 reml = T,
                 spatial_varying = ~ 0 + canopy_height,
                 control = control)

m.final
sanity(m.final)


m.final.400 = sdmTMB(formula = update(f1, . ~ . -I(tmax_warm) - log_clim_velocity + (1|ecoregion)),
                 data = dat,
                 mesh = mesh400km,
                 spatial = "on",
                 reml = T,
                 spatial_varying = ~ 0 + canopy_height,
                 control = sdmTMBcontrol(nlminb_loops = 2, start = list(ln_phi = 0.07)))

m.final.400
sanity(m.final.400)

pred.400 = predict(m.final.400)
ggplot(pred.400, aes(x = x*1e5, y = y*1e5, fill = zeta_s_canopy_height)) + 
  geom_tile() +
  scale_fill_continuous_divergingx("BrBG", na.value = "gray") +
  coord_sf(crs = "+proj=cea + datum=WGS84") +
  theme_classic()


m.final.100 = sdmTMB(formula = update(f1, . ~ . -I(tmax_warm) - log_clim_velocity + (1|ecoregion)),
                     data = dat,
                     mesh = mesh100km,
                     spatial = "on",
                     reml = T,
                     spatial_varying = ~ 0 + canopy_height,
                     control = sdmTMBcontrol(nlminb_loops = 2, start = list(ln_phi = -1)))
sdmTMB::run_extra_optimization(m.final.100)
m.final.100
sanity(m.final.100)

pred.100 = predict(m.final.100)
ggplot(pred.100, aes(x = x*1e5, y = y*1e5, fill = zeta_s_canopy_height)) + 
  geom_tile() +
  scale_fill_continuous_divergingx("BrBG", na.value = "gray") +
  coord_sf(crs = "+proj=cea + datum=WGS84") +
  theme_classic()



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

# need to subset predictions to calculate spatial autocorrelation - otherwise there are too many points
sub = sample(1:nrow(s_m.final), 10000)
r_m.final.sub = DHARMa::createDHARMa(
  simulatedResponse = s_m.final[sub,],
  observedResponse = dat$vert.mean.ses[sub],
  fittedPredictedResponse = pred$est_non_rf[sub]
)
testSpatialAutocorrelation(r_m.final.sub, x = dat$x[sub], y = dat$y[sub])



# predict models to future env data ---------------------------------------
dat.f$ecoregion = factor(dat.f$ecoregion)
dat.f$x = dat.f$x/1e5
dat.f$y = dat.f$y/1e5

pred.f = predict(m.final, newdata = dat.f, type = "response")
head(pred.f)

pred = predict(m.final, type = "response")

pred.f$est.dif = pred.f$est - pred$est

ggplot() +
  geom_tile(data = pred.f, aes(x*1e5, y*1e5, fill = est.dif)) +
  scale_fill_continuous_divergingx("spectral") +
  coord_sf(crs = "+proj=cea +datum=WGS84")

save(m.final, pred, pred.f, file = "results/sdmTMB_models/mammals_sesvert.RData")


# MEAN VERTICALITY --------------------------------------------------------

head(dat)
ls = ls()
rm(list = ls[-c(2,3)])

# plot relationship between mean verticality and env predictors
dat %>% 
  dplyr::select(vert.mean, biome:clim_velocity, elev, veg_den, veg_complexity, log_precip_dry, log_clim_velocity) %>%
  pivot_longer(cols = 3:17, names_to = "var", values_to = "val") %>% 
  ggplot(aes(x = val, y = vert.mean, color = biome)) +
  geom_point(pch = ".") +
  facet_wrap(~var, scales = "free") +
  theme_classic()



mesh300km = make_mesh(data = dat, xy_cols = c("x", "y"), cutoff = 3)
plot(mesh300km)

# set up five random folds for cross validation
folds = sample(1:5, size = nrow(dat), replace = T)

f1 = formula(vert.mean ~ canopy_height + veg_den + I(tmax_warm^2) +tmax_warm + tmin_cold + precip_wet + log_precip_dry + log_clim_velocity)

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
# gradient issues with 1,3,4

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
# issues with gradients in 1 and 2


mod.mesh.biome.cv = sdmTMB_cv(update(f1, ~ . +(1|biome)), 
                              data = dat,
                              mesh = mesh300km,
                              spatial = "on",
                              family = Beta(),
                              reml = T, 
                              fold_ids = folds,
                              k_folds = length(unique(folds)))

mod.mesh.biome.cv$converged
sanity(mod.mesh.biome.cv$models[[1]]) #ok
sanity(mod.mesh.biome.cv$models[[2]]) #ok
sanity(mod.mesh.biome.cv$models[[3]]) #ok
sanity(mod.mesh.biome.cv$models[[4]]) #ok
sanity(mod.mesh.biome.cv$models[[5]]) # gradient issues

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
# 1,2,5 have gradient issues

mod.mesh.ecoregion.cv = sdmTMB_cv(update(f1, ~ . +(1|ecoregion)), 
                                  data = dat,
                                  mesh = mesh300km,
                                  spatial = "on",
                                  family = Beta(),
                                  reml = T, 
                                  fold_ids = folds,
                                  k_folds = length(unique(folds)))

mod.mesh.ecoregion.cv$converged
lapply(mod.mesh.ecoregion.cv$models, sanity)
# 3 has gradient issues

# compare sum_logLik
mod.ecoregion.cv$sum_loglik
mod.mesh.cv$sum_loglik
mod.mesh.biome.cv$sum_loglik
mod.mesh.biorealm.cv$sum_loglik
mod.mesh.ecoregion.cv$sum_loglik #best

mod.mesh.spatialvar.cv = sdmTMB_cv(update(f1, ~ . +(1|ecoregion)), 
                                  data = dat,
                                  mesh = mesh300km,
                                  spatial = "on",
                                  family = Beta(),
                                  reml = T, 
                                  fold_ids = folds,
                                  spatial_varying = ~0 + canopy_height,
                                  k_folds = length(unique(folds)))

mod.mesh.spatialvar.cv$converged
lapply(mod.mesh.spatialvar.cv$models, sanity)

mod.mesh.ecoregion.cv$sum_loglik
mod.mesh.spatialvar.cv$sum_loglik

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
mod.mesh.spatialvar.summ = summ_coefs(mod.mesh.spatialvar.cv, "mod.mesh.spatialvar")

mods = rbind(mod.ecoregion.summ, mod.mesh.summ, mod.mesh.biome.summ, mod.mesh.biorealm.summ, mod.mesh.ecoregion.summ, mod.mesh.spatialvar.summ)
ggplot(mods %>% filter(term != "(Intercept)")) +
  geom_pointrange(aes(x = est, xmax = est+se*1.96, xmin = est-se*1.96, y = term, color = modname),
                  position = position_dodge2(width=0.5)) +
  geom_vline(xintercept = 0, linetype = "dashed")



# Compare models with different fixed effects -----------------------------

m.full = sdmTMB_cv(formula = update(f1, ~ . +(1|ecoregion)),
                   data = dat,
                   mesh = mesh300km,
                   spatial = "on",
                   family = Beta(),
                   reml = F, 
                   fold_ids = folds,
                   spatial_varying = ~0 + canopy_height,
                   k_folds = length(unique(folds)))
m.full$converged
lapply(m.full$models, sanity)
# all okay


m.nopoly = sdmTMB_cv(formula = update(f1, ~ . -I(tmax_warm^2) +(1|ecoregion)),
                     data = dat,
                     mesh = mesh300km,
                     spatial = "on",
                     family = Beta(),
                     reml = F, 
                     fold_ids = folds,
                     spatial_varying = ~ 0 + canopy_height,
                     k_folds = length(unique(folds)))
m.nopoly$converged
lapply(m.nopoly$models, sanity)
# all okay



m.noclimvel = sdmTMB_cv(formula = update(f1, ~ . -log_clim_velocity +(1|ecoregion)),
                        data = dat,
                        mesh = mesh300km,
                        spatial = "on",
                        family = Beta(),
                        reml = F, 
                        fold_ids = folds,
                        spatial_varying = ~ 0 + canopy_height,
                        k_folds = length(unique(folds)))
m.noclimvel$converged
lapply(m.noclimvel$models, sanity)

m.noclimvel_nopoly = sdmTMB_cv(formula = update(f1, ~ . -I(tmax_warm^2) -log_clim_velocity +(1|ecoregion)),
                        data = dat,
                        mesh = mesh300km,
                        spatial = "on",
                        family = Beta(),
                        reml = F, 
                        fold_ids = folds,
                        spatial_varying = ~ 0 + canopy_height,
                        k_folds = length(unique(folds)))
m.noclimvel_nopoly$converged
lapply(m.noclimvel_nopoly$models, sanity)

m.full$sum_loglik 
m.nopoly$sum_loglik
m.noclimvel$sum_loglik
m.noclimvel_nopoly$sum_loglik # performs best

test1 = m.noclimvel_nopoly$models[[1]]$tmb_obj$env$parList()
test3 = m.noclimvel_nopoly$models[[3]]$tmb_obj$env$parList()
test4 = m.noclimvel_nopoly$models[[4]]$tmb_obj$env$parList()
test5 = m.noclimvel_nopoly$models[[5]]$tmb_obj$env$parList()
test1$ln_phi
test3$ln_phi
test4$ln_phi
test5$ln_phi



# Refit best model using all the data ------------------------------------------
control = sdmTMBcontrol(nlminb_loops = 2)
m.final = sdmTMB(formula = update(f1, . ~ . -I(tmax_warm^2) -log_clim_velocity + (1|ecoregion)),
                 data = dat,
                 mesh = mesh300km,
                 family = Beta(),
                 spatial = "on",
                 spatial_varying = ~0 + canopy_height,
                 reml = T,
                 control = sdmTMBcontrol(start = list(ln_phi = 7.62), nlminb_loops = 2, iter.max = 2000, eval.max = 2000))

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
#testSpatialAutocorrelation(r_m.final, x = dat$x, y = dat$y)

# need to subset predictions to calculate spatial autocorrelation - otherwise there are too many points
sub = sample(1:nrow(s_m.final), 10000)
r_m.final.sub = DHARMa::createDHARMa(
  simulatedResponse = s_m.final[sub,],
  observedResponse = dat$vert.mean[sub],
  fittedPredictedResponse = pred$est_non_rf[sub]
)
testSpatialAutocorrelation(r_m.final.sub, x = dat$x[sub], y = dat$y[sub])


pred.f = predict(m.final, newdata = dat.f, type = "response")
head(pred.f)

# predict in response scale
pred = predict(m.final, type = "response")

pred.f$est.dif = pred.f$est - pred$est

ggplot() +
  geom_tile(data = pred.f, aes(x*1e5, y*1e5, fill = est.dif)) +
  scale_fill_continuous_divergingx("spectral") +
  coord_sf(crs = "+proj=cea +datum=WGS84")

save(m.final, pred, pred.f, file = "results/sdmTMB_models/mammals_meanvert.RData")

