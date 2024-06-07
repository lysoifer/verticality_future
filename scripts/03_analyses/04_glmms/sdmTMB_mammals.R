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


mod.mesh.ecoregion.cv2 = sdmTMB_cv(update(f1, ~ . +(1|ecoregion)), 
                                  data = dat,
                                  mesh = mesh300km,
                                  spatial = "on",
                                  reml = T, 
                                  fold_ids = folds,
                                  k_folds = length(unique(folds)),
                                  control = sdmTMBcontrol(nlminb_loops = 2, start = list(ln_tau_o = 0.364, ln_phi = -0.988)))

mod.mesh.ecoregion.cv2$converged
lapply(mod.mesh.ecoregion.cv2$models, sanity)
# only model with issues still is #5

test1 = mod.mesh.ecoregion.cv2$models[[1]]$tmb_obj$env$parList()
test2 = mod.mesh.ecoregion.cv2$models[[2]]$tmb_obj$env$parList()
test3 = mod.mesh.ecoregion.cv2$models[[3]]$tmb_obj$env$parList()
test4 = mod.mesh.ecoregion.cv2$models[[4]]$tmb_obj$env$parList()

test1$ln_phi
test2$ln_phi
test3$ln_phi
test4$ln_phi

mod.mesh.ecoregion.cv3 = sdmTMB_cv(update(f1, ~ . +(1|ecoregion)), 
                                   data = dat,
                                   mesh = mesh300km,
                                   spatial = "on",
                                   reml = T, 
                                   fold_ids = folds,
                                   k_folds = length(unique(folds)),
                                   control = sdmTMBcontrol(
                                     nlminb_loops = 2, 
                                     start = list(ln_tau_o = 0.364, ln_phi = -0.989, b_j = c(-0.42, 0.08, -0.0005, -0.03, 0.086,0.19,0.03,0.13,-0.0333))))
                                   

mod.mesh.ecoregion.cv3$converged
lapply(mod.mesh.ecoregion.cv3$models, sanity)
# mod3 still not okay
# no differences in sum_log_lik, so go with the best one and then focus on fitting the full model

# compare sum_logLik
mod.ecoregion.cv$sum_loglik
mod.mesh.cv$sum_loglik
mod.mesh.biome.cv$sum_loglik
mod.mesh.biorealm.cv$sum_loglik
mod.mesh.ecoregion.cv3$sum_loglik # performed best

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
                     k_folds = length(unique(folds)))
m.nopoly$converged
lapply(m.nopoly$models, sanity)

m.noclimvel = sdmTMB_cv(formula = update(f1, ~ . -log_clim_velocity +(1|ecoregion)),
                        data = dat,
                        mesh = mesh300km,
                        spatial = "on",
                        reml = F, 
                        fold_ids = folds,
                        k_folds = length(unique(folds)))
m.noclimvel$converged
lapply(m.noclimvel$models, sanity)

m.full$sum_loglik # performs best
m.nopoly$sum_loglik
m.noclimvel$sum_loglik

# H1: vegetation
# m.veg = sdmTMB_cv(formula = vert.mean.ses ~ canopy_height + veg_den + (1|ecoregion),
#                   data = dat,
#                   mesh = mesh300km,
#                   spatial = "on",
#                   reml = F, 
#                   fold_ids = folds,
#                   k_folds = length(unique(folds)))
# 
# m.clim =  sdmTMB_cv(formula = vert.mean.ses ~ tmax_warm + precip_dry + tmin_cold + precip_sea + precip_wet + (1|ecoregion),
#                     data = dat,
#                     mesh = mesh300km,
#                     spatial = "on",
#                     reml = F, 
#                     fold_ids = folds,
#                     k_folds = length(unique(folds)))
# 
# m.climvel = sdmTMB_cv(formula = vert.mean.ses ~ clim_velocity + (1|ecoregion),
#                       data = dat,
#                       mesh = mesh300km,
#                       spatial = "on",
#                       reml = F, 
#                       fold_ids = folds,
#                       k_folds = length(unique(folds)))
# 
# m.vegclim = sdmTMB_cv(formula = vert.mean.ses ~ canopy_height + veg_den + tmax_warm + precip_dry + tmin_cold + precip_sea + precip_wet + (1|ecoregion),
#                       data = dat,
#                       mesh = mesh300km,
#                       spatial = "on",
#                       reml = F, 
#                       fold_ids = folds,
#                       k_folds = length(unique(folds)))
# 
# m.vegclimvel = sdmTMB_cv(formula = vert.mean.ses ~ canopy_height + veg_den + clim_velocity + (1|ecoregion),
#                          data = dat,
#                          mesh = mesh300km,
#                          spatial = "on",
#                          reml = F, 
#                          fold_ids = folds,
#                          k_folds = length(unique(folds)))
# 
# m.climclimvel = sdmTMB_cv(formula = vert.mean.ses ~ tmax_warm + precip_dry + tmin_cold + precip_sea + precip_wet + clim_velocity + (1|ecoregion),
#                           data = dat,
#                           mesh = mesh300km,
#                           spatial = "on",
#                           reml = F, 
#                           fold_ids = folds,
#                           k_folds = length(unique(folds)))
# 
# m.all = sdmTMB_cv(formula = vert.mean.ses ~ canopy_height + veg_den + tmax_warm + precip_dry + tmin_cold + precip_sea + precip_wet + clim_velocity + (1|ecoregion),
#                   data = dat,
#                   mesh = mesh300km,
#                   spatial = "on",
#                   reml = F, 
#                   fold_ids = folds,
#                   k_folds = length(unique(folds)))
# 
# m.veg$sum_loglik
# m.clim$sum_loglik
# m.climvel$sum_loglik
# m.vegclim$sum_loglik # vegclim performed the best by a small margin
# m.vegclimvel$sum_loglik
# m.climclimvel$sum_loglik
# m.all$sum_loglik # best by a tiny margin
# 
# summ.veg = summ_coefs(m.veg, "m.veg")
# summ.clim = summ_coefs(m.clim, "m.clim")
# summ.climvel = summ_coefs(m.climvel, "m.climvel")
# summ.vegclim = summ_coefs(m.vegclim, "m.vegclim")
# summ.vegclimvel = summ_coefs(m.vegclimvel, "m.vegclimvel")
# summ.climclimvel = summ_coefs(m.climclimvel, "m.climclimvel")
# summ.all = summ_coefs(m.all, "m.all")
# 
# summ.m = rbind(summ.veg, summ.clim, summ.climvel,
#                summ.vegclim, summ.vegclimvel, summ.climclimvel, summ.all)
# 
# ggplot(summ.m %>% filter(term != "(Intercept)"), aes(est, xmin = est-se*1.96, xmax = est+se*1.96, y = term, color = modname)) +
#   geom_pointrange(position = position_dodge2(width = 0.5)) +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   scale_color_brewer(palette = "Dark2")
# 

# Refit best model using all the data ------------------------------------------
control = sdmTMBcontrol(nlminb_loops = 2)
m.final = sdmTMB(formula = update(f1, . ~ . + (1|ecoregion)),
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


dat$ecoregion = factor(dat$ecoregion)
dat$realm = factor(dat$realm)
dat$biome = factor(dat$biome)
dat$biorealm = factor(paste(dat$biome, dat$realm, spe = "_"))
#dat$x = dat$x/1e5
#dat$y = dat$y/1e5

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
                        k_folds = length(unique(folds)))
m.noclimvel$converged
lapply(m.noclimvel$models, sanity)

m.full$sum_loglik # performs best
m.nopoly$sum_loglik
m.noclimvel$sum_loglik

# # H1: vegetation
# m.veg = sdmTMB_cv(formula = vert.mean ~ canopy_height + veg_den + (1|ecoregion),
#                   data = dat,
#                   mesh = mesh300km,
#                   spatial = "on",
#                   family = Beta(),
#                   reml = F, 
#                   fold_ids = folds,
#                   k_folds = length(unique(folds)))
# 
# m.clim =  sdmTMB_cv(formula = vert.mean ~ tmax_warm + precip_dry + I(tmin_cold^2) + tmin_cold + precip_sea + precip_wet + (1|ecoregion),
#                     data = dat,
#                     mesh = mesh300km,
#                     spatial = "on",
#                     family = Beta(),
#                     reml = F, 
#                     fold_ids = folds,
#                     k_folds = length(unique(folds)))
# lapply(m.clim$models, sanity)
# 
# m.clim2 =  sdmTMB_cv(formula = vert.mean ~ tmax_warm + precip_dry + tmin_cold + precip_sea + precip_wet + (1|ecoregion),
#                      data = dat,
#                      mesh = mesh300km,
#                      spatial = "on",
#                      family = Beta(),
#                      reml = F, 
#                      fold_ids = folds,
#                      k_folds = length(unique(folds)))
# lapply(m.clim2$models, sanity)
# 
# m.clim$sum_loglik
# m.clim2$sum_loglik # better without quadratic term
# 
# m.climvel = sdmTMB_cv(formula = vert.mean ~ clim_velocity + (1|ecoregion),
#                       data = dat,
#                       mesh = mesh300km,
#                       spatial = "on",
#                       family = Beta(),
#                       reml = F, 
#                       fold_ids = folds,
#                       k_folds = length(unique(folds)))
# lapply(m.climvel$models, sanity)
# 
# m.vegclim = sdmTMB_cv(formula = vert.mean ~ canopy_height + veg_den + tmax_warm + precip_dry + tmin_cold + precip_sea + precip_wet + (1|ecoregion),
#                       data = dat,
#                       mesh = mesh300km,
#                       spatial = "on",
#                       family = Beta(),
#                       reml = F, 
#                       fold_ids = folds,
#                       k_folds = length(unique(folds)))
# lapply(m.vegclim$models, sanity)
# 
# m.vegclimvel = sdmTMB_cv(formula = vert.mean ~ canopy_height + veg_den + clim_velocity + (1|ecoregion),
#                          data = dat,
#                          mesh = mesh300km,
#                          spatial = "on",
#                          family = Beta(),
#                          reml = F, 
#                          fold_ids = folds,
#                          k_folds = length(unique(folds)))
# 
# m.climclimvel = sdmTMB_cv(formula = vert.mean ~ tmax_warm + precip_dry + tmin_cold + precip_sea + precip_wet + clim_velocity + (1|ecoregion),
#                           data = dat,
#                           mesh = mesh300km,
#                           spatial = "on",
#                           family = Beta(),
#                           reml = F, 
#                           fold_ids = folds,
#                           k_folds = length(unique(folds)))
# 
# m.all = sdmTMB_cv(formula = vert.mean ~ canopy_height + veg_den + tmax_warm + precip_dry + tmin_cold + precip_sea + precip_wet + clim_velocity + (1|ecoregion),
#                   data = dat,
#                   mesh = mesh300km,
#                   spatial = "on",
#                   family = Beta(),
#                   reml = F, 
#                   fold_ids = folds,
#                   k_folds = length(unique(folds)))
# 
# m.veg$sum_loglik
# m.clim$sum_loglik
# m.clim2$sum_loglik
# m.climvel$sum_loglik
# m.vegclim$sum_loglik 
# m.vegclimvel$sum_loglik
# m.climclimvel$sum_loglik
# m.all$sum_loglik # best model
# 
# summ.veg = summ_coefs(m.veg, "m.veg")
# summ.clim = summ_coefs(m.clim, "m.clim")
# summ.clim2 = summ_coefs(m.clim2, "m.clim2")
# summ.climvel = summ_coefs(m.climvel, "m.climvel")
# summ.vegclim = summ_coefs(m.vegclim, "m.vegclim")
# summ.vegclimvel = summ_coefs(m.vegclimvel, "m.vegclimvel")
# summ.climclimvel = summ_coefs(m.climclimvel, "m.climclimvel")
# summ.all = summ_coefs(m.all, "m.all")
# 
# summ.m = rbind(summ.veg, summ.clim2, summ.climvel,
#                summ.vegclim, summ.vegclimvel, summ.climclimvel, summ.all)
# 
# ggplot(summ.m %>% filter(term != "(Intercept)"), aes(est, xmin = est-se*1.96, xmax = est+se*1.96, y = term, color = modname)) +
#   geom_pointrange(position = position_dodge2(width = 0.5)) +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   scale_color_brewer(palette = "Dark2")


# Refit best model using all the data ------------------------------------------
control = sdmTMBcontrol(nlminb_loops = 2)
m.final = sdmTMB(formula = update(f1, . ~ . + (1|ecoregion)),
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

