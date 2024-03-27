source("scripts/00_functions/00_functions.R")
source("scripts/00_functions/00_functions_predict_sarlm_debug.R")
library(tidyverse)
library(vegan)
library(pgirmess)
library(ncf)
library(spdep)
library(sp)
library(DHARMa)
library(spatialreg)
library(MuMIn)
library(data.table)
library(parallel)
library(foreach)
library(doParallel)
library(terra)
library(colorspace)
library(patchwork)
library(broom)

options(na.action = "na.omit")

# VIF variable selection
env = read.csv("data/derivative_data/env_data.csv")

#env_vars = vif_func(in_frame = env[,c(4:18,22)], thresh = 3, trace = T)

# variable selection without hurs and vpd
#env_vars = vif_func(in_frame = env[c(4:12,17:18,22)], thresh = 3, trace = T)

env_vars = vif_func(in_frame = env[c(4, 7:12,17:18,22)], thresh = 5, trace = T)



# Birds elton --------------------------------------------------------------

birds = read.csv("data/derivative_data/gridcell_data/birds_comdat/birds_breedresident_comdat_elton.csv")

#  * - data check --------------------------------------------------------------

## SES vert model
sesvert = birds[,c("x", "y", "vert.mean.ses", "biome", "realm", env_vars)] %>% 
  drop_na()

hist(sesvert$vert.mean.ses)

# check linear association between predictor and response
as.data.frame(sesvert) %>% 
  pivot_longer(cols = 6:ncol(sesvert), names_to = "vars", values_to = "val") %>% 
  ggplot(aes(x = val, y = vert.mean.ses)) +
  geom_point(col = "grey", alpha = 0.2, size = 1) +
  geom_smooth(col = "blue", method = "lm") +
  facet_wrap(facets = vars(vars), scales = "free_x") +
  theme_classic()

# transform clim_velocity and temp seasonality and elevation
# sesvert.t = as.data.frame(sesvert) %>% 
#   mutate(log_clim_velocity = log10(clim_velocity),
#          log_precp_dry = log10(precip_dry),
#          log_temp_sea = log10(temp_sea))

sesvert.t = as.data.frame(sesvert) %>% 
  mutate(log_clim_velocity = log10(clim_velocity),
         log_precp_dry = log10(precip_dry),
         log_elev = log10(elev),
         log_precip_wet = log10(precip_wet),
         log_precip_sea = log10(precip_sea))


sesvert.t %>% 
  pivot_longer(cols = 6:ncol(sesvert.t), names_to = "vars", values_to = "val") %>% 
  ggplot(aes(x = val, y = vert.mean.ses)) +
  geom_point(col = "grey", alpha = 0.2, size = 1) +
  geom_smooth(col = "blue", method = "lm") +
  facet_wrap(facets = vars(vars), scales = "free_x") +
  theme_classic()

sesvert.t = sesvert.t %>% 
  dplyr::select(!c(clim_velocity, elev, precip_sea, log_precp_dry, log_precip_wet))


# scale predictor variables
sesvert.scale = sesvert.t %>% 
  mutate_at(.vars = 6:ncol(sesvert.t), .funs = function(x){scale(x)})
colnames(sesvert.scale) = colnames(sesvert.t)

# SES model
f = formula(vert.mean.ses ~ log_elev + log_clim_velocity + tmax_warm + tmin_cold + precip_dry + log_precip_sea + precip_wet + veg_complexity)
#f.realm = formula(vert.mean.ses ~ 0 + realm/(elev + log_clim_velocity + precip_dry + precip_sea + precip_wet + temp_diu + log_temp_sea + veg_complexity))

sesvert.mod = glm(f , data = sesvert.scale)
plot(sesvert.mod)
hist(sesvert.mod$residuals)
summary(sesvert.mod)

# sesvert.glm.realm = glm(f.realm, data = sesvert.scale)
# summary(sesvert.glm.realm)

# add residuals to dataframe to test for spatial autocorrelation
sesvert.scale$resid = sesvert.mod$residuals


# * - Spatial autocorrelation ---------------------------------------------
# see https://geo200cn.github.io/spatialreg.html (very helpful tutorial)

# examine spatial autocorrelation in response var with a correlogram
## take sample to look at spatial autocorrelation
# samp = sesvert.scale %>% sample_n(1000)
# sesvert.cor = ncf::correlog(x = samp$x, y = samp$y, z = samp$vert.mean.ses, increment = 111000, resamp = 99)
# ncf:::plot.correlog(sesvert.cor)
# ncf:::plot.correlog(sesvert.cor, xlim = c(0,5e6))
# 
# # make connectivity matrix (i.e., neighborhood matrix)
# 
# ## distance-based nearest neighbors (1st degree queen's case)
# neigh.160km = dnearneigh(x = sesvert.scale[,c("x", "y")], d1 = 0, d2 = 160000, longlat = F)
# 
# ## distance-based nearest neighbors (2nd degree queen's case)
# neigh.420km = dnearneigh(x = sesvert.scale[,c("x", "y")], d1 = 0, d2 = 420000, longlat = F)
# 
# ## distance-based nearest neighbors (up to 3000km - based on inspection of the correlogram)
# ## spatial autocorrelation is apparent until approximately 3000km (3million meters)
# neigh.3000km = dnearneigh(x = sesvert.scale[,c("x", "y")], d1 = 0, d2 = 3000000, longlat = F)
# 
# # make neighborhood weights matrix
# ## row standardized
# wts.160km = nb2listw(neighbours = neigh.160km, style = "W", zero.policy = T)
# wts.420km = nb2listw(neighbours = neigh.420km, style = "W", zero.policy = T)
# wts.3000km = nb2listw(neighbours = neigh.3000km, style = "W", zero.policy = T)
# 
# ## distance-based and row standardized weights matrix
# ### get distance between focal cell and neighbors within 2000km
# ### use 3000km because this is the distance at which spautocor appears to disappear based on correlogram
# dists.3000km = nbdists(nb = neigh.3000km, coords = sesvert.scale[,c("x", "y")], longlat = F)
# 
# ### inverse distances
# idists.3000km = lapply(dists.3000km, FUN = function(x){1/x})
# 
# ### inverse squared distances
# idists2.3000km = lapply(dists.3000km, FUN = function(x){1/(x^2)})
# 
# ### specify distances when creating weights matrix
# wts.idists.3000km = nb2listw(neighbours = neigh.3000km, glist = idists.3000km, style = "W", zero.policy = T)
# wts.idists2.3000km = nb2listw(neighbours = neigh.3000km, glist = idists2.3000km, style = "W", zero.policy = T)
# 
# # plot spatial autocorrelation in response var
# sesvert.scale.r = rast(sesvert.scale, type = "xyz", crs = "+proj=cea +datum=WGS84")
# plot(sesvert.scale.r$vert.mean.ses)
# 
# ## local spatial autocorrelation of ses.vert
# moran.plot(sesvert.scale$vert.mean.ses, listw = wts.160km)
# moran.plot(sesvert.scale$vert.mean.ses, listw = wts.420km)
# moran.plot(sesvert.scale$vert.mean.ses, listw = wts.3000km)
# moran.plot(sesvert.scale$vert.mean.ses, listw = wts.idists.3000km)
# moran.plot(sesvert.scale$vert.mean.ses, listw = wts.idists2.3000km)
# 
# ## local moran plots show that there is a positive correlation between vert.mean.ses and the average vert.mean.ses among neighbors
# ## i.e., there is clear spatial autocorrelation that must be accounted for in the models
# 
# 
# # plot spatial autocorrelation in residuals of glm
# plot(sesvert.scale.r$resid)
# 
# ## plot local spatial autocorrelatio in resids
# moran.plot(sesvert.scale$resid, listw = wts.160km)
# moran.plot(sesvert.scale$resid, listw = wts.420km)
# moran.plot(sesvert.scale$resid, listw = wts.3000km)
# moran.plot(sesvert.scale$resid, listw = wts.idists.3000km)
# moran.plot(sesvert.scale$resid, listw = wts.idists2.3000km)
# 
# ## clear spatial autocorrelation in residuals
# 
# # formally test spatial autocorrelation in residuals of glm using Moran's I
# moran.160km = moran.mc(sesvert.scale$resid, listw = wts.160km, nsim = 99, zero.policy = T)
# moran.160km
# 
# moran.420km = moran.mc(sesvert.scale$resid, listw = wts.420km, nsim = 99, zero.policy = T)
# moran.420km
# 
# moran.3000km = moran.mc(sesvert.scale$resid, listw = wts.3000km, nsim = 99, zero.policy = T)
# moran.3000km
# 
# moran.idists.3000km = moran.mc(sesvert.scale$resid, listw = wts.idists.3000km, nsim = 99, zero.policy = T)
# moran.idists.3000km
# 
# moran.idists2.3000km = moran.mc(sesvert.scale$resid, listw = wts.idists2.3000km, nsim = 99, zero.policy = T)
# moran.idists2.3000km
# ## make neighborhood weights matrix
# 
# ## test spatial autocorrelation in residuals using lm.morantest()
# lm.morantest(sesvert.mod, listw = wts.idists.3000km)


# * - selection of weights matrix ---------------------------------------------
## 1) fit full models with each weight matrix

### test neighborhood thresholds between 157km (queen's case 1st degree neighbors) and 
### 3000000m (queen's case 10th degree neighbors
### (expect spatial neighborhood to not exceed 2000km because of spautocor in response var)
### https://www.nature.com/articles/nature09329#Sec5
### simultaneously compare inverse distance and inverse distance^2 weighting
### https://www.sciencedirect.com/science/article/pii/S0006320717316890#s0010

dists = seq(157000, 3000000, 157000)
dists = dists[1:7] # try just dists 1-7 for now to speed up computation since dist 1 was the only one resulting in no more spatial autocorrelation in amphibians and mammals
sesvert.scale.sf = st_as_sf(sesvert.scale, coords = c("x","y"), crs = "+proj=cea +datum=WGS84")

cl = makeCluster(7)
registerDoParallel(cl)

# row standardized (SAR model with row standardized weights (Ver Hoef et al. 2018))
wts.rs <- foreach(i=1:length(dists), .packages = c("spdep")) %dopar% {
  n = dnearneigh(sesvert.scale[,c("x", "y")], d1 = 0, d2 = dists[i], longlat = F)
  nb2listw(neighbours = n, style = "W", zero.policy = T)
}

stopCluster(cl)

# inverse distance row standardized
cl = makeCluster(7)
registerDoParallel(cl)

# row standardized (SAR model with row standardized weights (Ver Hoef et al. 2018))
wts.idist <- foreach(i=1:length(dists), .packages = c("spdep")) %dopar% {
  n = dnearneigh(sesvert.scale[,c("x", "y")], d1 = 0, d2 = dists[i], longlat = F)
  nb2listwdist(neighbours = n, x = sesvert.scale.sf, style = "W", type = "idw", alpha = 1, zero.policy = T)
}

stopCluster(cl)


# inverse distance squared row standardized
cl = makeCluster(7)
registerDoParallel(cl)

wts.idist2 <- foreach(i=1:length(dists), .packages = c("spdep")) %dopar% {
  n = dnearneigh(sesvert.scale[,c("x", "y")], d1 = 0, d2 = dists[i], longlat = F)
  nb2listwdist(neighbours = n, x = sesvert.scale.sf, style = "W", type = "idw", alpha = 2, zero.policy = T)
}

stopCluster(cl)

### make lists to store models
### one list for each weight type

# row standardized weights
cl = makeCluster(7)
registerDoParallel(cl)

# row standardized (SAR model with row standardized weights (Ver Hoef et al. 2018))
mod.rs <- foreach(i=1:length(dists), .packages = c("spatialreg")) %dopar% {
  errorsarlm(f, data = sesvert.scale, listw = wts.rs[[i]], method = "LU", zero.policy = T, control = list(returnHcov=FALSE))
}

stopCluster(cl)

#save(mod.rs, file = "results/birdselton_sar_mod_rs.RData")
#load("results/birdselton_sar_mod_rs.RData")

# returnHcov=F; https://stat.ethz.ch/pipermail/r-sig-geo/2013-October/019450.html
#testmod = errorsarlm(f, data = sesvert.scale, listw = wts.rs[[1]], method = "LU", zero.policy = T, control = list(returnHcov=FALSE))
#testmod2 = spautolm(f, data = sesvert.scale, listw = wts.rs[[1]], method = "LU", family = "SAR", verbose = T)

# inverse distance row standardized weights
cl = makeCluster(7)
registerDoParallel(cl)

# row standardized (SAR model with row standardized weights (Ver Hoef et al. 2018))
mod.idist <- foreach(i=1:length(dists), .packages = c("spatialreg")) %dopar% {
  errorsarlm(f, data = sesvert.scale, listw = wts.idist[[i]], method = "LU", zero.policy = T, control = list(returnHcov=FALSE))
}

stopCluster(cl)

#save(mod.idist, file = "results/birdselton_sar_mod_idist.RData")


# inverse distance squared row standardized weights
cl = makeCluster(7)
registerDoParallel(cl)

# row standardized (SAR model with row standardized weights (Ver Hoef et al. 2018))
mod.idist2 <- foreach(i=1:length(dists), .packages = c("spatialreg")) %dopar% {
  errorsarlm(f, data = sesvert.scale, listw = wts.idist2[[i]], method = "LU", zero.policy = T, control = list(returnHCov=FALSE))
}

stopCluster(cl)

#save(mod.idist2, file = "results/birdselton_sar_mod_idist2.RData")


## 2) compare full models using AIC

### make lists to store models
### one list for each weight type

names(mod.rs) = paste0("mod.rs_", 1:length(mod.rs))
names(mod.idist) = paste0("mod.idist_", 1:length(mod.idist))
names(mod.idist2) = paste0("mod.idist2_", 1:length(mod.idist2))

sel = model.sel(c(mod.rs, mod.idist, mod.idist2), beta = "none", rank = "AIC")


# spatial autocorrelation check -------------------------------------------

# calculate Moran's I on residuals for each model
mod.rs.moran = foreach(m = 1:7, .combine = "c") %do% {
  moran = moran.mc(residuals(mod.rs[[m]]), listw = wts.rs[[m]], nsim = 99, zero.policy = T)
  moran$p.value
}
mod.rs.moran = data.frame(mod = paste0("mod.rs_", 1:7), moran.pval = mod.rs.moran)

mod.idist.moran = foreach(m = 1:length(mod.idist), .combine = "c") %do% {
  moran = moran.mc(residuals(mod.idist[[m]]), listw = wts.idist[[m]], nsim = 99, zero.policy = T)
  moran$p.value
}
mod.idist.moran = data.frame(mod = paste0("mod.idist_", 1:length(mod.idist)), moran.pval = mod.idist.moran)

mod.idist2.moran = foreach(m = 1:length(mod.idist2), .combine = "c") %do% {
  moran = moran.mc(residuals(mod.idist2[[m]]), listw = wts.idist2[[m]], nsim = 99, zero.policy = T)
  moran$p.value
}
mod.idist2.moran = data.frame(mod = paste0("mod.idist2_", 1:length(mod.idist2)), moran.pval = mod.idist2.moran)

# combine moran's I pvalue dataframes
# and filter to non-significant Moran's I
moran.df = bind_rows(mod.rs.moran, mod.idist.moran, mod.idist2.moran)
sel = sel %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "mod") %>% 
  left_join(moran.df, by = "mod") %>% 
  filter(moran.pval > 0.05) %>% 
  column_to_rownames(var = "mod")

sesvert.full = c(mod.rs, mod.idist, mod.idist2)[[rownames(sel)[1]]]
summary(sesvert.full)



# Model selection - dredge ------------------------------------------------

## set up cluster for parallel computation
## set up cluster for parallel computation
clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
clust <- try(makeCluster(getOption("cl.cores", 7), type = clusterType))

i = as.numeric(strsplit(rownames(sel)[1], split = "_")[[1]][2]) # index of best weights matrix in the list

clusterExport(clust, c("sesvert.scale", "sesvert.full", "wts.idist2", "f", "i"))
clusterEvalQ(clust, library(spatialreg))
options(na.action = "na.fail")

## run dredge to iterate through all possible variable combos
## tests all combinations of predictor variables (256 models)
sesvert.dredge = dredge(global.model = sesvert.full, beta = "none", evaluate = T, cluster = clust)

stopCluster(clust)

## run dredge to iterate through all possible variable combos
## tests all combinations of predictor variables (256 models)
#sesvert.dredge = dredge(sesvert.full, beta = "none", evaluate = T)

## calculate model average including models with delta AICc <= 2
#sesvert.avg = model.avg(sesvert.dredge, beta = "none", rank = "AICc", subset = delta <= 2, fit = T)


## second best model has deltaAICc of 4 - no need to average, just use top model (which is the full model)

## model average summary
#summary(sesvert.avg)

# save(sesvert.full, wts.rs, wts.idist, wts.idist2, mod.rs, mod.idist, mod.idist2, sesvert.dredge, file = "results/birdselton_sar_sesvert.RData")
# load("results/birdselton_sar_sesvert.RData")

save(sesvert.full, wts.rs, wts.idist, wts.idist2, mod.rs, mod.idist, mod.idist2, sesvert.dredge, file = "results/birdselton_sar_sesvert2.RData")
load("results/birdselton_sar_sesvert2.RData")


# predict -----------------------------------------------------------------

# future env data
env.future = list.files(path = "data/derivative_data/resampled_env_rasters/chelsa_future/2071_2100/ensemble/", pattern = ".tif", full.names = T)
env.future = lapply(env.future, rast)
env.future = rast(env.future)
names(env.future)[6] = "precip_ann"
env.future$temp_sea = env.future$temp_sea/100
env.future.df = terra::extract(env.future, sesvert.scale[,c('x','y')], xy = T, ID = F)

# transform variables that were previously transformed
# env.future.df = env.future.df %>% 
#   mutate(log_temp_sea = log10(temp_sea)) %>% 
#   dplyr::select(!temp_sea) %>% 
#   # scale using same scaling factors as original data
#   mutate(temp_diu = (temp_diu - mean(sesvert.t$temp_diu)) / sd(sesvert.t$temp_diu),
#          precip_wet = (precip_wet - mean(sesvert.t$precip_wet)) / sd(sesvert.t$precip_wet),
#          precip_dry = (precip_dry - mean(sesvert.t$precip_dry)) / sd(sesvert.t$precip_dry),
#          precip_sea = (precip_sea - mean(sesvert.t$precip_sea)) / sd(sesvert.t$precip_sea),
#          log_temp_sea = (log_temp_sea - mean(sesvert.t$log_temp_sea)) / sd(sesvert.t$log_temp_sea))

env.future.df = env.future.df %>% 
  mutate(log_precip_sea = log10(precip_sea)) %>% 
  dplyr::select(!precip_sea) %>% 
  # scale using same scaling factors as original data
  mutate(tmax_warm = (max_temp_warm - mean(sesvert.t$tmax_warm)) / sd(sesvert.t$tmax_warm),
         tmin_cold = (min_temp_cold - mean(sesvert.t$tmin_cold)) / sd(sesvert.t$tmin_cold),
         precip_wet = (precip_wet - mean(sesvert.t$precip_wet)) / sd(sesvert.t$precip_wet),
         precip_dry = (precip_dry - mean(sesvert.t$precip_dry)) / sd(sesvert.t$precip_dry),
         log_precip_sea = (log_precip_sea - mean(sesvert.t$log_precip_sea)) / sd(sesvert.t$log_precip_sea))


env.future.df = sesvert.scale %>% 
  dplyr::select(x, y, log_elev, veg_complexity, log_clim_velocity) %>% 
  right_join(env.future.df, by = c("x", "y")) %>% 
  dplyr::select(x,y, tmax_warm, tmin_cold, precip_wet, precip_dry, log_precip_sea, log_elev, veg_complexity, log_clim_velocity)


sesvert.dredge.sel = sesvert.dredge %>% 
  filter(delta <= 2)

# NA rows cause issue for prediction
# which rows are NA
nas = which(is.na(env.future.df$precip_dry))

# remove rows that are NA from original and prediction data frames
sesvert.scale.sub = sesvert.scale[-nas,]
env.future.df = env.future.df[-nas,]

fms <- foreach(r=1:nrow(sesvert.dredge.sel)) %do% {
  var = which(!is.na(sesvert.dredge.sel[r,2:9]))
  var = colnames(sesvert.dredge.sel[r,2:9])[var]
  ids = paste(var, collapse = " + ")
  fm = paste0("vert.mean.ses ~ ", ids)
}  

## https://stat.ethz.ch/pipermail/r-sig-geo/2009-September/006500.html

cl = makeCluster(1)
registerDoParallel(cl)

best_mods_pred <- foreach(r=1:nrow(sesvert.dredge.sel), .packages = c("spatialreg", "spdep", "broom")) %dopar% {
  
  # neighborhood weights matrix for data subset without nas
  n.sub = dnearneigh(sesvert.scale.sub[,c("x", "y")], d1 = 0, d2 = dists[i], longlat = F)
  wts.sub = nb2listwdist(neighbours = n.sub, x = sesvert.scale.sf, style = "W", type = "idw", alpha = 2, zero.policy = T)
  
  .GlobalEnv$fms <- fms
  .GlobalEnv$r <- r
  msar.err = errorsarlm(formula = fms[[r]], data = sesvert.scale, listw = wts.idist2[[i]], zero.policy = T, method = "LU", control = list(returnHCov=FALSE))
  
  # samp = sesvert.scale.sub %>%
  #   mutate(resid = resid(msar.err)) %>% 
  #   sample_n(1000)
  # sar.cor = ncf::correlog(x = samp$x, y = samp$y, z = samp$resid, increment = 111000, resamp = 99)
  # ncf:::plot.correlog(sar.cor)
  # ncf:::plot.correlog(sar.cor, xlim = c(0,5e6))
  
  moran_test = moran.mc(msar.err$residuals, listw = wts.idist2[[i]], nsim = 99, zero.policy = T)
  r2 = round(as.numeric(glance(msar.err)[1,1]),2)
  
  # in sample prediction (i.e., predicting within same spatial units that the model was fit on)
  pred = predict.sarlm.debug(msar.err, newdata = env.future.df, listw = wts.idist2[[i]], type = "TC")
  
  list(row = r, mod = msar.err, moran = moran_test, r2 = r2, pred = pred)
}

stopCluster(cl)

# save(sesvert.full, wts.rs, wts.idist, wts.idist2, mod.rs, mod.idist, mod.idist2, sesvert.dredge, best_mods_pred, file = "results/birdselton_sar_sesvert.RData")
# load("results/birdselton_sar_sesvert.RData")

save(sesvert.full, wts.rs, wts.idist, wts.idist2, mod.rs, mod.idist, mod.idist2, sesvert.dredge, best_mods_pred, file = "results/birdselton_sar_sesvert2.RData")
load("results/birdselton_sar_sesvert2.RData")


# check spatial autocorrelation
lapply(best_mods_pred, "[[", 3) # all okay
r2 = mean(sapply(best_mods_pred, "[[", 4))
sesvert.scale.sub$pred1 = as.vector(best_mods_pred[[1]][["pred"]])

# Plot results ------------------------------------------------------------

sum.avg = summary(sesvert.full)
sesvert.avg.full = sum.avg$Coef %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "var") %>% 
  rename(SE = "Std. Error", z_value = "z value", p_value = "Pr(>|z|)")

sesvert.avg.full %>% 
  filter(var != "(Intercept)") %>% 
  ggplot(aes(x = Estimate, y = reorder(var, Estimate), xmin = Estimate-SE, xmax = Estimate+SE, color = p_value<0.05)) +
  geom_pointrange() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("red3","black")) +
  scale_y_discrete("") +
  ggtitle("SES Mean Verticality: Birds elton") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))



# sesvert.scale.sub = sesvert.scale.sub %>% 
#   mutate(vert.mean.ses.future = rowMeans(select(sesvert.scale.sub, pred1:pred6)))

sesvert.scale.sub.r = rast(sesvert.scale.sub, type = "xyz", crs = "+proj=cea +datum=WGS84")
plot(sesvert.scale.sub.r$pred1)

sesvert.scale.sub = sesvert.scale.sub %>% 
  rename(vert.mean.ses.future = pred1) 

save(sesvert.full, wts.rs, wts.idist, wts.idist2, mod.rs, mod.idist, mod.idist2, 
     sesvert.dredge, best_mods_pred, sesvert.scale.sub, file = "results/birdselton_sar_sesvert2.RData")
load("results/birdselton_sar_sesvert2.RData")


mp = sesvert.scale.sub %>% 
  dplyr::select(x,y,vert.mean.ses, vert.mean.ses.future) %>% 
  pivot_longer(cols = c(vert.mean.ses, vert.mean.ses.future), names_to = "Time", values_to = "ses.vert.mean") %>% 
  ggplot() +
  geom_tile(aes(x = x, y = y, fill = ses.vert.mean)) +
  geom_sf(data = world, fill = NA) +
  coord_sf(crs = "+proj=cea +datum=WGS84") +
  scale_fill_continuous_divergingx(palette = "Spectral", rev = T, na.value = "white") +
  ggtitle("Birds elton") +
  facet_wrap(facets = ~Time, nrow = 2) +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5))

sesvertmeans = sesvert.scale.sub %>% 
  dplyr::select(x,y,vert.mean.ses.future, vert.mean.ses) %>% 
  pivot_longer(cols = c(vert.mean.ses, vert.mean.ses.future), names_to = "Time", values_to = "vert.mean.ses") %>% 
  group_by(Time) %>% 
  summarise(mean = mean(vert.mean.ses, na.rm = T))

hist = sesvert.scale.sub %>% 
  dplyr::select(x,y,vert.mean.ses.future, vert.mean.ses) %>%
  pivot_longer(cols = c(vert.mean.ses, vert.mean.ses.future), names_to = "Time", values_to = "ses.vert.mean") %>% 
  ggplot(aes(x = ses.vert.mean)) +
  geom_histogram() +
  geom_vline(data = sesvertmeans, aes(xintercept = mean), color = "red3", linetype = "dashed", size = 1) +
  scale_fill_continuous_divergingx(palette = "Spectral", rev = T, na.value = "white") +
  ggtitle("Birds elton") +
  facet_wrap(facets = ~Time, nrow = 2) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

hist + mp + plot_layout(design = "12
                                  12", widths = c(5,20), heights = c(0.5, 1))



# Birds Moura -------------------------------------------------------------

# VIF variable selection
env = read.csv("data/derivative_data/env_data.csv")

#env_vars = vif_func(in_frame = env[,c(4:18,22)], thresh = 3, trace = T)

# variable selection without hurs and vpd
env_vars = vif_func(in_frame = env[c(4:12,17:18,22)], thresh = 3, trace = T)

birds = read.csv("data/derivative_data/gridcell_data/birds_comdat/birds_breedresident_comdat_parallel.csv")

#  * - data check --------------------------------------------------------------

## SES vert model
sesvert = birds[,c("x", "y", "vert.mean.ses", "biome", "realm", env_vars)] %>% 
  drop_na()

hist(sesvert$vert.mean.ses)

# check linear association between predictor and response
as.data.frame(sesvert) %>% 
  pivot_longer(cols = 6:ncol(sesvert), names_to = "vars", values_to = "val") %>% 
  ggplot(aes(x = val, y = vert.mean.ses)) +
  geom_point(col = "grey", alpha = 0.2, size = 1) +
  geom_smooth(col = "blue", method = "lm") +
  facet_wrap(facets = vars(vars), scales = "free_x") +
  theme_classic()

# transform clim_velocity and temp seasonality and elevation
sesvert.t = as.data.frame(sesvert) %>% 
  mutate(log_clim_velocity = log10(clim_velocity),
         log_precp_dry = log10(precip_dry),
         log_temp_sea = log10(temp_sea))

sesvert.t %>% 
  pivot_longer(cols = 6:ncol(sesvert.t), names_to = "vars", values_to = "val") %>% 
  ggplot(aes(x = val, y = vert.mean.ses)) +
  geom_point(col = "grey", alpha = 0.2, size = 1) +
  geom_smooth(col = "blue", method = "lm") +
  facet_wrap(facets = vars(vars), scales = "free_x") +
  theme_classic()

sesvert.t = sesvert.t %>% 
  dplyr::select(!c(clim_velocity, log_precp_dry, temp_sea))


# scale predictor variables
sesvert.scale = sesvert.t %>% 
  mutate_at(.vars = 6:ncol(sesvert.t), .funs = function(x){scale(x)})
colnames(sesvert.scale) = colnames(sesvert.t)

# SES model
f = formula(vert.mean.ses ~ elev + log_clim_velocity + precip_dry + precip_sea + precip_wet + temp_diu + log_temp_sea + veg_complexity)

# * - selection of weights matrix ---------------------------------------------
## 1) fit full models with each weight matrix

### test neighborhood thresholds between 157km (queen's case 1st degree neighbors) and 
### 3000000m (queen's case 10th degree neighbors
### (expect spatial neighborhood to not exceed 2000km because of spautocor in response var)
### https://www.nature.com/articles/nature09329#Sec5
### simultaneously compare inverse distance and inverse distance^2 weighting
### https://www.sciencedirect.com/science/article/pii/S0006320717316890#s0010

dists = seq(157000, 3000000, 157000)
dists = dists[1:7] # try just dists 1-7 for now to speed up computation since dist 1 was the only one resulting in no more spatial autocorrelation in amphibians and mammals
sesvert.scale.sf = st_as_sf(sesvert.scale, coords = c("x","y"), crs = "+proj=cea +datum=WGS84")

cl = makeCluster(7)
registerDoParallel(cl)

# row standardized (SAR model with row standardized weights (Ver Hoef et al. 2018))
wts.rs <- foreach(i=1:length(dists), .packages = c("spdep")) %dopar% {
  n = dnearneigh(sesvert.scale[,c("x", "y")], d1 = 0, d2 = dists[i], longlat = F)
  nb2listw(neighbours = n, style = "W", zero.policy = T)
}

stopCluster(cl)

# inverse distance row standardized
cl = makeCluster(7)
registerDoParallel(cl)

# row standardized (SAR model with row standardized weights (Ver Hoef et al. 2018))
wts.idist <- foreach(i=1:length(dists), .packages = c("spdep")) %dopar% {
  n = dnearneigh(sesvert.scale[,c("x", "y")], d1 = 0, d2 = dists[i], longlat = F)
  nb2listwdist(neighbours = n, x = sesvert.scale.sf, style = "W", type = "idw", alpha = 1, zero.policy = T)
}

stopCluster(cl)


# inverse distance squared row standardized
cl = makeCluster(7)
registerDoParallel(cl)

wts.idist2 <- foreach(i=1:length(dists), .packages = c("spdep")) %dopar% {
  n = dnearneigh(sesvert.scale[,c("x", "y")], d1 = 0, d2 = dists[i], longlat = F)
  nb2listwdist(neighbours = n, x = sesvert.scale.sf, style = "W", type = "idw", alpha = 2, zero.policy = T)
}

stopCluster(cl)

### make lists to store models
### one list for each weight type

# row standardized weights
cl = makeCluster(7)
registerDoParallel(cl)

# row standardized (SAR model with row standardized weights (Ver Hoef et al. 2018))
mod.rs <- foreach(i=1:length(dists), .packages = c("spatialreg")) %dopar% {
  errorsarlm(f, data = sesvert.scale, listw = wts.rs[[i]], method = "LU", zero.policy = T, control = list(returnHcov=FALSE))
}

stopCluster(cl)

save(mod.rs, file = "results/birdsmoura_sar_mod_rs.RData")
load("results/birdsmoura_sar_mod_rs.RData")

# returnHcov=F; https://stat.ethz.ch/pipermail/r-sig-geo/2013-October/019450.html
#testmod = errorsarlm(f, data = sesvert.scale, listw = wts.rs[[1]], method = "LU", zero.policy = T, control = list(returnHcov=FALSE))
#testmod2 = spautolm(f, data = sesvert.scale, listw = wts.rs[[1]], method = "LU", family = "SAR", verbose = T)

# inverse distance row standardized weights
cl = makeCluster(7)
registerDoParallel(cl)

# row standardized (SAR model with row standardized weights (Ver Hoef et al. 2018))
mod.idist <- foreach(i=1:length(dists), .packages = c("spatialreg")) %dopar% {
  errorsarlm(f, data = sesvert.scale, listw = wts.idist[[i]], method = "LU", zero.policy = T, control = list(returnHcov=FALSE))
}

stopCluster(cl)

save(mod.idist, file = "results/birdsmoura_sar_mod_idist.RData")


# inverse distance squared row standardized weights
cl = makeCluster(7)
registerDoParallel(cl)

# row standardized (SAR model with row standardized weights (Ver Hoef et al. 2018))
mod.idist2 <- foreach(i=1:length(dists), .packages = c("spatialreg")) %dopar% {
  errorsarlm(f, data = sesvert.scale, listw = wts.idist2[[i]], method = "LU", zero.policy = T, control = list(returnHCov=FALSE))
}

stopCluster(cl)

save(mod.idist2, file = "results/birdsmoura_sar_mod_idist2.RData")


## 2) compare full models using AIC

### make lists to store models
### one list for each weight type

names(mod.rs) = paste0("mod.rs_", 1:length(mod.rs))
names(mod.idist) = paste0("mod.idist_", 1:length(mod.idist))
names(mod.idist2) = paste0("mod.idist2_", 1:length(mod.idist2))

sel = model.sel(c(mod.rs, mod.idist, mod.idist2), beta = "none", rank = "AIC")


# spatial autocorrelation check -------------------------------------------

# calculate Moran's I on residuals for each model
mod.rs.moran = foreach(m = 1:7, .combine = "c") %do% {
  moran = moran.mc(residuals(mod.rs[[m]]), listw = wts.rs[[m]], nsim = 99, zero.policy = T)
  moran$p.value
}
mod.rs.moran = data.frame(mod = paste0("mod.rs_", 1:7), moran.pval = mod.rs.moran)

mod.idist.moran = foreach(m = 1:length(mod.idist), .combine = "c") %do% {
  moran = moran.mc(residuals(mod.idist[[m]]), listw = wts.idist[[m]], nsim = 99, zero.policy = T)
  moran$p.value
}
mod.idist.moran = data.frame(mod = paste0("mod.idist_", 1:length(mod.idist)), moran.pval = mod.idist.moran)

mod.idist2.moran = foreach(m = 1:length(mod.idist2), .combine = "c") %do% {
  moran = moran.mc(residuals(mod.idist2[[m]]), listw = wts.idist2[[m]], nsim = 99, zero.policy = T)
  moran$p.value
}
mod.idist2.moran = data.frame(mod = paste0("mod.idist2_", 1:length(mod.idist2)), moran.pval = mod.idist2.moran)

# combine moran's I pvalue dataframes
# and filter to non-significant Moran's I
moran.df = bind_rows(mod.rs.moran, mod.idist.moran, mod.idist2.moran)
sel = sel %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "mod") %>% 
  left_join(moran.df, by = "mod") %>% 
  filter(moran.pval > 0.05) %>% 
  column_to_rownames(var = "mod")

sesvert.full = c(mod.rs, mod.idist, mod.idist2)[[rownames(sel)[1]]]
summary(sesvert.full)



# Model selection - dredge ------------------------------------------------

## set up cluster for parallel computation
## set up cluster for parallel computation
clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
clust <- try(makeCluster(getOption("cl.cores", 7), type = clusterType))

i = as.numeric(strsplit(rownames(sel)[1], split = "_")[[1]][2]) # index of best weights matrix in the list

clusterExport(clust, c("sesvert.scale", "sesvert.full", "wts.idist2", "f", "i"))
clusterEvalQ(clust, library(spatialreg))
options(na.action = "na.fail")

## run dredge to iterate through all possible variable combos
## tests all combinations of predictor variables (256 models)
sesvert.dredge = dredge(global.model = sesvert.full, beta = "none", evaluate = T, cluster = clust)

stopCluster(clust)

## run dredge to iterate through all possible variable combos
## tests all combinations of predictor variables (256 models)
#sesvert.dredge = dredge(sesvert.full, beta = "none", evaluate = T)

## calculate model average including models with delta AICc <= 2
sesvert.avg = model.avg(sesvert.dredge, beta = "none", rank = "AICc", subset = delta <= 2, fit = T)

## model average summary
summary(sesvert.avg)

save(sesvert.full, wts.rs, wts.idist, wts.idist2, mod.rs, mod.idist, mod.idist2, sesvert.dredge, file = "results/birdsmoura_sar_sesvert.RData")
load("results/birdsmoura_sar_sesvert.RData")


# Plot results ------------------------------------------------------------

sum.avg = summary(sesvert.full)
sesvert.avg.full = sum.avg$Coef %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "var") %>% 
  rename(SE = "Std. Error", z_value = "z value", p_value = "Pr(>|z|)")

sesvert.avg.full %>% 
  filter(var != "(Intercept)") %>% 
  ggplot(aes(x = Estimate, y = reorder(var, Estimate), xmin = Estimate-SE, xmax = Estimate+SE, color = p_value<0.05)) +
  geom_pointrange() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("black","red3")) +
  scale_y_discrete("") +
  ggtitle("SES Mean Verticality: Birds moura") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

# predict -----------------------------------------------------------------


# future env data
env.future = list.files(path = "data/derivative_data/resampled_env_rasters/chelsa_future/2071_2100/ensemble/", pattern = ".tif", full.names = T)
env.future = lapply(env.future, rast)
env.future = rast(env.future)
names(env.future)[6] = "precip_ann"
env.future$temp_sea = env.future$temp_sea/100
env.future.df = terra::extract(env.future, sesvert.scale[,c('x','y')], xy = T, ID = F)

# transform variables that were previously transformed
env.future.df = env.future.df %>% 
  mutate(log_temp_sea = log10(temp_sea)) %>% 
  dplyr::select(!temp_sea) %>% 
  # scale using same scaling factors as original data
  mutate(temp_diu = (temp_diu - mean(sesvert.t$temp_diu)) / sd(sesvert.t$temp_diu),
         precip_wet = (precip_wet - mean(sesvert.t$precip_wet)) / sd(sesvert.t$precip_wet),
         precip_dry = (precip_dry - mean(sesvert.t$precip_dry)) / sd(sesvert.t$precip_dry),
         precip_sea = (precip_sea - mean(sesvert.t$precip_sea)) / sd(sesvert.t$precip_sea),
         log_temp_sea = (log_temp_sea - mean(sesvert.t$log_temp_sea)) / sd(sesvert.t$log_temp_sea))

env.future.df = sesvert.scale %>% 
  dplyr::select(x, y, elev, veg_complexity, log_clim_velocity) %>% 
  right_join(env.future.df, by = c("x", "y")) %>% 
  dplyr::select(x,y,temp_diu, precip_wet, precip_dry, precip_sea, elev, veg_complexity, log_clim_velocity, log_temp_sea)


sesvert.dredge.sel = sesvert.dredge %>% 
  filter(delta <= 2)

# NA rows cause issue for prediction
# which rows are NA
nas = which(is.na(env.future.df$precip_dry))

# remove rows that are NA from original and prediction data frames
sesvert.scale.sub = sesvert.scale[-nas,]
env.future.df = env.future.df[-nas,]

fms <- foreach(r=1:nrow(sesvert.dredge.sel)) %do% {
  var = which(!is.na(sesvert.dredge.sel[r,2:9]))
  var = colnames(sesvert.dredge.sel[r,2:9])[var]
  ids = paste(var, collapse = " + ")
  fm = paste0("vert.mean.ses ~ ", ids)
}  

## https://stat.ethz.ch/pipermail/r-sig-geo/2009-September/006500.html

cl = makeCluster(6)
registerDoParallel(cl)

best_mods_pred <- foreach(r=1:nrow(sesvert.dredge.sel), .packages = c("spatialreg", "spdep")) %dopar% {
  
  # neighborhood weights matrix for data subset without nas
  n.sub = dnearneigh(sesvert.scale.sub[,c("x", "y")], d1 = 0, d2 = dists[i], longlat = F)
  wts.sub = nb2listwdist(neighbours = n.sub, x = sesvert.scale.sf, style = "W", type = "idw", alpha = 2, zero.policy = T)
  
  .GlobalEnv$fms <- fms
  .GlobalEnv$r <- r
  msar.err = errorsarlm(formula = fms[[r]], data = sesvert.scale, listw = wts.idist2[[i]], zero.policy = T, method = "LU", control = list(returnHCov=FALSE))
  
  # samp = sesvert.scale.sub %>%
  #   mutate(resid = resid(msar.err)) %>% 
  #   sample_n(1000)
  # sar.cor = ncf::correlog(x = samp$x, y = samp$y, z = samp$resid, increment = 111000, resamp = 99)
  # ncf:::plot.correlog(sar.cor)
  # ncf:::plot.correlog(sar.cor, xlim = c(0,5e6))
  
  moran_test = moran.mc(msar.err$residuals, listw = wts.idist2[[i]], nsim = 99, zero.policy = T)
  
  # in sample prediction (i.e., predicting within same spatial units that the model was fit on)
  pred = predict.sarlm.debug(msar.err, newdata = env.future.df, listw = wts.idist2[[i]], type = "TC")
  
  list(row = r, mod = msar.err, moran = moran_test, pred = pred)
}

stopCluster(cl)

save(sesvert.full, wts.rs, wts.idist, wts.idist2, mod.rs, mod.idist, mod.idist2, sesvert.dredge, best_mods_pred, file = "results/birdsmoura_sar_sesvert.RData")
load("results/birdsmoura_sar_sesvert.RData")

# check spatial autocorrelation
lapply(best_mods_pred, "[[", 3) # all okay
sesvert.scale.sub$pred1 = as.vector(best_mods_pred[[1]][["pred"]])
sesvert.scale.sub$pred2 = as.vector(best_mods_pred[[2]][["pred"]])
sesvert.scale.sub$pred3 = as.vector(best_mods_pred[[3]][["pred"]])
sesvert.scale.sub$pred4 = as.vector(best_mods_pred[[4]][["pred"]])

sesvert.scale.sub = sesvert.scale.sub %>% 
 mutate(vert.mean.ses.future = rowMeans(select(sesvert.scale.sub, pred1:pred4)))

sesvert.scale.sub.r = rast(sesvert.scale.sub, type = "xyz", crs = "+proj=cea +datum=WGS84")
plot(sesvert.scale.sub.r$vert.mean.ses.future)

mp = sesvert.scale.sub %>% 
  dplyr::select(x,y,vert.mean.ses.future) %>% 
  right_join(sesvert.scale, by = c("x", "y")) %>% 
  dplyr::select(x,y,vert.mean.ses, vert.mean.ses.future) %>% 
  pivot_longer(cols = c(vert.mean.ses, vert.mean.ses.future), names_to = "Time", values_to = "ses.vert.mean") %>% 
  ggplot() +
  geom_tile(aes(x = x, y = y, fill = ses.vert.mean)) +
  geom_sf(data = world, fill = NA) +
  coord_sf(crs = "+proj=cea +datum=WGS84") +
  scale_fill_continuous_divergingx(palette = "Spectral", rev = T, na.value = "white") +
  ggtitle("Birds moura") +
  facet_wrap(facets = ~Time, nrow = 2) +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5))


hist = sesvert.scale.sub %>% 
  dplyr::select(x,y,vert.mean.ses.future) %>%
  right_join(sesvert.scale, by = c("x", "y")) %>% 
  dplyr::select(x,y,vert.mean.ses, vert.mean.ses.future) %>% 
  pivot_longer(cols = c(vert.mean.ses, vert.mean.ses.future), names_to = "Time", values_to = "ses.vert.mean") %>% 
  ggplot(aes(x = ses.vert.mean)) +
  geom_histogram() +
  scale_fill_continuous_divergingx(palette = "Spectral", rev = T, na.value = "white") +
  ggtitle("Birds moura") +
  facet_wrap(facets = ~Time, nrow = 2) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

hist + mp + plot_layout(design = "12
                                  12", widths = c(5,20), heights = c(0.5, 1))



# Mean verticality Elton -----------------------------------------------------

# VIF variable selection
env = read.csv("data/derivative_data/env_data.csv")

#env_vars = vif_func(in_frame = env[,c(4:18,22)], thresh = 3, trace = T)

options(na.action = "na.omit")

# variable selection without hurs and vpd
env_vars = vif_func(in_frame = env[c(4, 7:12,17:18,22)], thresh = 5, trace = T)


birds = read.csv("data/derivative_data/gridcell_data/birds_comdat/birds_breedresident_comdat_elton.csv")

#  * - data check --------------------------------------------------------------

## vert mean model
vert = birds[,c("x", "y", "vert.mean", "biome", "realm", env_vars)] %>% 
  drop_na()

hist(vert$vert.mean)

# check linear association between predictor and response
as.data.frame(vert) %>% 
  pivot_longer(cols = 6:ncol(vert), names_to = "vars", values_to = "val") %>% 
  ggplot(aes(x = val, y = vert.mean)) +
  geom_point(col = "grey", alpha = 0.2, size = 1) +
  geom_smooth(col = "blue", method = "lm") +
  facet_wrap(facets = vars(vars), scales = "free_x") +
  theme_classic()

# transform clim_velocity and temp seasonality and elevation
vert.t = as.data.frame(vert) %>% 
  mutate(log_clim_velocity = log10(clim_velocity),
         log_elev = log10(elev),
         log_precip_dry = log10(precip_dry),
         log_precip_sea = log10(precip_sea))

vert.t %>% 
  pivot_longer(cols = 6:ncol(vert.t), names_to = "vars", values_to = "val") %>% 
  ggplot(aes(x = val, y = vert.mean)) +
  geom_point(col = "grey", alpha = 0.2, size = 1) +
  geom_smooth(col = "blue", method = "lm") +
  facet_wrap(facets = vars(vars), scales = "free_x") +
  theme_classic()

vert.t = vert.t %>% 
  dplyr::select(!c(clim_velocity, elev, precip_sea, log_precip_dry))


# scale predictor variables
vert.scale = vert.t %>% 
  mutate_at(.vars = 6:ncol(vert.t), .funs = function(x){scale(x)})
colnames(vert.scale) = colnames(vert.t)

# vert mean model
#f = formula(vert.mean ~ log_elev + log_clim_velocity + precip_dry + log_precip_sea + precip_wet + temp_diu + log_temp_sea + veg_complexity)
f = formula(vert.mean ~ log_elev + log_clim_velocity + precip_dry + log_precip_sea + precip_wet + tmax_warm + I(tmax_warm^2) + tmin_cold + veg_complexity)

vert.mod = glm(f , data = vert.scale)
plot(vert.mod)
hist(vert.mod$residuals)
summary(vert.mod)
library(car)
vif(vert.mod)

# add residuals to dataframe to test for spatial autocorrelation
vert.scale$resid = vert.mod$residuals


# * - Spatial autocorrelation ---------------------------------------------
# see https://geo200cn.github.io/spatialreg.html (very helpful tutorial)

# examine spatial autocorrelation in response var with a correlogram
## take sample to look at spatial autocorrelation
samp = vert.scale %>% sample_n(1000)
vert.cor = ncf::correlog(x = samp$x, y = samp$y, z = samp$vert.mean, increment = 111000, resamp = 99)
ncf:::plot.correlog(vert.cor)
ncf:::plot.correlog(vert.cor, xlim = c(0,5e6), ylim = c(-2,2))

# * - selection of weights matrix ---------------------------------------------
## 1) fit full models with each weight matrix

### test neighborhood thresholds between 157km (queen's case 1st degree neighbors) and 
### 3000000m (queen's case 10th degree neighbors
### (expect spatial neighborhood to not exceed 2000km because of spautocor in response var)
### https://www.nature.com/articles/nature09329#Sec5
### simultaneously compare inverse distance and inverse distance^2 weighting
### https://www.sciencedirect.com/science/article/pii/S0006320717316890#s0010

dists = seq(157000, 3000000, 157000)
# only check first 7 distances because all other taxa have used first distance and would like to save compute time
dists = dists[1:7]
vert.scale.sf = st_as_sf(vert.scale, coords = c("x","y"), crs = "+proj=cea +datum=WGS84")

cl = makeCluster(7)
registerDoParallel(cl)

# row standardized (SAR model with row standardized weights (Ver Hoef et al. 2018))
wts.rs <- foreach(i=1:length(dists), .packages = c("spdep")) %dopar% {
  n = dnearneigh(vert.scale[,c("x", "y")], d1 = 0, d2 = dists[i], longlat = F)
  nb2listw(neighbours = n, style = "W", zero.policy = T)
}

stopCluster(cl)

# inverse distance row standardized
cl = makeCluster(7)
registerDoParallel(cl)

# row standardized (SAR model with row standardized weights (Ver Hoef et al. 2018))
wts.idist <- foreach(i=1:length(dists), .packages = c("spdep")) %dopar% {
  n = dnearneigh(vert.scale[,c("x", "y")], d1 = 0, d2 = dists[i], longlat = F)
  nb2listwdist(neighbours = n, x = vert.scale.sf, style = "W", type = "idw", alpha = 1, zero.policy = T)
}

stopCluster(cl)


# inverse distance squared row standardized
cl = makeCluster(7)
registerDoParallel(cl)

wts.idist2 <- foreach(i=1:length(dists), .packages = c("spdep")) %dopar% {
  n = dnearneigh(vert.scale[,c("x", "y")], d1 = 0, d2 = dists[i], longlat = F)
  nb2listwdist(neighbours = n, x = vert.scale.sf, style = "W", type = "idw", alpha = 2, zero.policy = T)
}

stopCluster(cl)

### make lists to store models
### one list for each weight type

# row standardized weights
cl = makeCluster(7)
registerDoParallel(cl)

# row standardized (SAR model with row standardized weights (Ver Hoef et al. 2018))
mod.rs <- foreach(i=1:length(dists), .packages = c("spatialreg")) %dopar% {
  errorsarlm(f, data = vert.scale, listw = wts.rs[[i]], method = "LU", zero.policy = T, control = list(returnHcov=FALSE))
}

stopCluster(cl)

save(mod.rs, file = "results/birdselton_vertmean_sar_mod_rs.RData")

# returnHcov=F; https://stat.ethz.ch/pipermail/r-sig-geo/2013-October/019450.html
#testmod = errorsarlm(f, data = vert.scale, listw = wts.rs[[1]], method = "LU", zero.policy = T, control = list(returnHcov=FALSE))
#testmod2 = spautolm(f, data = vert.scale, listw = wts.rs[[1]], method = "LU", family = "SAR", verbose = T)

# inverse distance row standardized weights
cl = makeCluster(7)
registerDoParallel(cl)

# row standardized (SAR model with row standardized weights (Ver Hoef et al. 2018))
mod.idist <- foreach(i=1:length(dists), .packages = c("spatialreg")) %dopar% {
  errorsarlm(f, data = vert.scale, listw = wts.idist[[i]], method = "LU", zero.policy = T, control = list(returnHcov=FALSE))
}

stopCluster(cl)

save(mod.idist, file = "results/birdselton_vertmean_sar_mod_idist.RData")


# inverse distance squared row standardized weights
cl = makeCluster(7)
registerDoParallel(cl)

# row standardized (SAR model with row standardized weights (Ver Hoef et al. 2018))
mod.idist2 <- foreach(i=1:length(dists), .packages = c("spatialreg")) %dopar% {
  errorsarlm(f, data = vert.scale, listw = wts.idist2[[i]], method = "LU", zero.policy = T, control = list(returnHCov=FALSE))
}

stopCluster(cl)

save(mod.idist2, file = "results/birdselton_vertmean_sar_mod_idist2.RData")
load("results/birdselton_sar_vertmean.RData")


## 2) compare full models using AIC

### make lists to store models
### one list for each weight type

names(mod.rs) = paste0("mod.rs_", 1:length(mod.rs))
names(mod.idist) = paste0("mod.idist_", 1:length(mod.idist))
names(mod.idist2) = paste0("mod.idist2_", 1:length(mod.idist2))

sel = model.sel(c(mod.rs, mod.idist, mod.idist2), beta = "none", rank = "AIC")


# spatial autocorrelation check -------------------------------------------

# calculate Moran's I on residuals for each model
mod.rs.moran = foreach(m = 1:length(mod.rs), .combine = "c") %do% {
  moran = moran.mc(residuals(mod.rs[[m]]), listw = wts.rs[[m]], nsim = 99, zero.policy = T)
  moran$p.value
}
mod.rs.moran = data.frame(mod = paste0("mod.rs_", 1:length(mod.rs)), moran.pval = mod.rs.moran)

mod.idist.moran = foreach(m = 1:length(mod.idist), .combine = "c") %do% {
  moran = moran.mc(residuals(mod.idist[[m]]), listw = wts.idist[[m]], nsim = 99, zero.policy = T)
  moran$p.value
}
mod.idist.moran = data.frame(mod = paste0("mod.idist_", 1:length(mod.idist)), moran.pval = mod.idist.moran)

mod.idist2.moran = foreach(m = 1:length(mod.idist2), .combine = "c") %do% {
  moran = moran.mc(residuals(mod.idist2[[m]]), listw = wts.idist2[[m]], nsim = 99, zero.policy = T)
  moran$p.value
}
mod.idist2.moran = data.frame(mod = paste0("mod.idist2_", 1:length(mod.idist2)), moran.pval = mod.idist2.moran)

# spatial autocorrelation is still an issue - what to do

# combine moran's I pvalue dataframes
# and filter to non-significant Moran's I
moran.df = bind_rows(mod.rs.moran, mod.idist.moran, mod.idist2.moran)
sel = sel %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "mod") %>% 
  left_join(moran.df, by = "mod") %>% 
  filter(moran.pval > 0.05) %>% 
  column_to_rownames(var = "mod")

vert.full = c(mod.rs, mod.idist, mod.idist2)[[rownames(sel)[1]]]
summary(vert.full)



# Model selection - dredge ------------------------------------------------

## set up cluster for parallel computation
## set up cluster for parallel computation
clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
clust <- try(makeCluster(getOption("cl.cores", 7), type = clusterType))

i = as.numeric(strsplit(rownames(sel)[1], split = "_")[[1]][2]) # index of best weights matrix in the list

clusterExport(clust, c("vert.scale", "vert.full", "wts.rs", "f", "i"))
clusterEvalQ(clust, library(spatialreg))
options(na.action = "na.fail")

## run dredge to iterate through all possible variable combos
## tests all combinations of predictor variables (256 models)
vert.dredge = dredge(global.model = vert.full, beta = "none", evaluate = T, cluster = clust)

stopCluster(clust)

## run dredge to iterate through all possible variable combos
## tests all combinations of predictor variables (256 models)
#vert.dredge = dredge(vert.full, beta = "none", evaluate = T)

## calculate model average including models with delta AICc <= 2
#vert.avg = model.avg(vert.dredge, beta = "none", rank = "AICc", subset = delta <= 2, fit = T)

## model average summary
#summary(vert.avg)

# second best model is more than two AICc from best model
# vert.full is the best model

save(vert.full, wts.rs, wts.idist, wts.idist2, mod.rs, mod.idist, mod.idist2, vert.dredge, file = "results/birdselton_sar_vertmean.RData")


# Plot results ------------------------------------------------------------
#vert.avg = vert.full

sum.full = summary(vert.full)

vert.avg.full = tidy(vert.full)
r2 = round(as.numeric(glance(vert.full)[1,1]),2)

# vert.avg.full = sum.avg$coefmat.full %>% 
#   as.data.frame() %>% 
#   rownames_to_column(var = "var") %>% 
#   rename(SE = "Std. Error", z_value = "z value", p_value = "Pr(>|z|)")

vert.avg.full %>% 
  filter(term != "(Intercept)" & term != "lambda") %>% 
  ggplot(aes(x = estimate, y = reorder(term, estimate), xmin = estimate-std.error, xmax = estimate+std.error, color = p.value<0.05)) +
  geom_pointrange() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  annotate(geom = "text", label = paste0("R\u00b2 = ", r2), x = -0.003, y = 9) +
  scale_color_manual(values = c("red3")) +
  scale_y_discrete("") +
  scale_x_continuous("Estimate") +
  ggtitle("Mean Verticality: Birds (elton)") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

# predict -----------------------------------------------------------------


# future env data
env.future = list.files(path = "data/derivative_data/resampled_env_rasters/chelsa_future/2071_2100/ensemble/", pattern = ".tif", full.names = T)
env.future = lapply(env.future, rast)
env.future = rast(env.future)
names(env.future)[6] = "precip_ann"
env.future$temp_sea = env.future$temp_sea/100
env.future.df = terra::extract(env.future, vert.scale[,c('x','y')], xy = T, ID = F)

# transform variables that were previously transformed
env.future.df = env.future.df %>% 
  # scale using same scaling factors as original data
  mutate(log_precip_sea = log10(precip_sea),
         tmax_warm = (max_temp_warm - mean(vert.t$tmax_warm)) / sd(vert.t$tmax_warm),
         tmin_cold = (min_temp_cold - mean(vert.t$tmin_cold)) / sd(vert.t$tmin_cold),
         precip_wet = (precip_wet - mean(vert.t$precip_wet)) / sd(vert.t$precip_wet),
         precip_dry = (precip_dry - mean(vert.t$precip_dry)) / sd(vert.t$precip_dry),
         log_precip_sea = (log_precip_sea - mean(vert.t$log_precip_sea)) / sd(vert.t$log_precip_sea))

env.future.df = vert.scale %>% 
  dplyr::select(x, y, log_elev, veg_complexity, log_clim_velocity) %>% 
  right_join(env.future.df, by = c("x", "y")) %>% 
  dplyr::select(x,y,tmax_warm, tmin_cold, precip_wet, precip_dry, log_precip_sea, log_elev, veg_complexity, log_clim_velocity)


vert.dredge.sel = vert.dredge %>% 
  filter(delta <= 2)

# NA rows cause issue for prediction
# which rows are NA
nas = which(is.na(env.future.df$precip_dry))

# remove rows that are NA from original and prediction data frames
vert.scale.sub = vert.scale[-nas,]
env.future.df = env.future.df[-nas,]

fms <- foreach(r=1:nrow(vert.dredge.sel)) %do% {
  var = which(!is.na(vert.dredge.sel[r,2:9]))
  var = colnames(vert.dredge.sel[r,2:9])[var]
  ids = paste(var, collapse = " + ")
  fm = paste0("vert.mean ~ ", ids)
}  

## https://stat.ethz.ch/pipermail/r-sig-geo/2009-September/006500.html

cl = makeCluster(1)
registerDoParallel(cl)

best_mods_pred <- foreach(r=1:nrow(vert.dredge.sel), .packages = c("spatialreg", "spdep", "broom")) %dopar% {
  
  # neighborhood weights matrix for data subset without nas
  n.sub = dnearneigh(vert.scale.sub[,c("x", "y")], d1 = 0, d2 = dists[i], longlat = F)
  wts.sub = nb2listwdist(neighbours = n.sub, x = vert.scale.sf, style = "W", type = "idw", alpha = 2, zero.policy = T)
  
  .GlobalEnv$fms <- fms
  .GlobalEnv$r <- r
  msar.err = errorsarlm(formula = fms[[r]], data = vert.scale, listw = wts.rs[[i]], zero.policy = T, method = "LU", control = list(returnHCov=FALSE))
  
  # samp = vert.scale.sub %>%
  #   mutate(resid = resid(msar.err)) %>% 
  #   sample_n(1000)
  # sar.cor = ncf::correlog(x = samp$x, y = samp$y, z = samp$resid, increment = 111000, resamp = 99)
  # ncf:::plot.correlog(sar.cor)
  # ncf:::plot.correlog(sar.cor, xlim = c(0,5e6))
  
  moran_test = moran.mc(msar.err$residuals, listw = wts.rs[[i]], nsim = 99, zero.policy = T)
  r2 = round(as.numeric(glance(msar.err)[1,1]),2)
  
  # in sample prediction (i.e., predicting within same spatial units that the model was fit on)
  pred = predict.sarlm.debug(msar.err, newdata = env.future.df, listw = wts.idist2[[i]], type = "TC")
  
  list(row = r, mod = msar.err, moran = moran_test, r2 = r2, pred = pred)
}

stopCluster(cl)

save(vert.full, wts.rs, wts.idist, wts.idist2, mod.rs, mod.idist, mod.idist2, vert.dredge, best_mods_pred, file = "results/birdselton_sar_vertmean.RData")
load("results/birdselton_sar_vertmean.RData")

# check spatial autocorrelation
lapply(best_mods_pred, "[[", 3) # all okay
r2 = sapply(best_mods_pred, "[[", 4)
vert.scale.sub$pred1 = as.vector(best_mods_pred[[1]][["pred"]])

vert.scale.sub = vert.scale.sub %>% 
  mutate(vert.mean.future = pred1)

vert.scale.sub.r = rast(vert.scale.sub, type = "xyz", crs = "+proj=cea +datum=WGS84")
plot(vert.scale.sub.r$vert.mean.future)

vert.scale.sub = vert.scale.sub %>% 
  dplyr::select(x,y,vert.mean.future) %>% 
  right_join(vert.scale, by = c("x", "y"))

mp = vert.scale.sub %>%
  pivot_longer(cols = c(vert.mean, vert.mean.future), names_to = "Time", values_to = "vert.mean") %>% 
  ggplot() +
  geom_tile(aes(x = x, y = y, fill = vert.mean)) +
  geom_sf(data = world, fill = NA) +
  coord_sf(crs = "+proj=cea +datum=WGS84") +
  scale_fill_continuous_sequential(palette = "Viridis", rev = T, na.value = "white") +
  #ggtitle("Amphibians mean verticality (moura)") +
  facet_wrap(facets = ~Time, nrow = 2) +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5))

vertmeans = vert.scale.sub %>% 
  dplyr::select(x,y,vert.mean.future) %>% 
  right_join(vert.scale, by = c("x", "y")) %>% 
  dplyr::select(x,y,vert.mean, vert.mean.future) %>% 
  pivot_longer(cols = c(vert.mean, vert.mean.future), names_to = "Time", values_to = "vert.mean") %>% 
  group_by(Time) %>% 
  summarise(mean = mean(vert.mean, na.rm = T))

hist = vert.scale.sub %>% 
  dplyr::select(x,y,vert.mean, vert.mean.future) %>% 
  pivot_longer(cols = c(vert.mean, vert.mean.future), names_to = "Time", values_to = "vert.mean") %>% 
  ggplot(aes(x = vert.mean)) +
  geom_histogram() +
  geom_vline(data = vertmeans, aes(xintercept = mean), color = "red3", linetype = "dashed", size = 1) +
  scale_fill_continuous_sequential(palette = "YlGnBu", rev = T, na.value = "white") +
  #ggtitle("Amphibians (moura)") +
  facet_wrap(facets = ~Time, nrow = 2) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

hist + mp + plot_layout(design = "12
                                  12", widths = c(5,20), heights = c(0.5, 1)) +
  plot_annotation(title = "Birds: mean verticality (elton)")


save(vert.full, wts.rs, wts.idist, wts.idist2, mod.rs, mod.idist, mod.idist2, vert.dredge,
     best_mods_pred, vert.scale.sub, file = "results/birdselton_sar_vertmean.RData")
