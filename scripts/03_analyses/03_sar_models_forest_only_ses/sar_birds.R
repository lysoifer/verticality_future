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
library(doParallel)
library(terra)
library(colorspace)
library(patchwork)
library(broom)

options(na.action = "na.omit")


# VIF variable selection
env = read.csv("data/derivative_data/env_data.csv")

# variable selection 
env_vars = vif_func(in_frame = env[c(4, 7:12,17:18,22)], thresh = 5, trace = T)



# Birds --------------------------------------------------------------

birds = read.csv("data/derivative_data/gridcell_data/birds_comdat/birds_comdat_parallel_elton_forestsOnly.csv")

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

sesvert.t = as.data.frame(sesvert) %>% 
  mutate(log_clim_velocity = log10(clim_velocity),
         log_elev = log10(elev),
         log_precip_sea = log10(precip_sea),
         log_precip_dry = log10(precip_dry+1))

sesvert.t %>% 
  pivot_longer(cols = 6:ncol(sesvert.t), names_to = "vars", values_to = "val") %>% 
  ggplot(aes(x = val, y = vert.mean.ses)) +
  geom_point(col = "grey", alpha = 0.2, size = 1) +
  geom_smooth(col = "blue", method = "lm") +
  facet_wrap(facets = vars(vars), scales = "free_x") +
  theme_classic()

# keep log transformation for climate velocity and elevation
sesvert.t = sesvert.t %>% 
  dplyr::select(!c(clim_velocity, elev, precip_sea, precip_dry))

# scale predictor variables
sesvert.scale = sesvert.t %>% 
  mutate_at(.vars = 6:ncol(sesvert.t), .funs = function(x){scale(x)})
colnames(sesvert.scale) = colnames(sesvert.t)

sesvert.scale$id = 1:nrow(sesvert.scale)

# SES model
# log_precip_dry = log10(precip_dry + 1)
f = formula(vert.mean.ses ~ log_elev + log_clim_velocity + I(tmax_warm^2) + tmax_warm + 
              tmin_cold + log_precip_dry + log_precip_sea + precip_wet + veg_complexity)

# * - selection of weights matrix ---------------------------------------------
## 1) fit full models with each weight matrix

### test neighborhood thresholds between 157km (queen's case 1st degree neighbors) and 
### 5570km (queen's case 10th degree neighbors
### (expect spatial neighborhood to not exceed 2000km because of spautocor in response var)
### https://www.nature.com/articles/nature09329#Sec5
### simultaneously compare inverse distance and inverse distance^2 weighting
### https://www.sciencedirect.com/science/article/pii/S0006320717316890#s0010

dists = seq(157000, 2000000, 157000)
dists = dists[1:7]
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

# row standardized weights
cl = makeCluster(7)
registerDoParallel(cl)

# row standardized (SAR model with row standardized weights (Ver Hoef et al. 2018))
mod.rs <- foreach(i=1:length(dists), .packages = c("spatialreg")) %dopar% {
  errorsarlm(f, data = sesvert.scale, listw = wts.rs[[i]], method = "LU", 
             zero.policy = T, control = list(pWOrder=1000, returnHcov=TRUE))
}

stopCluster(cl)


# inverse distance row standardized weights
cl = makeCluster(7)
registerDoParallel(cl)

# row standardized (SAR model with row standardized weights (Ver Hoef et al. 2018))
mod.idist <- foreach(i=1:length(dists), .packages = c("spatialreg")) %dopar% {
  errorsarlm(f, data = sesvert.scale, listw = wts.idist[[i]], method = "LU", 
             zero.policy = T, control = list(pWOrder=1000, returnHcov=TRUE))
}

stopCluster(cl)

#save(mod.idist, file = "results/sar_mods_forest_only/amphibians/amphibians_sar_mod_idist.RData")


# inverse distance squared row standardized weights
cl = makeCluster(7)
registerDoParallel(cl)

mod.idist2 <- foreach(i=1:length(dists), .packages = c("spatialreg")) %dopar% {
  errorsarlm(f, data = sesvert.scale, listw = wts.idist2[[i]], method = "LU", 
             zero.policy = T, control = list(pWOrder=1000, returnHcov=TRUE))
}

stopCluster(cl)

#save(mod.idist2, file = "results/sar_mods_forest_only/amphibians/amphibians_sar_mod_idist2.RData")


## 2) compare full models using AIC

### make lists to store models
### one list for each weight type

names(mod.rs) = paste0("mod.rs_", 1:length(mod.rs))
names(mod.idist) = paste0("mod.idist_", 1:length(mod.idist))
names(mod.idist2) = paste0("mod.idist2_", 1:length(mod.idist2))

sel = model.sel(c(mod.rs, mod.idist, mod.idist2), beta = "none", rank = "AICc")


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

# fit all of the best models so inspect for lambdas and spatial autocorrelation
sesvert.dredge.sel = sesvert.dredge %>% 
  filter(delta <= 2)

fms <- foreach(r=1:nrow(sesvert.dredge.sel)) %do% {
  var = which(!is.na(sesvert.dredge.sel[r,2:9]))
  var = colnames(sesvert.dredge.sel[r,2:9])[var]
  ids = paste(var, collapse = " + ")
  fm = paste0("vert.mean.ses ~ ", ids)
}  

cl = makeCluster(7)
registerDoParallel(cl)

# these are the best models identified by dredge function with deltaAICc <= 2
best_mods <- foreach(r=1:nrow(sesvert.dredge.sel), .packages = c("spatialreg", "spdep", "broom")) %dopar% {
  
  # neighborhood weights matrix for data subset without nas
  .GlobalEnv$fms <- fms
  .GlobalEnv$r <- r
  msar.err = errorsarlm(formula = fms[[r]], data = sesvert.scale, listw = wts.idist2[[i]], zero.policy = T, method = "LU", control = list(pWOrder=500, returnHcov=TRUE))
  
  return(msar.err)
}

stopCluster(cl)


## calculate model average including models with delta AICc <= 2
sesvert.avg = model.avg(sesvert.dredge, beta = "none", rank = "AICc", subset = delta <= 2, fit = T)

## model average summary
summary(sesvert.avg)

# predict model
pred = predict(sesvert.avg)

# residuals
resid = pred - sesvert.scale$vert.mean.ses
moran.mc(resid, wts.idist2[[i]], nsim = 100)

# pseudo-r2 (see https://onlinelibrary.wiley.com/doi/full/10.1111/j.1466-8238.2007.00334.x)
r2 = cor(pred, sesvert.scale$vert.mean.ses)^2

# predict model with only the trend
pred.trend = predict(sesvert.avg, newdata = sesvert.scale, type = "response", 
                     pred.type = "TS", listw = wts.idist2[[2]])
r2.trend = cor(pred.trend, sesvert.scale$vert.mean.ses)

save.image(file = "results/sar_mods_forestOnly_forestSES/birds/birds_sar_sesvert.RData")
load("results/sar_mods_forestOnly_forestSES/birds/birds_sar_sesvert.RData")


plot(sesvert.scale$vert.mean.ses, pred)
plot(pred, sesvert.scale$vert.mean.ses-pred)

plot(sesvert.scale$vert.mean.ses, pred.trend)
plot(pred.trend, sesvert.scale$vert.mean.ses-pred.trend)

# predict -----------------------------------------------------------------

# future env data
env.future = list.files(path = "data/derivative_data/resampled_env_rasters/chelsa_future/2071_2100/ensemble/", pattern = ".tif", full.names = T)
env.future = lapply(env.future, rast)
env.future = rast(env.future)
names(env.future)[6] = "precip_ann"
env.future$temp_sea = env.future$temp_sea/100
env.future.df = terra::extract(env.future, sesvert.scale[,c('x','y')], xy = T, ID = F)

env.future.df = env.future.df %>% 
  mutate(log_precip_dry = log10(precip_dry + 1),
         log_precip_sea = log10(precip_sea)) %>% 
  dplyr::select(!c(precip_dry, precip_sea)) %>%  
  # scale using same scaling factors as original data
  mutate(tmax_warm = (max_temp_warm - mean(sesvert.t$tmax_warm)) / sd(sesvert.t$tmax_warm),
         tmin_cold = (min_temp_cold - mean(sesvert.t$tmin_cold)) / sd(sesvert.t$tmin_cold),
         precip_wet = (precip_wet - mean(sesvert.t$precip_wet)) / sd(sesvert.t$precip_wet),
         log_precip_dry = ((log_precip_dry) - mean(sesvert.t$log_precip_dry)) / sd(sesvert.t$log_precip_dry),
         log_precip_sea = (log_precip_sea - mean(sesvert.t$log_precip_sea)) / sd(sesvert.t$log_precip_sea))
env.future.df$id = 1:nrow(env.future.df)


env.future.df = sesvert.scale %>% 
  dplyr::select(x, y, id, log_elev, veg_complexity, log_clim_velocity) %>% 
  right_join(env.future.df, by = c("x", "y", "id")) %>% 
  dplyr::select(id, x,y, tmin_cold, tmax_warm, precip_wet, log_precip_dry, log_precip_sea,
                log_elev, veg_complexity, log_clim_velocity)

# NA rows cause issue for prediction
# which rows are NA
nas = which(is.na(env.future.df$log_precip_dry))
env.future.df = env.future.df[-nas,]
rownames(env.future.df) = env.future.df$id
env.future.sf = st_as_sf(env.future.df, coords = c("x", "y"), crs = "+proj=cea +datum=WGS84")

# Predict averaged model
neigh = dnearneigh(env.future.df[,c("x", "y")], d1 = 0, d2 = dists[i], longlat = F, row.names = rownames(env.future.df))
wts.pred = nb2listwdist(neighbours = neigh, x = env.future.sf, style = "W", type = "idw", alpha = 2, zero.policy = T)

pred.future = MuMIn:::predict.averaging(sesvert.avg, newdata = env.future.df, type = "response", 
                                        pred.type = "TS", listw = wts.pred)

env.future.df$pred.future = pred.future
sesvert.scale$pred.pres.trend = pred.trend
sesvert.future = env.future.df %>% dplyr::select(x,y,pred.future)
sesvert.pres = sesvert.scale %>% dplyr::select(x,y,vert.mean.ses, pred.pres.trend)

pred.df = left_join(sesvert.pres, sesvert.future, by = c('x', "y"))

## https://stat.ethz.ch/pipermail/r-sig-geo/2009-September/006500.html

save.image(file = "results/sar_mods_forestOnly_forestSES/birds/birds_sar_sesvert.RData")
load("results/sar_mods_forestOnly_forestSES/birds/birds_sar_sesvert.RData")


# Plot results ------------------------------------------------------------

sum.avg = summary(sesvert.avg)
sesvert.avg.full = sum.avg$coefmat.full %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "var") %>% 
  rename(SE = "Std. Error", z_value = "z value", p_value = "Pr(>|z|)")

sesvert.avg.full %>% 
  filter(var != "(Intercept)") %>% 
  ggplot(aes(x = Estimate, y = reorder(var, Estimate), xmin = Estimate-SE, xmax = Estimate+SE, color = p_value<0.05)) +
  geom_pointrange() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("black", "red3")) +
  scale_y_discrete("") +
  ggtitle("SES Mean Verticality: Birds") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

pred.r = rast(pred.df, type = "xyz", crs = "+proj=cea +datum=WGS84")
plot(pred.r$pred.future)
plot(pred.r$pred.future - pred.r$pred.pres.trend)

# Mean verticality -----------------------------------------------------

# VIF variable selection
env = read.csv("data/derivative_data/env_data.csv")

#env_vars = vif_func(in_frame = env[,c(4:18,22)], thresh = 3, trace = T)
options(na.action = "na.omit")

# variable selection without hurs and vpd
#env_vars = vif_func(in_frame = env[c(4:12,17:18,22)], thresh = 3, trace = T)
env_vars = vif_func(in_frame = env[c(4, 7:12,17:18,22)], thresh = 5, trace = T)


birds = read.csv("data/derivative_data/gridcell_data/birds_comdat/birds_comdat_parallel_elton_forestsOnly.csv")


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

vert.t = as.data.frame(vert) %>% 
  mutate(log_clim_velocity = log10(clim_velocity),
         log_elev = log10(elev),
         log_precip_dry = log10(precip_dry+1),
         log_precip_wet = log10(precip_wet),
         log_veg_complexity = log10(veg_complexity),
         log_precip_sea = log10(precip_sea))


vert.t %>% 
  pivot_longer(cols = 6:ncol(vert.t), names_to = "vars", values_to = "val") %>% 
  ggplot(aes(x = val, y = vert.mean)) +
  geom_point(col = "grey", alpha = 0.2, size = 1) +
  geom_smooth(col = "blue", method = "lm") +
  facet_wrap(facets = vars(vars), scales = "free_x") +
  theme_classic()

vert.t = vert.t %>% 
  dplyr::select(!c(clim_velocity, elev, log_precip_wet, precip_dry, log_veg_complexity, precip_sea))

# scale predictor variables
vert.scale = vert.t %>% 
  mutate_at(.vars = 6:ncol(vert.t), .funs = function(x){scale(x)})
colnames(vert.scale) = colnames(vert.t)

vert.scale$id = 1:nrow(vert.scale)

# vert mean model# vertvert.scale mean model
f = formula(vert.mean ~ log_elev + log_clim_velocity + I(tmax_warm^2) + tmax_warm + tmin_cold + log_precip_dry + log_precip_sea + precip_wet + veg_complexity)

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
  errorsarlm(f, data = vert.scale, listw = wts.rs[[i]], method = "LU", zero.policy = T, 
             control = list(pWOrder=1000, returnHcov=T))
}

stopCluster(cl)

# inverse distance row standardized weights
cl = makeCluster(7)
registerDoParallel(cl)

# row standardized (SAR model with row standardized weights (Ver Hoef et al. 2018))
mod.idist <- foreach(i=1:length(dists), .packages = c("spatialreg")) %dopar% {
  errorsarlm(f, data = vert.scale, listw = wts.idist[[i]], method = "LU", zero.policy = T, 
             control = list(pWOrder=1000, returnHcov=TRUE))
}

stopCluster(cl)

# inverse distance squared row standardized weights
cl = makeCluster(7)
registerDoParallel(cl)

# row standardized (SAR model with row standardized weights (Ver Hoef et al. 2018))
mod.idist2 <- foreach(i=1:length(dists), .packages = c("spatialreg")) %dopar% {
  errorsarlm(f, data = vert.scale, listw = wts.idist2[[i]], method = "LU", zero.policy = T, 
             control = list(pWOrder=1000, returnHcov=TRUE))
}

stopCluster(cl)

## 2) compare full models using AICc

### make lists to store models
### one list for each weight type

names(mod.rs) = paste0("mod.rs_", 1:length(mod.rs))
names(mod.idist) = paste0("mod.idist_", 1:length(mod.idist))
names(mod.idist2) = paste0("mod.idist2_", 1:length(mod.idist2))

sel = model.sel(c(mod.rs, mod.idist, mod.idist2), beta = "none", rank = "AICc")


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

clusterExport(clust, c("vert.scale", "vert.full", "wts.idist2", "f", "i"))
clusterEvalQ(clust, library(spatialreg))
options(na.action = "na.fail")

## run dredge to iterate through all possible variable combos
## tests all combinations of predictor variables (256 models)
vert.dredge = dredge(global.model = vert.full, beta = "none", evaluate = T, cluster = clust)

stopCluster(clust)

# fit all of the best models so inspect for lambdas and spatial autocorrelation
vert.dredge.sel = vert.dredge %>% 
  filter(delta <= 2)

fms <- foreach(r=1:nrow(vert.dredge.sel)) %do% {
  var = which(!is.na(vert.dredge.sel[r,2:9]))
  var = colnames(vert.dredge.sel[r,2:9])[var]
  ids = paste(var, collapse = " + ")
  fm = paste0("vert.mean ~ ", ids)
}  

cl = makeCluster(7)
registerDoParallel(cl)

# these are the best models identified by dredge function with deltaAICc <= 2
best_mods <- foreach(r=1:nrow(vert.dredge.sel), .packages = c("spatialreg", "spdep", "broom")) %dopar% {
  
  # neighborhood weights matrix for data subset without nas
  .GlobalEnv$fms <- fms
  .GlobalEnv$r <- r
  msar.err = errorsarlm(formula = fms[[r]], data = vert.scale, listw = wts.idist2[[i]], zero.policy = T, method = "LU", control = list(pWOrder=500, returnHcov=TRUE))
  
  return(msar.err)
}

stopCluster(cl)

## calculate model average including models with delta AICc <= 2
vert.avg = model.avg(vert.dredge, beta = "none", rank = "AICc", subset = delta <= 2, fit = T)

## model average summary
summary(vert.avg)

# predict model
pred = predict(vert.avg)

# residuals
resid = pred - vert.scale$vert.mean
moran.mc(resid, wts.idist2[[i]], nsim = 100)

# pseudo-r2 (see https://onlinelibrary.wiley.com/doi/full/10.1111/j.1466-8238.2007.00334.x)
r2 = cor(pred, vert.scale$vert.mean)^2

# predict model with only the trend
pred.trend = predict(vert.avg, newdata = vert.scale, type = "response", 
                     pred.type = "TS", listw = wts.idist2[[2]])
r2.trend = cor(pred.trend, vert.scale$vert.mean)

save.image(file = "results/sar_mods_forestOnly_forestSES/birds/birds_sar_meanvert.RData")

plot(vert.scale$vert.mean, pred)
plot(pred, vert.scale$vert.mean-pred)

plot(vert.scale$vert.mean, pred.trend)
plot(pred.trend, vert.scale$vert.mean-pred.trend)


# predict -----------------------------------------------------------------

# future env data
env.future = list.files(path = "data/derivative_data/resampled_env_rasters/chelsa_future/2071_2100/ensemble/", pattern = ".tif", full.names = T)
env.future = lapply(env.future, rast)
env.future = rast(env.future)
names(env.future)[6] = "precip_ann"
env.future$temp_sea = env.future$temp_sea/100
env.future.df = terra::extract(env.future, vert.scale[,c('x','y')], xy = T, ID = F)


env.future.df = env.future.df %>% 
  mutate(log_precip_dry = log10(precip_dry+1),
         log_precip_sea = log10(precip_sea)) %>% 
  # scale using same scaling factors as original data
  mutate(tmax_warm = (max_temp_warm - mean(vert.t$tmax_warm)) / sd(vert.t$tmax_warm),
         tmin_cold = (min_temp_cold - mean(vert.t$tmin_cold)) / sd(vert.t$tmin_cold),
         precip_wet = (precip_wet - mean(vert.t$precip_wet)) / sd(vert.t$precip_wet),
         log_precip_dry = (log_precip_dry - mean(vert.t$log_precip_dry)) / sd(vert.t$log_precip_dry),
         log_precip_sea = (log_precip_sea - mean(vert.t$log_precip_sea)) / sd(vert.t$log_precip_sea))
env.future.df$id = 1:nrow(env.future.df)


env.future.df = vert.scale %>% 
  dplyr::select(x, y, log_elev, veg_complexity, log_clim_velocity) %>% 
  right_join(env.future.df, by = c("x", "y")) %>% 
  dplyr::select(id, x,y, tmin_cold, tmax_warm, precip_wet, log_precip_dry, log_precip_sea, 
                log_elev, veg_complexity, log_clim_velocity)

# NA rows cause issue for prediction
# which rows are NA
nas = which(is.na(env.future.df$log_precip_dry))
env.future.df = env.future.df[-nas,]
rownames(env.future.df) = env.future.df$id
env.future.sf = st_as_sf(env.future.df, coords = c("x", "y"), crs = "+proj=cea +datum=WGS84")

# Predict averaged model
neigh = dnearneigh(env.future.df[,c("x", "y")], d1 = 0, d2 = dists[i], longlat = F, row.names = rownames(env.future.df))
wts.pred = nb2listwdist(neighbours = neigh, x = env.future.sf, style = "W", type = "idw", alpha = 2, zero.policy = T)

pred.future = MuMIn:::predict.averaging(vert.avg, newdata = env.future.df, type = "response", 
                                        pred.type = "TS", listw = wts.pred)

# make dataframe with predicted values
env.future.df$pred.future = pred.future
vert.scale$pred.pres.trend = pred.trend # predictions to present data - only trend (i.e., does not include spatial signal)
vert.future = env.future.df %>% dplyr::select(x,y,pred.future)
vert.pres = vert.scale %>% dplyr::select(x,y,vert.mean, pred.pres.trend)


pred.df = left_join(vert.pres, vert.future, by = c('x', "y"))


save.image(file = "results/sar_mods_forestOnly_forestSES/birds/birds_sar_meanvert.RData")
load("results/sar_mods_forestOnly_forestSES/birds/birds_sar_meanvert.RData")


# Plot results ------------------------------------------------------------

sum.avg = summary(vert.avg)

vert.avg.full = sum.avg$coefmat.full %>%
  as.data.frame() %>%
  rownames_to_column(var = "var") %>%
  rename(SE = "Std. Error", z_value = "z value", p_value = "Pr(>|z|)")


vert.avg.full %>% 
  filter(var != "(Intercept)") %>% 
  ggplot(aes(x = Estimate, y = reorder(var, Estimate), xmin = Estimate-SE, xmax = Estimate+SE, color = p_value<0.05)) +
  geom_pointrange() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  annotate(geom = "text", x = -0.01, y = 7, label = paste0("R\u00b2 = ", r2)) +
  scale_color_manual(values = c("black", "red3")) +
  scale_y_discrete("") +
  ggtitle("Mean Verticality: Amphibians (moura)") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))


pred.r = rast(pred.df, type = "xyz", crs = "+proj=cea +datum=WGS84")
plot(pred.r$pred.future)
plot(pred.r$pred.future - pred.r$pred.pres.trend)

