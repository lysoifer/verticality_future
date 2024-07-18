library(sdmTMB)
library(foreach)
library(doParallel)

load("results/sdmTMB_models/amphibians_meanvert.RData")
visreg(m.final, "tmax_warm")

dat = m.final$data
dat %>% 
  mutate(tt = ifelse(abs(y) < 23.5, "Tropical", "Temperate")) %>% 
  ggplot(aes(x = vert.mean.ses, y = rich, color = tt)) +
  geom_point(, size = 0.05) +
  geom_smooth()

conditional_effects = function(m, random_int, ncore = 7, fpath) {
  # m: sdmTMB model object
  # random_int: name of variable used as a random intercept in the model
  # ncore: number of cores for parallelization
  # fpath: file path to save predictions to
  
  dat = m$data
  terms = as.character(m$formula[[1]])[3]
  terms = strsplit(terms, split = "\\+")[[1]]
  terms = grep("\\(", terms, invert = T, value = T)
  terms = str_trim(terms, side = "both")
  
  dat = dat[,which(colnames(dat) %in% terms)]
  
  # make matrix with median values
  newdat = apply(dat, 2, median)
  newdat = matrix(newdat, nrow = 100, ncol = length(newdat), byrow = T)
  colnames(newdat) = colnames(dat)
  
  if(!is.null(random_int)){
    randint = unique(m$data[,random_int])
  }
  
  # replace variables one by one to get conditional predictions
  cl = makeCluster(ncore)
  registerDoParallel(cl)
  cond_resp = foreach(i = 1:ncol(dat)) %dopar% {
    var = seq(min(dat[,i]), max(dat[,i]), length.out = 100)
    n = newdat
    n[,i] = var
    n = as.data.frame(n)
    p = foreach(r = randint, .combine = "rbind") %do% {
      n[,random_int] = rep(r, 100)
      
      # predict to newdata
      pred = sdmTMB:::predict.sdmTMB(m, newdata = n, type = "response", re_form = NA, se_fit = F)
    }
    
  }
  stopCluster(cl)
  names(cond_resp) = colnames(dat)
  save(cond_resp, file = fpath)
}


# Unscale the data --------------------------------------------------------

# amphibian forest only 50km resolution
raw = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/amph_comdat.csv")

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

dat.t = raw %>% 
  mutate(log_clim_velocity = log10(clim_velocity),
         log_tmin_cold = log10(tmin_cold + 50),
         log_precip_dry = log10(precip_dry),
         log_precip_wet = log10(precip_wet))


ggplot(cond_resp[[3]]) +
  geom_line(aes(x = tmax_warm, y = est, colour = ecoregion)) +
  theme(legend.position = "none")

##RESCALE X AXIS#
ggplot(cond_resp[[7]]) +
  geom_line(aes(x = 10^(log_precip_dry*), y = est, colour = ecoregion)) +
  theme(legend.position = "none")

ggplot(cond_resp[[1]]) +
  geom_line(aes(x = canopy_height, y = est, colour = ecoregion)) +
  theme(legend.position = "none")
















