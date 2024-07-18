library(terra)
library(tidyverse)

env = read.csv("data/derivative_data/env_data.csv")

env.future = list.files(path = "data/derivative_data/resampled_env_rasters/chelsa_future/2071_2100/ensemble/", pattern = ".tif", full.names = T)
env.future = lapply(env.future, rast)
env.future = rast(env.future)
names(env.future)[6] = "precip_ann"
env.future$temp_sea = env.future$temp_sea/100
env.future = as.data.frame(env.future, xy = T)
env.future = env.future %>% 
  rename_at(.vars = 3:11, .funs = function(x) paste0(x, "_future"))

load("results/amphibians_sar_vertmean2.RData")
amph.vert = vert.scale.sub %>% 
  dplyr::select(x,y,vert.mean, vert.mean.future) %>% 
  left_join(env, by = c("x", "y")) %>% 
  left_join(env.future, by = c("x", "y"))

ggplot(amph.vert) +
  geom_point(aes(x = precip_sea, y = precip_sea_future, color = vert.mean.future - vert.mean), size = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  #facet_wrap(~realm) +
  scale_color_continuous_divergingx("spectral") +
  theme_classic()

ggplot(amph.vert) +
  geom_point(aes(x = tmin_cold, y = min_temp_cold_future, color = vert.mean.future - vert.mean), size = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  facet_wrap(~realm) +
  scale_color_continuous_divergingx("spectral") +
  theme_classic()
