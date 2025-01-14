library(terra)
library(tidyverse)
library(tidyterra)

# future env data
env.f = rast(list.files(path = "data/derivative_data/resampled_env_rasters_50km/chelsa_future/2071_2100/ensemble/",
                        pattern = ".tif", full.names = T))
plot(env.f$precip_dry)

# current env data
env = read.csv("data/derivative_data/env_data_50km_forest.csv")
env = env %>% 
  dplyr::select(x,y,precip_dry, precip_wet, tmax_warm, tmin_cold, precip_sea) %>% 
  rast(crs = "+proj=cea +datum=WGS84")
plot(env$precip_dry)

v = vect("data/original/rnaturalearth_world.shp")
v = project(v, "+proj=cea +datum=WGS84")

env.f = crop(env.f, env)

tmax.change = env.f$max_temp_warm - env$tmax_warm
tmin.change = env.f$min_temp_cold - env$tmin_cold
precipdry.change = env.f$precip_dry - env$precip_dry
precipwet.change = env.f$precip_wet - env$precip_wet

plot(tmax.change)
plot(tmin.change)
plot(precipdry.change)
plot(precipwet.change)

