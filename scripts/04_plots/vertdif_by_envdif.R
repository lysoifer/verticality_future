library(terra)
library(tidyterra)
library(tidyverse)

load("results/amphibians_sar_sesvert2.RData")
amph = sesvert.scale.sub

amph.pdat = amph %>% 
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  #project("epsg:4326") %>% 
  as.data.frame(xy = T) %>% 
  mutate(sesvert.dif = vert.mean.ses.future-vert.mean.ses) %>% 
  dplyr::select(x,y,sesvert.dif) %>% 
  rename(sesvert.dif.amph = sesvert.dif)


load("results/mammals_sar_sesvert2.RData")
mammals = sesvert.scale.sub

mammals.pdat = mammals %>% 
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  #project("epsg:4326") %>% 
  as.data.frame(xy = T) %>% 
  mutate(sesvert.dif = vert.mean.ses.future-vert.mean.ses) %>% 
  dplyr::select(x,y,sesvert.dif) %>% 
  rename(sesvert.dif.mammals = sesvert.dif)



load("results/birdselton_sar_sesvert2.RData")
birds = sesvert.scale.sub

birds.pdat = birds %>% 
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  #project("epsg:4326") %>% 
  as.data.frame(xy = T) %>% 
  mutate(sesvert.dif = vert.mean.ses.future-vert.mean.ses) %>% 
  dplyr::select(x,y,sesvert.dif) %>% 
  rename(sesvert.dif.birds = sesvert.dif)


load("results/reptiles_sar_sesvert2.RData")
repts = sesvert.scale.sub

repts.pdat = repts %>% 
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  #project("epsg:4326") %>% 
  as.data.frame(xy = T) %>% 
  mutate(sesvert.dif = vert.mean.ses.future-vert.mean.ses) %>% 
  dplyr::select(x,y,sesvert.dif) %>% 
  rename(sesvert.dif.repts = sesvert.dif)


# env data ----------------------------------------------------------------

env.future = list.files(path = "data/derivative_data/resampled_env_rasters/chelsa_future/2071_2100/ensemble/", pattern = ".tif", full.names = T)
env.future = lapply(env.future, rast)
env.future = rast(env.future)
names(env.future)[6] = "precip_ann"
env.future$temp_sea = env.future$temp_sea/100
env.future = as.data.frame(env.future, xy = T)
env.future = env.future %>% 
  rename_at(.vars = 3:11, .funs = function(x) paste0(x, "_future"))

env = read.csv("data/derivative_data/env_data.csv")

env = left_join(env.future, env, by = c('x', 'y'))

envdif = env %>% 
  mutate(tmin_cold_dif = min_temp_cold_future - tmin_cold,
         tmax_warm_dif = max_temp_warm_future - tmax_warm,
         precip_dry_dif = precip_dry_future - precip_dry,
         precip_wet_dif = precip_wet_future - precip_wet,
         precip_sea_dif = precip_sea_future - precip_sea) %>% 
  dplyr::select(x, y, tmin_cold_dif, tmax_warm_dif, precip_dry_dif, precip_wet_dif, precip_sea_dif, biome, realm)

difs = left_join(envdif, birds.pdat, by = c("x", "y")) %>% 
  left_join(mammals.pdat, by = c("x", "y")) %>% 
  left_join(repts.pdat, by = c("x", "y")) %>% 
  left_join(amph.pdat, by = c("x", "y")) %>% 
  pivot_longer(cols = 3:7, names_to = "climvar", values_to = "clim.dif") %>% 
  pivot_longer(cols = 5:8, names_to = "class", values_to = "sesvert.dif") %>% 
  drop_na()

ggplot(difs, aes(x = clim.dif, y = sesvert.dif)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(climvar ~ class, scales = "free", ncol = 4) +
  theme_classic()

ggplot(difs, aes(x = clim.dif, y = sesvert.dif, color = class)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  facet_wrap(biome~climvar, scale = "free", nrow = 5, strip.position = "left") +
  theme_classic()

plts = foreach(i = unique(difs$biome)) %do% {
 plt = ggplot(difs %>% filter(biome == i), aes(x = clim.dif, y = sesvert.dif, color = class)) +
    geom_point(size = 0.25, alpha = 0.5) +
    geom_smooth(method = "lm") +
    facet_wrap(~climvar, scale = "free", ncol = 5) +
    scale_color_discrete_qualitative() +
    ggtitle(i) +
    theme_classic()
}

for(i in 1:length(plts)) {print(plts[[i]])}

ggplot(difs, aes(x = sesvert.dif)) +
  geom_histogram() +
  geom_vline(xintercept = 0) +
  facet_grid(rows = vars(biome), cols = vars(class)) +
  theme_classic()

difs %>% 
  distinct(x, y, biome, realm, climvar, .keep_all = T) %>% 
  ggplot(aes(x = clim.dif)) +
  geom_histogram() +
  geom_vline(xintercept = 0) +
  facet_grid(rows = vars(biome), cols = vars(climvar), scales = "free") +
  theme_classic()





