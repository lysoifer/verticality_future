library(terra)
library(tidyverse)
library(tidyterra)
library(data.table)

env = fread("data/derivative_data/env_data.csv")

ggplot(data = env, aes(x = x, y = y, fill = biome)) +
  geom_raster() +
  coord_sf(crs = "+proj=cea +datum=WGS84")

