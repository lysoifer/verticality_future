library(tidyverse)
library(png)


# function -----------------------------------------------------------
plt_pres_future = function(r, v, plot.title, vars, scale) {
  mp = ggplot() +
    geom_spatraster(data = r) +
    geom_sf(data = wd, fill = NA) +
    coord_sf(crs = 4326) +
    facet_wrap(~lyr, nrow = 2, strip.position = "left") +
    theme_classic() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "right",
          legend.margin = margin(0,0,0,0),
          legend.key = element_blank(),
          legend.background = element_blank(),
          #legend.key.width = unit(10, units = "mm"),
          legend.box.margin = margin(0,0,0,0),
          legend.box.spacing = unit(1,"mm"),
          legend.box.background = element_blank(),
          legend.title = element_blank(),
          strip.background = element_blank(),
          strip.text.y = element_blank())
  
  if(scale == "seq") {
    mp = mp + scale_fill_continuous_sequential(palette = "viridis", rev = T, na.value = "white") 
  } else {
    mp = mp + scale_fill_continuous_divergingx(palette = "spectral", rev = T, na.value = "white")
  }
  
    
    
  means = r %>% 
    as.data.frame(xy = T) %>% 
    dplyr::select(any_of(c("x", "y", vars))) %>% 
    pivot_longer(cols = c(.data[[vars[1]]], .data[[vars[2]]]), names_to = "Time", values_to = "val") %>% 
    group_by(Time) %>% 
    summarise(mean = mean(val, na.rm = T))
  
  
  hist.plt = r %>% 
    as.data.frame(xy = T) %>% 
    dplyr::select(any_of(c("x", "y", vars))) %>% 
    pivot_longer(cols = c(.data[[vars[1]]], .data[[vars[2]]]), names_to = "Time", values_to = "val") %>% 
    ggplot(aes(x = val)) +
    geom_histogram() +
    geom_vline(data = means, aes(xintercept = mean), color = "red3", linetype = "dashed", size = 1) +
    facet_wrap(facets = ~Time, nrow = 2, strip.position = "left") +
    scale_x_continuous(expand = c(0,0)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "cm"))

  plt.final = hist.plt + mp + plot_layout(widths = c(1,3)) + plot_annotation(title = plot.title)
  
  return(plt.final)
}


# plot ses vert mean ------------------------------------------------------

load("results/amphibians_sar_sesvert2.RData")
amph = sesvert.scale.sub

load("results/mammals_sar_sesvert2.RData")
mammals = sesvert.scale.sub

load("results/birdselton_sar_sesvert2.RData")
birds = sesvert.scale.sub

load("results/reptiles_sar_sesvert2.RData")
repts = sesvert.scale.sub

wd = vect("data/original/rnaturalearth_world.shp") %>% 
  project("epsg:4326") %>% 
  crop(ext(-180,180,-60,83.64))

amph.r = amph %>% 
  dplyr::select(x,y,vert.mean.ses, vert.mean.ses.future) %>%
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326")

amph.plt = plt_pres_future(amph.r, v = wd, plot.title = "Amphibians", vars = c("vert.mean.ses", "vert.mean.ses.future"))

repts.r = repts %>% 
  dplyr::select(x,y,vert.mean.ses, vert.mean.ses.future) %>%
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326")

repts.plt = plt_pres_future(repts.r, v = wd, plot.title = "Reptiles", vars = c("vert.mean.ses", "vert.mean.ses.future"))

mammals.r = mammals %>% 
  dplyr::select(x,y,vert.mean.ses, vert.mean.ses.future) %>%
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326")

mammals.plt = plt_pres_future(mammals.r, v = wd, plot.title = "Mammals", vars = c("vert.mean.ses", "vert.mean.ses.future"))

birds.r = birds %>% 
  dplyr::select(x,y,vert.mean.ses, vert.mean.ses.future) %>%
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326")

birds.plt = plt_pres_future(birds.r, v = wd, plot.title = "Birds", vars = c("vert.mean.ses", "vert.mean.ses.future"))





# plot vert mean ------------------------------------------------------

load("results/amphibians_sar_vertmean2.RData")
amph = vert.scale.sub

load("results/mammals_sar_vertmean2.RData")
mammals = vert.scale.sub

load("results/birdselton_sar_vertmean.RData")
birds = vert.scale.sub

load("results/reptiles_sar_vertmean2.RData")
repts = vert.scale.sub

wd = vect("data/original/rnaturalearth_world.shp") %>% 
  project("epsg:4326") %>% 
  crop(ext(-180,180,-60,83.64))

amph.r = amph %>% 
  dplyr::select(x,y,vert.mean, vert.mean.future) %>%
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326")

amph.plt = plt_pres_future(amph.r, v = wd, plot.title = "Amphibians", vars = c("vert.mean", "vert.mean.future"), scale = "seq")

repts.r = repts %>% 
  dplyr::select(x,y,vert.mean, vert.mean.future) %>%
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326")

repts.plt = plt_pres_future(repts.r, v = wd, plot.title = "Reptiles", vars = c("vert.mean", "vert.mean.future"), scale = "seq")

mammals.r = mammals %>% 
  dplyr::select(x,y,vert.mean, vert.mean.future) %>%
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326")

mammals.plt = plt_pres_future(mammals.r, v = wd, plot.title = "Mammals", vars = c("vert.mean", "vert.mean.future"), scale = "seq")

birds.r = birds %>% 
  dplyr::select(x,y,vert.mean, vert.mean.future) %>%
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326")

birds.plt = plt_pres_future(birds.r, v = wd, plot.title = "Birds", vars = c("vert.mean", "vert.mean.future"), scale = "seq")


amph.plt
repts.plt
mammals.plt
birds.plt


# plot ses vert mean FOREST ONLY ------------------------------------------------------

load("results/sar_mods_forest_only/amphibians/amphibians_sar_sesvert2.RData")
amph = sesvert.scale.sub

load("results/sar_mods_forest_only/mammals/mammals_sar_sesvert2.RData")
mammals = sesvert.scale.sub

load("results/sar_mods_forest_only/birds/birdselton_sar_sesvert2.RData")
birds = sesvert.scale.sub

load("results/sar_mods_forest_only/reptiles/reptiles_sar_sesvert2.RData")
repts = sesvert.scale.sub

wd = vect("data/original/rnaturalearth_world.shp") %>% 
  project("epsg:4326") %>% 
  crop(ext(-180,180,-60,83.64))

amph.r = amph %>% 
  dplyr::select(x,y,vert.mean.ses, vert.mean.ses.future) %>%
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326")

amph.plt = plt_pres_future(amph.r, v = wd, plot.title = "Amphibians", 
                           vars = c("vert.mean.ses", "vert.mean.ses.future"), scale = "div")

repts.r = repts %>% 
  dplyr::select(x,y,vert.mean.ses, vert.mean.ses.future) %>%
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326")

repts.plt = plt_pres_future(repts.r, v = wd, plot.title = "Reptiles", 
                            vars = c("vert.mean.ses", "vert.mean.ses.future"), scale = "div")

mammals.r = mammals %>% 
  dplyr::select(x,y,vert.mean.ses, vert.mean.ses.future) %>%
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326")

mammals.plt = plt_pres_future(mammals.r, v = wd, plot.title = "Mammals", 
                              vars = c("vert.mean.ses", "vert.mean.ses.future"), scale = "div")

birds.r = birds %>% 
  dplyr::select(x,y,vert.mean.ses, vert.mean.ses.future) %>%
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326")

birds.plt = plt_pres_future(birds.r, v = wd, plot.title = "Birds", 
                            vars = c("vert.mean.ses", "vert.mean.ses.future"), scale = "div")

png("figures/sesvert_future/forest_only/repts.png", width = 200, height = 100, res = 300, units = "mm")
repts.plt
dev.off()

png("figures/sesvert_future/forest_only/birds.png", width = 200, height = 100, res = 300, units = "mm")
birds.plt
dev.off()

png("figures/sesvert_future/forest_only/mammals.png", width = 200, height = 100, res = 300, units = "mm")
mammals.plt
dev.off()

png("figures/sesvert_future/forest_only/amphibians.png", width = 200, height = 100, res = 300, units = "mm")
amph.plt
dev.off()




# plot vert mean FOREST ONLY ------------------------------------------------------

load("results/sar_mods_forest_only/amphibians/amphibians_sar_vertmean2.RData")
amph = vert.scale.sub

load("results/sar_mods_forest_only/mammals/mammals_sar_vertmean2.RData")
mammals = vert.scale.sub

load("results/sar_mods_forest_only/birds/birdselton_sar_vertmean.RData")
birds = vert.scale.sub

load("results/sar_mods_forest_only/reptiles/reptiles_sar_vertmean2.RData")
repts = vert.scale.sub

wd = vect("data/original/rnaturalearth_world.shp") %>% 
  project("epsg:4326") %>% 
  crop(ext(-180,180,-60,83.64))

amph.r = amph %>% 
  dplyr::select(x,y,vert.mean, vert.mean.future) %>%
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326")

amph.plt = plt_pres_future(amph.r, v = wd, plot.title = "Amphibians", vars = c("vert.mean", "vert.mean.future"), scale = "seq")

repts.r = repts %>% 
  dplyr::select(x,y,vert.mean, vert.mean.future) %>%
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326")

repts.plt = plt_pres_future(repts.r, v = wd, plot.title = "Reptiles", vars = c("vert.mean", "vert.mean.future"), scale = "seq")

mammals.r = mammals %>% 
  dplyr::select(x,y,vert.mean, vert.mean.future) %>%
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326")

mammals.plt = plt_pres_future(mammals.r, v = wd, plot.title = "Mammals", vars = c("vert.mean", "vert.mean.future"), scale = "seq")

birds.r = birds %>% 
  dplyr::select(x,y,vert.mean, vert.mean.future) %>%
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326")

birds.plt = plt_pres_future(birds.r, v = wd, plot.title = "Birds", vars = c("vert.mean", "vert.mean.future"), scale = "seq")


amph.plt
repts.plt
mammals.plt
birds.plt


png("figures/vert_future/forest_only/repts.png", width = 200, height = 100, res = 300, units = "mm")
repts.plt
dev.off()

png("figures/vert_future/forest_only/birds.png", width = 200, height = 100, res = 300, units = "mm")
birds.plt
dev.off()

png("figures/vert_future/forest_only/mammals.png", width = 200, height = 100, res = 300, units = "mm")
mammals.plt
dev.off()

png("figures/vert_future/forest_only/amphibians.png", width = 200, height = 100, res = 300, units = "mm")
amph.plt
dev.off()
