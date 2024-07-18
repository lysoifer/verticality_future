library(tidyverse)
library(png)
library(ggplot2)
library(patchwork)
library(terra)
library(data.table)
library(tidyterra)
library(colorspace)


# function -----------------------------------------------------------
#plt_pres_future = function(r, v, plot.title, vars, scale) {
  
  mp.pres = ggplot() +
    geom_spatraster(data = r[[1]]) +
    geom_sf(data = wd, fill = NA) +
    coord_sf(crs = 4326) +
    theme_classic() +
    ggtitle(plot.title) +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
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
          strip.text.y = element_blank(),
          plot.background = element_rect(color = "black", fill = NA)) 
  
  if(scale == "seq") {
    mp.pres = mp.pres + scale_fill_continuous_sequential(palette = "viridis", rev = T, na.value = "white") 
  } else {
    mp.pres = mp.pres + scale_fill_continuous_divergingx(palette = "spectral", rev = T, na.value = "white")
  }
  
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
  
 
  return(list(mp.pres, plt.final))



# plot ses vert mean ------------------------------------------------------
# 
# load("results/sar_mods_forestOnly_forestSES/amphibians/amphibians_sar_sesvert.RData")
# amph = pred.df
# amph2 = fread("data/derivative_data/gridcell_data/amphibians_comdat/amph_comdat_parallel_forestsOnly.csv") %>% 
#   dplyr::select(x, y, p.arb, p.fos, p.ter, vert.mean, vert.mean.ses)
# 
# load("results/sar_mods_forestOnly_forestSES/mammals/mammals_sar_sesvert.RData")
# mammals = pred.df
# mammals2 = fread("data/derivative_data/gridcell_data/mammals_comdat/mammals_comdat_parallel_forestsOnly.csv") %>% 
#   dplyr::select(x, y, p.arb, p.fos, p.ter, vert.mean, vert.mean.ses)
# 
# load("results/sar_mods_forestOnly_forestSES/birds/birds_sar_sesvert.RData")
# birds = pred.df
# birds2 = fread("data/derivative_data/gridcell_data/birds_comdat/birds_comdat_parallel_elton_forestsOnly.csv") %>% 
#   dplyr::select(x, y, p.arb, p.fos, p.ter, vert.mean, vert.mean.ses)
# 
# 
# load("results/sar_mods_forestOnly_forestSES/reptiles/reptiles_sar_sesvert.RData")
# repts = pred.df
# repts2 = fread("data/derivative_data/gridcell_data/reptiles_comdat/rept_comdat_parallel_forestsOnly.csv") %>% 
#   dplyr::select(x, y, p.arb, p.fos, p.ter, vert.mean, vert.mean.ses)
# 
# 
# wd = vect("data/original/rnaturalearth_world.shp") %>% 
#   project("epsg:4326") %>% 
#   crop(ext(-180,180,-60,83.64))
# 
# amph.r = amph %>% 
#   #dplyr::select(x,y,vert.mean.ses, vert.mean.ses.future) %>%
#   rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
#   project("epsg:4326")
# 
# amph.plt = plt_pres_future(r = amph.r, v = wd, plot.title = "Amphibians", vars = c("vert.mean.ses", "vert.mean.ses.future"), scale = "div")
# amph.plt.vert = plt_vert_margin(amph.r, v = wd, prop = amph2, plot.title = "Amphibians", margin.title = "Proportion")
# 
# 
# repts.r = repts %>% 
#   dplyr::select(x,y,vert.mean.ses, vert.mean.ses.future) %>%
#   rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
#   project("epsg:4326")
# 
# repts.plt = plt_pres_future(repts.r, v = wd, plot.title = "Reptiles", vars = c("vert.mean.ses", "vert.mean.ses.future"), scale = "div")
# repts.plt.vert = plt_vert_margin(repts.r, v = wd, prop = repts2, plot.title = "Reptiles", margin.title = "Proportion")
# 
# mammals.r = mammals %>% 
#   dplyr::select(x,y,vert.mean.ses, vert.mean.ses.future) %>%
#   rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
#   project("epsg:4326")
# 
# mammals.plt = plt_pres_future(mammals.r, v = wd, plot.title = "Mammals", vars = c("vert.mean.ses", "vert.mean.ses.future"), scale = "div")
# mammals.plt.vert = plt_vert_margin(mammals.r, v = wd, prop = mammals2, plot.title = "Mammals", margin.title = "Proportion")
# 
# 
# birds.r = birds %>% 
#   dplyr::select(x,y,vert.mean.ses, vert.mean.ses.future) %>%
#   rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
#   project("epsg:4326")
# 
# birds.plt = plt_pres_future(birds.r, v = wd, plot.title = "Birds", vars = c("vert.mean.ses", "vert.mean.ses.future"), scale = "div")
# birds.plt.vert = plt_vert_margin(birds.r, v = wd, prop = birds2, plot.title = "Birds", margin.title = "Proportion")
# 
# 
# mp.pres = birds.plt[[1]] + mammals.plt[[1]] + repts.plt[[1]] + amph.plt[[1]] + plot_annotation(tag_levels = "A")
# png('figures/sesvert_maps/sesvert_vertebrates.png', width = 3000, height = 1500, res = 300)
# mp.pres
# dev.off()
# 
# # plot histograms separately from maps
# amph$class = "Amphibians"
# repts$class = "Reptiles"
# birds$class = "Birds"
# mammals$class = "Mammals"
# 
# verts = bind_rows(amph, repts, birds, mammals)
# 
# library(lemon)
# 
# verts.vline = verts %>% 
#   dplyr::select(vert.mean.ses, vert.mean.ses.future, class) %>% 
#   rename("Present SES verticality" = vert.mean.ses, "Future SES verticality" = vert.mean.ses.future) %>% 
#   pivot_longer(1:2, names_to = "time", values_to = "ses.vert") %>% 
#   group_by(class, time) %>% 
#   summarise(mean.sesvert = mean(ses.vert, na.rm = T)) %>% 
#   mutate(class = factor(class, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")),
#          time = factor(time, levels = c( "Future SES verticality", "Present SES verticality")))
# 
# sesvert.hists = verts %>% 
#   dplyr::select(vert.mean.ses, vert.mean.ses.future, class) %>% 
#   rename("Present SES verticality" = vert.mean.ses, "Future SES verticality" = vert.mean.ses.future) %>% 
#   pivot_longer(1:2, names_to = "time", values_to = "ses.vert") %>% 
#   mutate(class = factor(class, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")),
#          time = factor(time, levels = c( "Present SES verticality", "Future SES verticality"))) %>% 
#   ggplot(aes(x = ses.vert)) +
#   geom_histogram() +
#   geom_vline(data = verts.vline, aes(xintercept = mean.sesvert), color = "red3", linetype = "dashed", linewidth = 1) +
#   facet_rep_grid(rows = vars(time), cols = vars(class)) +
#   scale_x_continuous("SES mean verticality") +
#   theme_classic()
# 
# png("figures/sesvert_future_dif/sesvert_histograms.png", width = 2200, height = 1200, res = 300)
# sesvert.hists
# dev.off()
# 
# 
# # plot vert mean ------------------------------------------------------
# 
# load("results/amphibians_sar_vertmean2.RData")
# amph = vert.scale.sub
# 
# load("results/mammals_sar_vertmean2.RData")
# mammals = vert.scale.sub
# 
# load("results/birdselton_sar_vertmean.RData")
# birds = vert.scale.sub
# 
# load("results/reptiles_sar_vertmean2.RData")
# repts = vert.scale.sub
# 
# wd = vect("data/original/rnaturalearth_world.shp") %>% 
#   project("epsg:4326") %>% 
#   crop(ext(-180,180,-60,83.64))
# 
# amph.r = amph %>% 
#   dplyr::select(x,y,vert.mean, vert.mean.future) %>%
#   rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
#   project("epsg:4326")
# 
# amph.plt = plt_pres_future(amph.r, v = wd, plot.title = "Amphibians", vars = c("vert.mean", "vert.mean.future"), scale = "seq")
# 
# repts.r = repts %>% 
#   dplyr::select(x,y,vert.mean, vert.mean.future) %>%
#   rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
#   project("epsg:4326")
# 
# repts.plt = plt_pres_future(repts.r, v = wd, plot.title = "Reptiles", vars = c("vert.mean", "vert.mean.future"), scale = "seq")
# 
# mammals.r = mammals %>% 
#   dplyr::select(x,y,vert.mean, vert.mean.future) %>%
#   rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
#   project("epsg:4326")
# 
# mammals.plt = plt_pres_future(mammals.r, v = wd, plot.title = "Mammals", vars = c("vert.mean", "vert.mean.future"), scale = "seq")
# 
# birds.r = birds %>% 
#   dplyr::select(x,y,vert.mean, vert.mean.future) %>%
#   rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
#   project("epsg:4326")
# 
# birds.plt = plt_pres_future(birds.r, v = wd, plot.title = "Birds", vars = c("vert.mean", "vert.mean.future"), scale = "seq")
# 
# 
# amph.plt
# repts.plt
# mammals.plt
# birds.plt


# plot ses vert mean FOREST ONLY ------------------------------------------------------

# load("results/sar_mods_forest_only/amphibians/amphibians_sar_sesvert2.RData")
# amph = sesvert.scale.sub
# 
# load("results/sar_mods_forest_only/mammals/mammals_sar_sesvert2.RData")
# mammals = sesvert.scale.sub
# 
# load("results/sar_mods_forest_only/birds/birdselton_sar_sesvert2.RData")
# birds = sesvert.scale.sub
# 
# load("results/sar_mods_forest_only/reptiles/reptiles_sar_sesvert2.RData")
# repts = sesvert.scale.sub
# 
# wd = vect("data/original/rnaturalearth_world.shp") %>% 
#   project("epsg:4326") %>% 
#   crop(ext(-180,180,-60,83.64))
# 
# amph.r = amph %>% 
#   dplyr::select(x,y,vert.mean.ses, vert.mean.ses.future) %>%
#   rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
#   project("epsg:4326")
# 
# amph.plt = plt_pres_future(amph.r, v = wd, plot.title = "Amphibians", 
#                            vars = c("vert.mean.ses", "vert.mean.ses.future"), scale = "div")
# 
# repts.r = repts %>% 
#   dplyr::select(x,y,vert.mean.ses, vert.mean.ses.future) %>%
#   rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
#   project("epsg:4326")
# 
# repts.plt = plt_pres_future(repts.r, v = wd, plot.title = "Reptiles", 
#                             vars = c("vert.mean.ses", "vert.mean.ses.future"), scale = "div")
# 
# mammals.r = mammals %>% 
#   dplyr::select(x,y,vert.mean.ses, vert.mean.ses.future) %>%
#   rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
#   project("epsg:4326")
# 
# mammals.plt = plt_pres_future(mammals.r, v = wd, plot.title = "Mammals", 
#                               vars = c("vert.mean.ses", "vert.mean.ses.future"), scale = "div")
# 
# birds.r = birds %>% 
#   dplyr::select(x,y,vert.mean.ses, vert.mean.ses.future) %>%
#   rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
#   project("epsg:4326")
# 
# birds.plt = plt_pres_future(birds.r, v = wd, plot.title = "Birds", 
#                             vars = c("vert.mean.ses", "vert.mean.ses.future"), scale = "div")
# 
# png("figures/sesvert_future/forest_only/repts.png", width = 200, height = 100, res = 300, units = "mm")
# repts.plt[[2]]
# dev.off()
# 
# png("figures/sesvert_future/forest_only/birds.png", width = 200, height = 100, res = 300, units = "mm")
# birds.plt[[2]]
# dev.off()
# 
# png("figures/sesvert_future/forest_only/mammals.png", width = 200, height = 100, res = 300, units = "mm")
# mammals.plt[[2]]
# dev.off()
# 
# png("figures/sesvert_future/forest_only/amphibians.png", width = 200, height = 100, res = 300, units = "mm")
# amph.plt[[2]]
# dev.off()
# 
# mp.pres.forest = birds.plt[[1]] + mammals.plt[[1]] + repts.plt[[1]] + amph.plt[[1]] + plot_annotation(tag_levels = "A")
# png('figures/sesvert_maps/sesvert_vertebrates_forest.png', width = 3000, height = 1500, res = 300)
# mp.pres.forest
# dev.off()



# plot vert mean FOREST ONLY ------------------------------------------------------

# load("results/sar_mods_forest_only/amphibians/amphibians_sar_vertmean2.RData")
# amph = vert.scale.sub
# 
# load("results/sar_mods_forest_only/mammals/mammals_sar_vertmean2.RData")
# mammals = vert.scale.sub
# 
# load("results/sar_mods_forest_only/birds/birdselton_sar_vertmean.RData")
# birds = vert.scale.sub
# 
# load("results/sar_mods_forest_only/reptiles/reptiles_sar_vertmean2.RData")
# repts = vert.scale.sub

wd = vect("data/original/rnaturalearth_world.shp") %>% 
  project("epsg:4326") %>% 
  crop(ext(-180,180,-60,83.64))

# amph.r = amph %>% 
#   dplyr::select(x,y,vert.mean, vert.mean.future) %>%
#   rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
#   project("epsg:4326")
# 
# amph.plt = plt_pres_future(amph.r, v = wd, plot.title = "Amphibians", vars = c("vert.mean", "vert.mean.future"), scale = "seq")
# 
# repts.r = repts %>% 
#   dplyr::select(x,y,vert.mean, vert.mean.future) %>%
#   rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
#   project("epsg:4326")
# 
# repts.plt = plt_pres_future(repts.r, v = wd, plot.title = "Reptiles", vars = c("vert.mean", "vert.mean.future"), scale = "seq")
# 
# mammals.r = mammals %>% 
#   dplyr::select(x,y,vert.mean, vert.mean.future) %>%
#   rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
#   project("epsg:4326")
# 
# mammals.plt = plt_pres_future(mammals.r, v = wd, plot.title = "Mammals", vars = c("vert.mean", "vert.mean.future"), scale = "seq")
# 
# birds.r = birds %>% 
#   dplyr::select(x,y,vert.mean, vert.mean.future) %>%
#   rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
#   project("epsg:4326")
# 
# birds.plt = plt_pres_future(birds.r, v = wd, plot.title = "Birds", vars = c("vert.mean", "vert.mean.future"), scale = "seq")
# 
# 
# amph.plt
# repts.plt
# mammals.plt
# birds.plt
# 
# 
# png("figures/vert_future/forest_only/repts.png", width = 200, height = 100, res = 300, units = "mm")
# repts.plt
# dev.off()
# 
# png("figures/vert_future/forest_only/birds.png", width = 200, height = 100, res = 300, units = "mm")
# birds.plt
# dev.off()
# 
# png("figures/vert_future/forest_only/mammals.png", width = 200, height = 100, res = 300, units = "mm")
# mammals.plt
# dev.off()
# 
# png("figures/vert_future/forest_only/amphibians.png", width = 200, height = 100, res = 300, units = "mm")
# amph.plt
# dev.off()



# add proportion plots ----------------------------------------------------

plt_vert_margin = function(v, df, plot.title, tags = list(c("A", "B", "C"))) {
  # r: spatraster for plotting
  # v: spatvector for outline
  # plot.title = plot title
  # margin.titl = y-axis margin title
  
  # get minmax values from vector
  ymin = terra::ext(v)[3]
  ymax = terra::ext(v)[4]
  xmax = terra::ext(v)[2]
  
  # get color scale limits
  # col.lims = round(minmax(r),2)
  col.lims = c(round(min(df$vert.mean.ses),2), round(max(df$vert.mean.ses),2))
  
  r = df %>% rast(crs = "+proj=cea +datum=WGS84") %>% 
    project("epsg:4326") 
  
    # create main plot
  plot_main = ggplot() +
    # Overlay world
    geom_spatvector(data = wd, color = alpha("black", 0.7), fill = "gray70", linewidth = .1) +
    geom_spatraster(data = r$vert.mean.ses, maxcell = Inf) +
    geom_spatvector(data = wd, color = "black", fill = NA, linewidth = .1) +
    scale_fill_continuous_divergingx("spectral", na.value = NA, limits = col.lims, rev = T) +
    guides(fill = guide_colorbar(title = "")) +
   # theme_classic() +
    theme(legend.position = "inside",
          legend.position.inside = c(0.6,0.45),
          legend.justification = c(0,1),
          legend.margin = margin(0,0,0.2,0.3),
          legend.key.height = unit(3, units = "mm"),
          legend.background = element_blank(),
          legend.key.width = unit(2, units = "mm"),
          legend.box.margin = margin(0,0,0,0),
          legend.box.spacing = unit(0,"mm"),
          legend.box.background = element_blank(),
          plot.background = element_blank(),
          plot.margin = unit(c(0,0,0,0), units = "mm"),
          panel.spacing = unit(0,"mm"),
          #panel.background = element_rect(color = "black", fill = NA),
          panel.background = element_blank(),
          #legend.title = element_blank(),
          #legend.title.position = "top",
          legend.title = element_text(hjust = 0.5, size = 8),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank())

  
  plot_main2 = ggplot() +
    # Overlay world
    geom_spatvector(data = wd, fill = "gray70", color = "black", linewidth = .1) +
    geom_spatraster(data = r$vert.mean, maxcell = Inf) +
    scale_fill_continuous_divergingx("spectral", na.value = NA, limits = c(0,1), mid = 0.5, rev = T, breaks = c(0,0.5,1)) +
    geom_spatvector(data = wd, fill = NA, color = "black", linewidth = .1) +
    guides(fill = guide_colorbar(title = "")) +
    # theme_classic() +
    theme(legend.position = "inside",
          legend.position.inside = c(0.6,0.45),
          legend.justification = c(0,1),
          legend.margin = margin(0,0,0.2,0.3),
          #legend.key = element_blank(),
          legend.background = element_blank(),
          legend.key.height = unit(3, units = "mm"),
          legend.key.width = unit(2, units = "mm"),
          legend.box.margin = margin(0,0,0,0),
          legend.box.spacing = unit(1,"mm"),
          legend.box.background = element_blank(),
          plot.background = element_blank(),
          plot.margin = unit(c(0,0,0,0), units = "mm"),
          panel.spacing = unit(0,"mm"),
          #panel.background = element_rect(color = "black", fill = NA),
          panel.background = element_blank(),
          #legend.title.position = "top",
          #legend.title = element_text(hjust = 0.5, size = 8),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank())

  prop = r %>% 
    as.data.frame(xy = T) %>% 
    dplyr::select(x,y,p.arb, p.fos, p.ter) %>% 
    rename(Arboreal = p.arb, Fossorial = p.fos, Terrestrial = p.ter) %>% 
    pivot_longer(3:5, names_to = "micro", values_to = "prop") %>% 
    group_by(y, micro) %>% 
    summarize(prop.mean = mean(prop), prop.sd = sd(prop)) %>% 
    mutate(micro = factor(micro, levels = c("Arboreal", "Terrestrial", "Fossorial")))
  
  # testing
  prop.plt = prop %>% 
    ggplot() +
    geom_line(aes(x = y, y = prop.mean, color = micro), size=1) +
    geom_ribbon(aes(x = y, ymin = prop.mean - prop.sd, ymax = prop.mean + prop.sd, fill = micro), alpha = 0.2) +
    scale_y_continuous("", breaks = seq(0,1,0.5)) +
    scale_x_continuous("Latitude", breaks = seq(-30,60,30), limits = c(-55.61183,83.64)) +
    scale_color_discrete_divergingx("zissou1", guide = guide_legend(title = ""), rev = T) +
    scale_fill_discrete_divergingx("zissou1", guide = guide_legend(title = ""), rev = T) +
    coord_flip() +
    theme_classic() +
    theme(axis.title.x = element_blank())
  
  
  row_name = wrap_elements(panel = textGrob(plot.title, rot = 90))

  p = row_name + plot_main + plot_main2 + prop.plt + plot_layout(nrow = 1, width = unit(c(0.4, -1, -1, 3.5), c("cm", "null", "null", "cm"))) +
    plot_annotation(tag_levels = tags, theme = theme(plot.title = element_text(hjust = 0.5)))
  p[[2]] <- p[[2]] + theme(legend.title = element_blank(), legend.position.inside = c(0.1,0.5))
  p[[3]] <- p[[3]] + theme(legend.title = element_blank(), legend.position.inside = c(0.1,0.5))

  return(p)
}


# plot forest only with proportion margin

wd = vect("data/original/rnaturalearth_world.shp") %>% 
  project("epsg:4326") %>% 
  crop(ext(-180,180,-60,83.64))

amph = fread("data/derivative_data/gridcell_data/amphibians_comdat/amph_comdat_parallel_forestsOnly.csv") %>% 
  dplyr::select(x, y, p.arb, p.fos, p.ter, vert.mean, vert.mean.ses, biome)

vert.amph.forest = plt_vert_margin(v = wd, df = amph, plot.title = "Amphibians", tags = list(c("", "J", "K", "L")))

png("figures/vert_maps_forest/vert_present/amphibians.png", height = 800, width = 2500, res = 300)
vert.amph.forest
dev.off()

mammals = fread("data/derivative_data/gridcell_data/mammals_comdat/mammals_comdat_parallel_forestsOnly.csv") %>% 
  dplyr::select(x, y, p.arb, p.fos, p.ter, vert.mean, vert.mean.ses, biome)

vert.mammals.forest = plt_vert_margin(v = wd, df = mammals, plot.title = "Mammals", tags = list(c("D", "E", "F")))

png("figures/vert_maps_forest/vert_present/mammals.png", height = 800, width = 2500, res = 300)
vert.mammals.forest
dev.off()


birds = fread("data/derivative_data/gridcell_data/birds_comdat/birds_comdat_parallel_elton_forestsOnly.csv") %>% 
  dplyr::select(x, y, p.arb, p.fos, p.ter, vert.mean, vert.mean.ses, biome)

vert.birds.forest = plt_vert_margin(v = wd, df = birds, plot.title = "Birds", tags = list(c("A", "B", "C")))

png("figures/vert_maps_forest/vert_present/birds.png", height = 800, width = 2500, res = 300)
vert.birds.forest
dev.off()


repts = fread("data/derivative_data/gridcell_data/reptiles_comdat/rept_comdat_parallel_forestsOnly.csv") %>% 
  dplyr::select(x, y, p.arb, p.fos, p.ter, vert.mean, vert.mean.ses, biome)

vert.repts.forest = plt_vert_margin(v = wd, df = repts, plot.title = "Reptiles", tags = list(c("G", "H", "I")))

png("figures/vert_maps_forest/vert_present/reptiles.png", height = 800, width = 2500, res = 300)
vert.repts.forest
dev.off()

p = vert.birds.forest /
  vert.mammals.forest / 
  vert.repts.forest / 
  vert.amph.forest
p[[1]][[2]] <- p[[1]][[2]] + ggtitle("SES Verticality") + theme(plot.title = element_text(hjust = 0.5))
p[[1]][[3]] <- p[[1]][[3]] + ggtitle("Mean Verticality") + theme(plot.title = element_text(hjust = 0.5))
p[[4]][[4]] <- p[[4]][[4]] + scale_y_continuous("Proportion", breaks = seq(0,1,0.5)) + theme(axis.title.x = element_text())


png("figures/main_figs/vert_maps.png", height = 300, width = 500, res = 300, units = "mm")
p
dev.off()
