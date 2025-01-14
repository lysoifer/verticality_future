library(tidyverse)
library(png)
library(ggplot2)
library(patchwork)
library(terra)
library(data.table)
library(tidyterra)
library(colorspace)
library(grid)


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
          #legend.key.height = unit(3, units = "mm"),
          legend.key.height = unit(2, units = "mm"),
          legend.background = element_blank(),
          legend.key.width = unit(2, units = "mm"),
          legend.box.margin = margin(0,0,0,0),
          legend.box.spacing = unit(0,"mm"),
          legend.box.background = element_blank(),
          legend.text = element_text(size = 6),
          plot.background = element_blank(),
          plot.margin = unit(c(0,0,0,0), units = "mm"),
          panel.spacing = unit(0,"mm"),
          #panel.background = element_rect(color = "black", fill = NA),
          panel.background = element_rect(fill = "transparent", color = NA),
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
          #legend.key.height = unit(3, units = "mm"),
          legend.key.height = unit(2, units = "mm"),
          legend.key.width = unit(2, units = "mm"),
          legend.box.margin = margin(0,0,0,0),
          legend.box.spacing = unit(1,"mm"),
          legend.box.background = element_blank(),
          legend.text = element_text(size = 6),
          plot.background = element_rect(fill = "transparent", color = NA),
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
  
  
  row_name = wrap_elements(panel = textGrob(plot.title, rot = 90, gp=gpar(fontsize = 10)))

  p = row_name + plot_main + plot_main2 + prop.plt + plot_layout(nrow = 1, width = unit(c(0.4, -1, -1, 2), c("cm", "null", "null", "cm"))) +
    plot_annotation(tag_levels = tags, theme = theme(plot.title = element_text(hjust = 0.5)))
  p[[2]] <- p[[2]] + theme(legend.title = element_blank(), legend.position.inside = c(0.1,0.5))
  p[[3]] <- p[[3]] + theme(legend.title = element_blank(), legend.position.inside = c(0.1,0.5))
  p = p & theme(panel.background = element_rect(fill = NA, color = NA))

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
p[[1]][[2]] <- p[[1]][[2]] + ggtitle("SES Verticality") + theme(plot.title = element_text(hjust = 0.5, size = 10))
p[[1]][[3]] <- p[[1]][[3]] + ggtitle("Mean Verticality") + theme(plot.title = element_text(hjust = 0.5, size = 10))
p[[4]][[4]] <- p[[4]][[4]] + scale_y_continuous("Proportion", breaks = seq(0,1,0.5)) + theme(axis.title.x = element_text())


png("figures/main_figs/vert_maps.png", height = 140, width = 180, res = 300, units = "mm")
p
dev.off()
