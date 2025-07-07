library(tidyverse)
library(png)
library(ggplot2)
library(patchwork)
library(terra)
library(data.table)
library(tidyterra)
library(colorspace)
library(grid)
#library(tricolore)
library(lemon)


# load data ----------------------------------------------------

wd = vect("data/original/rnaturalearth_world.shp") %>% 
  crop(ext(-180,180,-60,83.64)) %>% 
  project("+proj=robin +datum=WGS84") 

wd = vect("data/original/rnaturalearth_world.shp") %>% 
  crop(ext(-180,180,-60,83.64)) %>% 
  project("+proj=cea +datum=WGS84") 

ymin = ext(wd)[3]
ymax = ext(wd)[4]

amph = fread("data/derivative_data/gridcell_data/env_forest/50_km/amph_comdat.csv") %>% 
  dplyr::select(x, y, rich, p.arb, p.fos, p.ter, vert.mean, vert.mean.ses, biome) %>% 
  filter(rich >= 5) %>% 
  filter(!is.na(vert.mean.ses)) %>% 
  drop_na() %>% 
  rast(crs = "+proj=cea +datum=WGS84") %>% 
  project("+proj=robin +datum=WGS84") %>% 
  as.data.frame(xy = T) %>% 
  mutate(taxa = "Amphibians")

mammals = fread("data/derivative_data/gridcell_data/env_forest/50_km/mammals_comdat.csv") %>% 
  dplyr::select(x, y, rich, p.arb, p.fos, p.ter, vert.mean, vert.mean.ses, biome) %>% 
  filter(rich >= 5) %>% 
  filter(!is.na(vert.mean.ses)) %>% 
  drop_na() %>% 
  rast(crs = "+proj=cea +datum=WGS84") %>% 
  project("+proj=robin +datum=WGS84") %>% 
  as.data.frame(xy = T) %>% 
  mutate(taxa = "Mammals")

birds = fread("data/derivative_data/gridcell_data/env_forest/50_km/birds_comdat.csv") %>% 
  dplyr::select(x, y, rich, p.arb, p.fos, p.ter, vert.mean, vert.mean.ses, biome) %>% 
  filter(rich >= 5) %>% 
  filter(!is.na(vert.mean.ses)) %>% 
  drop_na()%>% 
  rast(crs = "+proj=cea +datum=WGS84") %>% 
  project("+proj=robin +datum=WGS84") %>% 
  as.data.frame(xy = T) %>%  
  mutate(taxa = "Birds")

repts = fread("data/derivative_data/gridcell_data/env_forest/50_km/reptiles_comdat.csv") %>% 
  dplyr::select(x, y, rich, p.arb, p.fos, p.ter, vert.mean, vert.mean.ses, biome) %>% 
  filter(rich >= 5) %>% 
  filter(!is.na(vert.mean.ses)) %>% 
  drop_na()%>% 
  rast(crs = "+proj=cea +datum=WGS84") %>% 
  project("+proj=robin +datum=WGS84") %>% 
  as.data.frame(xy = T) %>%  
  mutate(taxa = "Reptiles")

all = bind_rows(amph, mammals, birds, repts) %>% 
  drop_na(vert.mean.ses) %>% 
  mutate(taxa = factor(taxa, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")))
min = min(all$vert.mean.ses, na.rm = T)
max = max(all$vert.mean.ses, na.rm = T)
col.lims = range(all$vert.mean.ses, na.rm = T)



# PLOT VERTICALITY --------------------------------------------------------

vertmaps = list()
for(i in levels(all$taxa)) {
  
  d = all %>% filter(taxa == i)
  
  map = ggplot() +
    geom_spatvector(data = wd, fill = "black", color = NA) +
    geom_raster(data = d, aes(x,y,fill = vert.mean.ses)) +
    #geom_spatvector(data = wd, fill = NA, color = "black") +
    scale_x_continuous("SES Verticality", expand = c(0,0)) +
    scale_fill_continuous_divergingx("spectral", na.value = NA, limits = col.lims,
                                     rev = T, guide = guide_colorbar(title = "Verticality")) +
    #facet_wrap(vars(taxa), ncol = 1, strip.position = "left") +
    coord_sf(crs = "+proj=robin +datum=WGS84") +
    theme_classic() +
    theme(strip.text = element_blank(),
          legend.position = "bottom",
          legend.title.position = "top",
          legend.title = element_text(hjust = 0.5, size = 10),
          legend.key.height = unit(3, units = "mm"),
          legend.key.width = unit(10, units = "mm"),
          axis.title = element_blank(),
          axis.line = element_line(color = "white"),
          axis.text = element_blank(),
          axis.ticks = element_line(color = "white"))
  
  vert.legend = ggpubr::get_legend(map)
  #vert.legend = as_ggplot(vert.legend)
  vertmaps[[i]] = map + theme(legend.position = "none")
  
  
}


# PLOT VERT-RICHNESS CORRELATIONS -------------------------------------------------------

corplot = list()

for(i in levels(all$taxa)) {
  d = all %>% filter(taxa==i)
  
  p = ggplot(d, aes(vert.mean.ses, rich)) +
    geom_point(pch = ".", alpha = 0.05) +
    geom_smooth(method = "lm", color = "green3") +
    #facet_rep_wrap(~taxa, scales = "free_y", ncol = 1) +
    scale_y_continuous("Richness") +
    scale_x_continuous("Verticality") +
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          axis.title = element_blank())
  
  corplot[[i]] = p
  
}


# add r2 values
# try flipping y-axis to right side of the plot


# PLOT PROPORTIONS ACROSS LATITUDES ---------------------------------------

pal = divergingx_hcl(palette = "zissou", n = 3)

cola = pal[1]
colt = pal[2]
colf = pal[3]


propplots = list()

for(i in levels(all$taxa)) {
  p = all %>% 
    filter(taxa == i) %>% 
    dplyr::select(p.arb, p.fos, p.ter, taxa, y) %>% 
    rename(Arboreal = p.arb, Fossorial = p.fos, Terrestrial = p.ter) %>% 
    pivot_longer(1:3, names_to = "micro", values_to = "prop") %>% 
    group_by(y, micro, taxa) %>% 
    summarize(prop.mean = mean(prop), prop.sd = sd(prop)) %>% 
    mutate(micro = factor(micro, levels = c("Arboreal", "Terrestrial", "Fossorial"))) %>% 
    #plot
    ggplot() +
    geom_line(aes(x = y, y = prop.mean, color = micro), size=1) +
    geom_ribbon(aes(x = y, ymin = prop.mean - prop.sd, ymax = prop.mean + prop.sd, fill = micro), alpha = 0.2) +
    scale_y_continuous("Proportion", breaks = seq(0,1,0.5)) +
    # scale_x_continuous("Latitude", breaks = seq(-60,60,30)) +
    #scale_color_discrete_divergingx("zissou1", guide = guide_legend(title = ""), rev = T) +
    scale_color_manual(values = c(cola, colt, colf), guide = guide_legend(title = "")) +
    scale_fill_manual(values = c(cola, colt, colf), guide = guide_legend(title = "")) +
    #scale_fill_discrete_divergingx("zissou1", guide = guide_legend(title = ""), rev = T) +
    #facet_rep_wrap(vars(taxa), ncol = 1, strip.position = "left") +
    coord_flip(xlim = c(ymin,ymax), clip = "off") +
    theme(axis.title.y = element_blank(),
          legend.position = "none",
          strip.placement = "outside",
          strip.background = element_rect(color = NA, fill = NA),
          axis.line = element_line(color = "black"),
          panel.background = element_blank(),
          strip.text = element_text(size = 10),
          axis.title = element_blank())
  
  propplots[[i]] = p
  
}



# Combine plots -----------------------------------------------------------

library(gridExtra)
library(grid)
library(gtable)

tight_theme <- theme(
  plot.margin = margin(0, 0, 0, 0),
  panel.spacing = unit(0, "pt"),  # for facet_wrap/grid
  #axis.ticks.length = unit(1, "pt"),
  #strip.placement = "outside"
)

# Apply to all plots
propplots <- lapply(propplots, `+`, tight_theme)
vertmaps  <- lapply(vertmaps,  `+`, tight_theme)
corplot  <- lapply(corplot,  `+`, tight_theme)

rows = list()

for (i in levels(all$taxa)) {
  # Convert each ggplot to grobs
  g_prop = ggplotGrob(propplots[[i]])
  g_map  = ggplotGrob(vertmaps[[i]])
  g_cor  = ggplotGrob(corplot[[i]])
  
  # Match heights to avoid vertical misalignment
  max_height = grid::unit.pmax(g_prop$heights, g_map$heights, g_cor$heights)
  g_prop$heights <- g_map$heights <- g_cor$heights <- max_height
  g_map$widths <- g_map$widths*3
  
  # Combine into one row (left to right: prop | map | cor)
  row = gtable_cbind(g_prop, g_map, g_cor, size = "first")
  rows[[i]] = row
}

# Combine all rows vertically
final_gtable = do.call(gtable_rbind, rows)

# Draw the full aligned plot
grid.newpage()
grid.draw(final_gtable)



# BIVARIATE MAPS ----------------------------------------------------------

library(biscale)

for(i in levels(all$taxa)) {
  d = all %>% filter(taxa == i)
  
  bclass = bi_class(d, x = vert.mean.ses, y = rich,
                    style = "quantile", dim = 4)
  
 map =  ggplot() +
    geom_raster(data = bclass, aes(x,y,fill = bi_class)) +
    bi_scale_fill(pal = "PinkGrn", dim = 4) +
    coord_quickmap() +
    bi_theme(base_size = 16) +
    theme(legend.position = "none")
  
  legend = bi_legend(pal = "PinkGrn", dim = 4,
                     xlab = "Verticality",
                     ylab = "Richness",
                     size = 6)
  
  finalPlot <- ggdraw() +
    draw_plot(map, 0, 0, 1, 1) +  # Draw the main map plot
    draw_plot(legend, 0.05, 0.05, 0.28, 0.28)

}













