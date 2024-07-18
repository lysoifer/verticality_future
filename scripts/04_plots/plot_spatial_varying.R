library(terra)
library(ggplot2)
library(colorspace)
library(tidyverse)
library(patchwork)
library(cowplot)

wd = vect("data/original/rnaturalearth_world.shp") %>% 
  project("epsg:4326") %>% 
  crop(ext(-180,180,-60,83.64))

plot_spatial_varying = function(mod, var, v, legend_title) {
  coefs = tidy(mod)
  coef.var = as.numeric(coefs[which(coefs$term == var), "estimate"])
  pred = predict(mod, type = "response")
  dat = pred %>% 
    mutate(x = x*1e5, y = y*1e5,
           canopy_height_est_slope = zeta_s_canopy_height + coef.var) %>% 
    relocate(x, y, .before = rich) %>% 
    rast(crs = "+proj=cea + datum=WGS84") %>% 
    project("epsg:4326")
    
    
  ggplot() +
    geom_spatvector(data = v, color = "black", fill = "gray70") +
    geom_spatraster(data = dat, aes(fill = canopy_height_est_slope)) +
    geom_spatvector(data = v, color = "black", fill = NA) +
    scale_fill_continuous_divergingx("BrBG", na.value = NA, guide = guide_colorbar(legend_title)) +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(color = "black", fill = NA),
          panel.grid = element_blank())

}

# plot canopy height coefficienby canopy height standard deviation
plot_coef_height = function(pred) {
  height = pred %>% 
    dplyr::select(x,y,canopy_height) %>% 
    mutate(x = x*1e5, y = y*1e5) %>% 
    rast(crs = "+proj=cea +datum=WGS84") %>% 
    focal(w = 5, fun = sd, na.rm = T) %>% 
    as.data.frame(xy = T) %>% 
    mutate(x = round(x/1e5,5), y = round(y/1e5,5)) %>% 
    inner_join(pred, by = c("x", "y")) %>% 
    dplyr::select(x,y,biome, ecoregion, biorealm, zeta_s_canopy_height, focal_sd) %>% 
    rename(canopy_height_sd = focal_sd)
  
  p = height %>% 
    ggplot(aes(x = canopy_height_sd, y = zeta_s_canopy_height)) +
    geom_point(pch = ".") +
    geom_smooth(method = "lm") +
    scale_x_continuous("SD canopy height", breaks = seq(0,1.25,0.25)) +
    scale_y_continuous("Canopy height coefficient") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p)
    
}



# Amphibians

load("results/sdmTMB_models/amphibians_sesvert.RData")
mod = m.final
amph.sesvert = plot_spatial_varying(mod = m.final, var = "canopy_height", v = wd, legend_title = "coefficient")
coef.a.ses = plot_coef_height(pred)

load("results/sdmTMB_models/amphibians_meanvert.RData")
mod = m.final
amph.meanvert = plot_spatial_varying(mod = m.final, var = "canopy_height", v = wd, legend_title = "coefficient")
coef.a.mean = plot_coef_height(pred)


# Reptiles

load("results/sdmTMB_models/reptiles_sesvert.RData")
mod = m.final
rept.sesvert = plot_spatial_varying(mod = m.final, var = "canopy_height", v = wd, legend_title = "coefficient")
coef.r.ses = plot_coef_height(pred)

load("results/sdmTMB_models/reptiles_meanvert.RData")
mod = m.final
rept.meanvert = plot_spatial_varying(mod = m.final, var = "canopy_height", v = wd, legend_title = "coefficient")
coef.r.mean = plot_coef_height(pred)


# Mammals

load("results/sdmTMB_models/mammals_sesvert.RData")
mod = m.final
mammals.sesvert = plot_spatial_varying(mod = m.final, var = "canopy_height", v = wd, legend_title = "coefficient")
coef.m.ses = plot_coef_height(pred)

load("results/sdmTMB_models/mammals_meanvert.RData")
mod = m.final
mammals.meanvert = plot_spatial_varying(mod = m.final, var = "canopy_height", v = wd, legend_title = "coefficient")
coef.m.mean = plot_coef_height(pred)


# Birds

load("results/sdmTMB_models/birds_sesvert.RData")
mod = m.final
birds.sesvert = plot_spatial_varying(mod = m.final, var = "canopy_height", v = wd, legend_title = "coefficient")
coef.b.ses = plot_coef_height(pred)

load("results/sdmTMB_models/birds_meanvert.RData")
mod = m.final
birds.meanvert = plot_spatial_varying(mod = m.final, var = "canopy_height", v = wd, legend_title = "coefficient")
coef.b.mean = plot_coef_height(pred)


# Full figure

design = "
AF
BG
CH
DI
EJ"

library(grid)
col1 = wrap_elements(panel = textGrob("SES Mean Verticality"))
col2 = wrap_elements(panel = textGrob("Mean Verticality"))

p = col1 + birds.sesvert+ mammals.sesvert + rept.sesvert + amph.sesvert +
  col2 + birds.meanvert + mammals.meanvert + rept.meanvert + amph.meanvert + 
  plot_layout(design = design, heights = c(0.01, -1, -1, -1, -1)) +
  plot_annotation(tag_levels = list(c("", "A", "B", "C", "D", "", "E", "F", "G", "H")))

png("figures/supp_figs/canopy_height_coefs.png", height = 300, width = 300, res = 300, units = "mm")
p
dev.off()

row1 = wrap_elements(panel = textGrob("Birds", rot = 90))
row2 = wrap_elements(panel = textGrob("Mammals", rot = 90))
row3 = wrap_elements(panel = textGrob("Reptiles", rot = 90))
row4 = wrap_elements(panel = textGrob("Amphibians", rot = 90))

design = "
AFK
BGL
CHM
DIN
EJO"

p = 
  plot_spacer() + row1 + row2 + row3 + row4 +
  col1 + coef.b.ses + coef.m.ses + coef.r.ses + coef.a.ses +
  col2 + coef.b.mean + coef.m.mean + coef.r.mean + coef.a.mean + 
  plot_layout(design = design, heights = c(0.2,1,1,1,1), widths = c(0.2,1,1)) +
  plot_annotation(tag_levels = list(c("","","","","", "A", "B", "C", "D", "", "E", "F", "G", "H")))

p[[12]] <- p[[12]] + theme(axis.title = element_blank())
p[[13]] <- p[[13]] + theme(axis.title = element_blank())
p[[14]] <- p[[14]] + theme(axis.title = element_blank())
p[[15]] <- p[[15]] + theme(axis.title.y = element_blank())
p[[7]] <- p[[7]] + theme(axis.title.x = element_blank())
p[[8]] <- p[[8]] + theme(axis.title.x = element_blank())
p[[9]] <- p[[9]] + theme(axis.title.x = element_blank())

png("figures/supp_figs/canopyHeightCoefs-canopyHeightSD.png", height = 300, width = 200, res = 300, units = "mm")
p
dev.off()
