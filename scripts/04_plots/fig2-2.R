library(tidyverse)
library(terra)
library(foreach)
library(colorspace)
library(tidyterra)
library(MASS)
library(ggpubr)
library(patchwork)
library(grid)


thm = theme(legend.key.height = unit(3, units = "mm"),
            legend.key.width = unit(10, units = "mm"),
            plot.background = element_blank(),
            plot.margin = unit(c(0,0,0,0), units = "mm"),
            panel.spacing = unit(0,"mm"),
            panel.background = element_rect(fill = "transparent", color = NA),
            strip.background = element_rect(fill = "white", color = "black"),
            panel.grid = element_blank(),
            legend.title = element_text(hjust = 0.5, size = 10))

wd = vect("data/original/rnaturalearth_world.shp") %>% 
  project("epsg:4326") %>% 
  crop(ext(-180,180,-60,83.64))

wd.cea = project(wd, "+proj=cea +datum=WGS84")

amph = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/amph_comdat.csv")
birds = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/birds_comdat.csv")
rept = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/reptiles_comdat.csv")
mammals = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/mammals_comdat.csv")

amph$taxa = "Amphibians"
birds$taxa = "Birds"
rept$taxa = "Reptiles"
mammals$taxa = "Mammals"

df = rbind(amph, birds, rept, mammals)
df = df %>% 
  mutate(biome2 = case_when(biome == "Boreal Forests/Taiga" ~ "Temperate & arctic",
                            biome == "Temperate Conifer Forests" ~ "Temperate & arctic",
                            biome == "Temperate Broadleaf & Mixed Forests" ~ "Temperate & arctic",
                            biome == "Tropical & Subtropical Moist Broadleaf Forests" ~ "Tropical & subtropical moist",
                            biome == "Tropical & Subtropical Coniferous Forests" ~ "Tropical & subtropical conifer",
                            biome == "Tropical & Subtropical Dry Broadleaf Forests" ~ "Tropical & subtropical dry",
                            biome == "Tropical & Subtropical Grasslands, Savannas & Shrublands" ~ "Tropical & subtropical savannas & shrublands",
                            biome == "Tundra" ~ "Temperate & arctic",
                            .default = "Temperate savanna, shrub, & scrub"),
         biome3 = ifelse(grepl("Tropical", biome2), "Tropical", "Temperate/arctic")) %>% 
  # drop_na(biome2) %>% 
  filter(rich > 5)


# GLOBAL RICHNESS ---------------------------------------------------------

df_combined = df %>% 
  dplyr::select(x,y,rich, canopy_height, taxa, biome2) %>% 
  pivot_wider(names_from = taxa, values_from = rich, names_prefix = "richness_") %>% 
  mutate(log_ch = log(canopy_height)) %>% 
  filter(canopy_height > 0)
df_combined$total_richness = apply(df_combined[,5:8], 1, sum, na.rm = T)

r_combined = df_combined %>% 
  rast(crs = "+proj=cea +datum=WGS84")

plot(r_combined$total_richness)
plot(r_combined$canopy_height)

df_combined = df_combined %>% 
  mutate(color = case_when(biome2 == "Tropical & subtropical moist" ~ "darkgreen",
                           biome2 == "Tropical & subtropical dry" ~ "greenyellow",
                           biome2 == "Tropical & subtropical conifer" ~ "olivedrab4",
                           biome2 == "Tropical & subtropical savannas & shrublands" ~ "darkorange3",
                           biome2 == "Temperate savanna, shrub, & scrub" ~ "orange1",
                           biome2 == "Temperate & arctic" ~ "deepskyblue3"),
         biome2 = factor(biome2, levels = c("Temperate & arctic",
                                            "Temperate savanna, shrub, & scrub",
                                            "Tropical & subtropical savannas & shrublands",
                                            "Tropical & subtropical dry",
                                            "Tropical & subtropical conifer",
                                            "Tropical & subtropical moist")),
         color = factor(color, levels = c("deepskyblue3", "orange1", "darkorange3",
                                          "greenyellow", "olivedrab4", "darkgreen")))


# GLOBAL RICH ~ CANOPY HEIGHT RELATIONSHIP --------------------------------

ggplot(df_combined, aes(canopy_height, log(total_richness))) +
  geom_point(pch = ".") +
  geom_smooth(method = "lm")

ggplot(df_combined, aes(log(canopy_height), log(total_richness))) +
  geom_point(pch = ".") +
  geom_smooth(method = "lm")


m1 = glm(total_richness ~ log_ch,
        data = df_combined,
        family = poisson())

m2 = MASS::glm.nb(total_richness ~ log_ch, data = df_combined)

m3 = glm(total_richness ~ canopy_height,
         data = df_combined,
         family = poisson())

m4 = MASS::glm.nb(total_richness ~ canopy_height, data = df_combined)

AIC(m1, m2, m3, m4) # model with log_canopy height due to power law relationship for SARs

nd = data.frame(canopy_height = seq(min(df_combined$canopy_height), max(df_combined$canopy_height), 0.1))
nd$pred = predict(m4, newdata = nd, type = "response")
ggplot() +
  geom_point(data = df_combined, aes(canopy_height, total_richness), pch = 20, alpha = 0.2) +
  geom_line(data = nd, aes(canopy_height, pred))

nd$log_ch = log(nd$canopy_height)
nd$pred = predict(m2, newdata = nd, type = "response")
ggplot() +
  geom_point(data = df_combined, aes(canopy_height, total_richness), pch = 20, alpha = 0.2) +
  geom_line(data = nd, aes(canopy_height, pred))

resid.deviance = residuals(m2, type = "deviance")
df_combined$resid = resid.deviance
r_combined = df_combined %>% 
  rast(crs = "+proj=cea +datum=WGS84")
plot(r_combined$resid)

ggplot() +
  geom_point(data = df_combined, aes(canopy_height, total_richness, color = color), pch = 20, alpha = 0.2) +
  geom_line(data = nd, aes(canopy_height, pred), linewidth = 1.5) +
  scale_color_identity(guide = "legend",
                       labels = levels(df_combined$biome2),
                       breaks = levels(df_combined$color)) +
  scale_y_continuous("Species Richness") +
  scale_x_continuous("Canopy Height (m)") +
  guides(color = guide_legend(nrow = 3, ncol = 2)) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "none")

resid_combined = ggplot() +
  geom_spatvector(data = wd.cea) +
  geom_spatraster(data = r_combined, aes(fill = resid)) +
  geom_spatvector(data = wd.cea, fill = NA, color = "black") +
  scale_fill_continuous_divergingx("RdYlBu", na.value = NA, rev = T) +
  scale_x_continuous(expand = c(0,0)) +
  labs(title = "Residuals") +
  theme_bw() + thm +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0,0),
        legend.justification = c(0, 0),
        legend.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.key.width = unit(3, "mm"),
        panel.background = element_blank(),
        plot.background = element_blank())

richness_combined = ggplot() +
  geom_spatvector(data = wd.cea) +
  geom_spatraster(data = r_combined, aes(fill = total_richness)) +
  geom_spatvector(data = wd.cea, fill = NA, color = "black") +
  scale_fill_continuous_sequential("Sunset", na.value = NA) +
  scale_x_continuous(expand = c(0,0)) +
  labs(title = "Vertebrate Richness") +
  theme_bw() + thm +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0,0),
        legend.justification = c(0, 0),
        legend.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.key.width = unit(3, "mm"),
        panel.background = element_blank(),
        plot.background = element_blank())

canopy_height = ggplot() +
  geom_spatvector(data = wd.cea) +
  geom_spatraster(data = r_combined, aes(fill = canopy_height)) +
  geom_spatvector(data = wd.cea, fill = NA, color = "black") +
  scale_fill_continuous_sequential("BluYl", na.value = NA) +
  scale_x_continuous(expand = c(0,0)) +
  labs(title = "Canopy Height") +
  theme_bw() + thm +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0,0),
        legend.justification = c(0, 0),
        legend.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.key.width = unit(3, "mm"),
        panel.background = element_blank(),
        plot.background = element_blank())


maps_combined = richness_combined + resid_combined + plot_layout(ncol = 2) & theme(panel.border = element_blank(), plot.margin = margin(0, 10, 0, 10))

richness_combined + canopy_height + resid_combined +
  plot_layout(nrow = 3) & 
  theme(panel.border = element_blank(), plot.margin = margin(0, 10, 0, 10))


# RICHNESS BY CANOPY HEIGHT RESIDUALS -------------------------------------

richness.ch.resids = foreach(i = unique(df$taxa)) %do% {
  d = df %>% filter(taxa == i) %>% 
    drop_na(vert.mean.ses, rich, canopy_height) %>% 
    filter(rich > 5 & canopy_height > 0) %>% 
    mutate(log_ch = log(canopy_height))
  
  # model richness by canopy height
  m1 = glm(rich ~ log_ch,
          data = d,
          family = poisson())
  
  m2 = MASS::glm.nb(rich ~ log_ch, data = d)
  
  AIC(m1, m2)
  
  # Deviance residuals
  resid.deviance = residuals(m2, type = "deviance")
  d$resid.deviance = resid.deviance
  
  r = d %>% 
    #dplyr::select(x,y,resid.deviance, vert.mean.ses, p.arb, vert.mean, biome2) %>%
    # mutate(rowid = 1:nrow(d)) %>% 
    # rast(crs = "+proj=cea +datum=WGS84") %>% 
    # project("epsg:4326") %>% 
    # as.data.frame(xy = T) %>% 
    mutate(taxa = i)
  r
}

richness.ch.resids = bind_rows(richness.ch.resids) %>% 
  mutate(taxa = factor(taxa, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")),
         color = case_when(biome2 == "Tropical & subtropical moist" ~ "darkgreen",
                           biome2 == "Tropical & subtropical dry" ~ "greenyellow",
                           biome2 == "Tropical & subtropical conifer" ~ "olivedrab4",
                           biome2 == "Tropical & subtropical savannas & shrublands" ~ "darkorange3",
                           biome2 == "Temperate savanna, shrub, & scrub" ~ "orange1",
                           biome2 == "Temperate & arctic" ~ "deepskyblue3"),
         biome2 = factor(biome2, levels = c("Temperate & arctic",
                                            "Temperate savanna, shrub, & scrub",
                                            "Tropical & subtropical savannas & shrublands",
                                            "Tropical & subtropical dry",
                                            "Tropical & subtropical conifer",
                                            "Tropical & subtropical moist")),
         color = factor(color, levels = c("deepskyblue3", "orange1", "darkorange3",
                                          "greenyellow", "olivedrab4", "darkgreen")))


# reproject
richness.ch.resids.wgs.a = richness.ch.resids %>% 
  dplyr::select(x,y,resid.deviance, vert.mean.ses, taxa) %>% 
  filter(taxa == "Amphibians") %>% 
  rast(crs = "+proj=cea +datum=WGS84") %>%
  project("epsg:4326") %>%
  as.data.frame(xy = T) %>% 
  mutate(taxa = "Amphibians")
richness.ch.resids.wgs.b = richness.ch.resids %>% 
  dplyr::select(x,y,resid.deviance, vert.mean.ses, taxa) %>% 
  filter(taxa == "Birds") %>% 
  rast(crs = "+proj=cea +datum=WGS84") %>%
  project("epsg:4326") %>%
  as.data.frame(xy = T) %>% 
  mutate(taxa = "Birds")
richness.ch.resids.wgs.m = richness.ch.resids %>% 
  dplyr::select(x,y,resid.deviance, vert.mean.ses, taxa) %>% 
  filter(taxa == "Mammals") %>% 
  rast(crs = "+proj=cea +datum=WGS84") %>%
  project("epsg:4326") %>%
  as.data.frame(xy = T) %>% 
  mutate(taxa = "Mammals")
richness.ch.resids.wgs.r = richness.ch.resids %>% 
  dplyr::select(x,y,resid.deviance, vert.mean.ses, taxa) %>% 
  filter(taxa == "Reptiles") %>% 
  rast(crs = "+proj=cea +datum=WGS84") %>%
  project("epsg:4326") %>%
  as.data.frame(xy = T) %>% 
  mutate(taxa = "Reptiles")

richness.ch.resids.wgs = bind_rows(richness.ch.resids.wgs.a,
                                   richness.ch.resids.wgs.b,
                                   richness.ch.resids.wgs.m,
                                   richness.ch.resids.wgs.r) %>% 
  mutate(taxa = factor(taxa, levels = c("Birds", "Mammals", "Reptiles","Amphibians")))


p2 = ggplot() +
  geom_raster(data = richness.ch.resids.wgs, aes(x,y, fill = resid.deviance)) +
  geom_spatvector(data = wd,fill = NA, color = "black") +
  coord_sf(crs = 4326) +
  scale_fill_continuous_divergingx(palette = "RdBu", rev = T,
                                   guide = guide_colorbar(title = "Residuals")) +
  facet_wrap(vars(taxa), ncol = 1, strip.position = "left") +
  scale_x_continuous("Residuals", expand = c(0,0)) +
  thm +
  theme(
    legend.position = "bottom",
    legend.title.position = "top",
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.title.y = element_blank())
p2.legend = ggpubr::get_legend(p2)
p2.legend = as_ggplot(p2.legend)
p2 = p2 + theme(legend.position = "none")


# Biome level richness models ------------------------------------------------------

richness.ch.mods = data.frame()
richness.ch.mods.resid = data.frame()
for(i in unique(df$taxa)) {
  for(b in unique(df$biome2)) {
    d = df %>% filter(taxa == i & biome2 == b) %>% 
      drop_na(vert.mean.ses, rich, canopy_height) %>% 
      filter(rich > 5 & canopy_height > 0) %>% 
      mutate(log_ch = log(canopy_height))
    
    # model richness by canopy height
    m1 = glm(rich ~ log_ch,
             data = d,
             family = poisson())
    
    m2 = MASS::glm.nb(rich ~ log_ch, data = d)
    
    
    AIC(m1, m2)
    
    resid.deviance = residuals(m2, type = "deviance")
    d$resid.deviance = resid.deviance
    
    nd = data.frame(canopy_height = seq(min(d$canopy_height), max(d$canopy_height), 0.1))
    nd$log_ch = log(nd$canopy_height)
    p = predict(m2, newdata = nd, type = "response")
    nd$richness.pred = p
    nd$biome2 = b
    nd$taxa = i
    
    richness.ch.mods = rbind(richness.ch.mods, nd)
    richness.ch.mods.resid = rbind(richness.ch.mods.resid, d)
  }
}

richness.ch.mods = richness.ch.mods %>% 
  mutate(color = case_when(biome2 == "Tropical & subtropical moist" ~ "darkgreen",
                           biome2 == "Tropical & subtropical dry" ~ "greenyellow",
                           biome2 == "Tropical & subtropical conifer" ~ "olivedrab4",
                           biome2 == "Tropical & subtropical savannas & shrublands" ~ "darkorange3",
                           biome2 == "Temperate savanna, shrub, & scrub" ~ "orange1",
                           biome2 == "Temperate & arctic" ~ "deepskyblue3"),
         taxa = factor(taxa, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")),
         biome2 = factor(biome2, levels = c("Temperate & arctic",
                                            "Temperate savanna, shrub, & scrub",
                                            "Tropical & subtropical savannas & shrublands",
                                            "Tropical & subtropical dry",
                                            "Tropical & subtropical conifer",
                                            "Tropical & subtropical moist")),
         color = factor(color, levels = c("deepskyblue3", "orange1", "darkorange3",
                                          "greenyellow", "olivedrab4", "darkgreen")))


df = df %>% 
  mutate(color = case_when(biome2 == "Tropical & subtropical moist" ~ "darkgreen",
                           biome2 == "Tropical & subtropical dry" ~ "greenyellow",
                           biome2 == "Tropical & subtropical conifer" ~ "olivedrab4",
                           biome2 == "Tropical & subtropical savannas & shrublands" ~ "darkorange3",
                           biome2 == "Temperate savanna, shrub, & scrub" ~ "orange1",
                           biome2 == "Temperate & arctic" ~ "deepskyblue3"),
         taxa = factor(taxa, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")),
         biome2 = factor(biome2, levels = c("Temperate & arctic",
                                            "Temperate savanna, shrub, & scrub",
                                            "Tropical & subtropical savannas & shrublands",
                                            "Tropical & subtropical dry",
                                            "Tropical & subtropical conifer",
                                            "Tropical & subtropical moist")),
         color = factor(color, levels = c("deepskyblue3", "orange1", "darkorange3",
                                          "greenyellow", "olivedrab4", "darkgreen"))) %>% 
  drop_na(vert.mean.ses, rich, canopy_height) %>% 
  filter(rich > 5 & canopy_height > 0)

library(ggnewscale)
p1 = ggplot() +
  geom_point(data = df, aes(x = canopy_height, y = rich, color = vert.mean.ses), alpha = 1, pch = ".") +
  scale_color_continuous_divergingx(palette = "RdGy",
                                    guide = guide_colorbar(position = "right",
                                                           theme = theme(legend.key.width = unit(3, "mm"),
                                                                         legend.key.height = unit(20, "mm"))),
                                    rev = T)+
  new_scale_color() +
  geom_line(data = richness.ch.mods, aes(x = canopy_height, y = richness.pred, color = color),
            linewidth = 1) +
  facet_rep_wrap(~taxa, scales = "free_y", ncol = 4) +
  scale_color_identity(guide = "legend",
                       labels = levels(richness.ch.mods$biome2),
                       breaks = levels(richness.ch.mods$color)) +
  scale_x_continuous("Canopy height (m)") +
  scale_y_continuous("Species richness") +
  guides(color = guide_legend(nrow = 3, ncol = 2)) +
  thm +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        legend.position = "bottom")

# p1.legend = ggpubr::get_legend(p1)
# p1.legend = as_ggplot(p1.legend)
#p1 = p1 + theme(legend.position = "none")

p = maps_combined / p1 + plot_layout(heights = c(1,0.5))

ggsave("figures/main_figs/total_richness_maps.svg", width = 180, height = 150, units = "mm")

# PLOT SES VERTICALITY ----------------------------------------------------

# p3 = ggplot(richness.ch.resids.wgs) +
#   geom_spatvector(data = wd) +
#   geom_raster(aes(x,y,fill = vert.mean.ses)) +
#   geom_spatvector(data = wd, fill = NA, color = "black") +
#   scale_x_continuous("SES Verticality", expand = c(0,0)) +
#   scale_fill_continuous_divergingx("spectral", na.value = NA, rev = T, guide = guide_colorbar(title = "SES Verticality")) +
#   facet_wrap(vars(taxa), ncol = 1, strip.position = "left") +
#   coord_sf(crs = "epsg:4326") +
#   thm + theme(
#     legend.position = "bottom",
#     legend.title.position = "top",
#     strip.background = element_blank(),
#     strip.text = element_blank(),
#     axis.text = element_blank(),
#     axis.title.y = element_blank(),
#     axis.ticks = element_blank())
# p3.legend = ggpubr::get_legend(p3)
# p3.legend = as_ggplot(p3.legend)
# p3 = p3 + theme(legend.position = "none")
# 


# COMBINE PLOTS -----------------------------------------------------------

# p2.grob = ggplotGrob(p2)
# p1.grob = ggplotGrob(p1)
# p3.grob = ggplotGrob(p3)
# 
# max_height <- unit.pmax(p2.grob$heights, p1.grob$heights, p3.grob$heights)
# p2.grob$heights <- max_height
# p1.grob$heights <- max_height
# p3.grob$heights <- max_height
# p1.grob$widths = p1.grob$widths
# p2.grob$widths = 2*p2.grob$widths
# p3.grob$widths = 2*p3.grob$widths
# 
# # Use grid.arrange to combine them while preserving aspect ratio for the maps
# grid.newpage()
# pfinal = cbind(p1.grob, p2.grob, p3.grob,  size = "first")
# grid.draw(pfinal)
# 
# 
# svg("figures/main_figs/fig2_vert_maps2/fig.svg", width = 5, height = 6)
# p2 + p3
# dev.off()
# 
# svg("figures/main_figs/fig2_vert_maps2/fig-2.svg", width = 2, height = 6)
# p1
# dev.off()
# 
# 
# svg("figures/main_figs/fig2_vert_maps2/sesvert_legend.svg", width = 6, height = 2)
# p3.legend
# dev.off()
# 
# svg("figures/main_figs/fig2_vert_maps2/resid_legend.svg", width = 6, height = 2)
# p2.legend
# dev.off()
# 
# svg("figures/main_figs/fig2_vert_maps2/biome_legend.svg", width = 6, height = 6)
# p1.legend
# dev.off()
# 
# 
# 
