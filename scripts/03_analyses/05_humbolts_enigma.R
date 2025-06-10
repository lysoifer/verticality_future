library(tidyverse)
library(Data4Ecologists)
library(terra)
library(foreach)
library(colorspace)
library(tidyterra)
library(MASS)

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
                            .default = "Temperate savanna, shrub, & scrub")) %>% 
 # drop_na(biome2) %>% 
  filter(rich > 5)

ggplot(df, aes(canopy_height, rich)) +
  geom_point(pch = ".", alpha = 0.2) +
  geom_smooth(method = "lm") +
  facet_wrap(~taxa, scales = "free_y") +
  theme_classic()

ggplot(df, aes(canopy_height, rich, color = biome2)) +
  geom_point(pch = ".", alpha = 0.2) +
  geom_smooth(method = "lm") +
  facet_wrap(~taxa, scales = "free_y") +
  theme_classic()

ggplot(df, aes(canopy_height, vert.mean.ses, color = biome)) +
  geom_point(pch = ".", alpha = 0.2) +
  geom_smooth(formula = y~log(x)) +
  facet_wrap(~taxa, scales = "free_y") +
  theme_classic()

richness.ch.resids = foreach(i = unique(df$taxa)) %do% {
  d = df %>% filter(taxa == i) %>% 
      drop_na(vert.mean.ses, rich, canopy_height) %>% 
      filter(rich > 5 & canopy_height > 0) %>% 
      mutate(log_ch = log(canopy_height))
  
  # model richness by canopy height
  # m1 = glm(rich ~ log_ch,
  #         data = d,
  #         family = poisson())
  
  m2 = MASS::glm.nb(rich ~ log_ch, data = d)
  
  # AIC(m1, m2)
  
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



wd = vect("data/original/rnaturalearth_world.shp") %>% 
  project("epsg:4326") %>% 
  crop(ext(-180,180,-60,83.64))

richness.ch.resids.wgs.a = richness.ch.resids %>% 
  dplyr::select(x,y,resid.deviance, taxa) %>% 
  filter(taxa == "Amphibians") %>% 
  rast(crs = "+proj=cea +datum=WGS84") %>%
  project("epsg:4326") %>%
  as.data.frame(xy = T) %>% 
  mutate(taxa = "Amphibians")
richness.ch.resids.wgs.b = richness.ch.resids %>% 
  dplyr::select(x,y,resid.deviance, taxa) %>% 
  filter(taxa == "Birds") %>% 
  rast(crs = "+proj=cea +datum=WGS84") %>%
  project("epsg:4326") %>%
  as.data.frame(xy = T) %>% 
  mutate(taxa = "Birds")
richness.ch.resids.wgs.m = richness.ch.resids %>% 
  dplyr::select(x,y,resid.deviance, taxa) %>% 
  filter(taxa == "Mammals") %>% 
  rast(crs = "+proj=cea +datum=WGS84") %>%
  project("epsg:4326") %>%
  as.data.frame(xy = T) %>% 
  mutate(taxa = "Mammals")
richness.ch.resids.wgs.r = richness.ch.resids %>% 
  dplyr::select(x,y,resid.deviance, taxa) %>% 
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

thm = theme(legend.key.height = unit(3, units = "mm"),
            legend.key.width = unit(10, units = "mm"),
            plot.background = element_blank(),
            plot.margin = unit(c(0,0,0,0), units = "mm"),
            panel.spacing = unit(0,"mm"),
            panel.background = element_rect(fill = "transparent", color = NA),
            strip.background = element_rect(fill = "white", color = "black"),
            panel.grid = element_blank(),
            legend.title = element_text(hjust = 0.5, size = 10),
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.text = element_blank())

p2 = ggplot() +
  geom_raster(data = richness.ch.resids.wgs, aes(x,y, fill = resid.deviance)) +
  geom_spatvector(data = wd,fill = NA, color = "black") +
  coord_sf(crs = 4326) +
  scale_fill_continuous_divergingx(palette = "RdBu", rev = T,
                                   guide = guide_colorbar(title = "Residuals")) +
  facet_wrap(vars(taxa), ncol = 1, strip.position = "left") +
  theme_bw() +
  thm + theme(
    legend.position = "bottom",
    legend.title.position = "top",
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.text.x = element_blank())
p2.legend = ggpubr::get_legend(p2)
p2.legend = as_ggplot(p2.legend)
p2 = p2 + theme(legend.position = "none")

ggsave("figures/main_figs/richness~CH_residuals.png", width = 180, height = 120, units = "mm", dpi = 300)


ggplot(data = richness.ch.resids, aes(p.arb, resid.deviance)) +
  geom_point(aes(color = color), pch = ".", alpha = 0.2) +
  geom_smooth(method = "lm") +
  scale_color_identity(guide = "legend",
                       labels = levels(richness.ch.mods$biome2),
                       breaks = levels(richness.ch.mods$color)) +
  facet_wrap(~taxa) +
  theme_bw()

ggplot(data = richness.ch.resids, aes(vert.mean, resid.deviance)) +
  geom_point(aes(color = color), pch = ".", alpha = 0.2) +
  geom_smooth(method = "lm") +
  scale_color_identity(guide = "legend",
                       labels = levels(richness.ch.mods$biome2),
                       breaks = levels(richness.ch.mods$color)) +
  facet_wrap(~taxa) +
  theme_bw()

ggplot(data = richness.ch.resids, aes(vert.mean.ses, resid.deviance)) +
  geom_point(aes(color = color), pch = ".", alpha = 0.2) +
  geom_smooth(method = "lm", color = "black") +
  scale_color_identity(guide = "legend",
                       labels = levels(richness.ch.mods$biome2),
                       breaks = levels(richness.ch.mods$color)) +
  scale_y_continuous("Residuals") +
  scale_x_continuous("Verticality") +
  facet_wrap(~taxa) +
  theme_bw() +
  theme(legend.position = "none")

richness.ch.resids %>% 
  mutate(biome3 = ifelse(grepl("Tropical", biome2), "Tropical", "Temperate/Arctic")) %>% 
  ggplot(aes(log10(precip_dry +1), resid.deviance, color = color)) +
  geom_point(pch = ".", alpha = 0.2) +
  geom_smooth(method = "lm") +
  scale_color_identity(guide = "legend",
                       labels = levels(richness.ch.mods$biome2),
                       breaks = levels(richness.ch.mods$color)) +
  facet_wrap(~taxa) +
  theme_bw()

richness.ch.resids %>% 
  mutate(biome3 = ifelse(grepl("Tropical", biome2), "Tropical", "Temperate/Arctic")) %>% 
  ggplot(aes(tmin_cold, resid.deviance, color = biome3)) +
  geom_point(pch = ".", alpha = 0.2) +
  geom_smooth(method = "lm") +
  # scale_color_identity(guide = "legend",
  #                      labels = levels(richness.ch.mods$biome2),
  #                      breaks = levels(richness.ch.mods$color)) +
  facet_wrap(~taxa) +
  theme_bw()


ggplot(data = richness.ch.resids, aes(tmin_cold, resid.deviance)) +
  geom_point(aes(color = color), pch = ".", alpha = 0.2) +
  geom_smooth(method = "lm", color = "black") +
  scale_color_identity(guide = "legend",
                       labels = levels(richness.ch.mods$biome2),
                       breaks = levels(richness.ch.mods$color)) +
  scale_y_continuous("Residuals") +
  scale_x_continuous("Minimum temperature of coldest month") +
  scale_color_identity(guide = "legend",
                       labels = levels(df_combined$biome2),
                       breaks = levels(df_combined$color)) +
  guides(color = guide_legend(nrow = 3, ncol = 2)) +
  facet_rep_wrap(~taxa, nrow = 1) +
  theme_classic() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_line(color = "black"))
  

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

richness.ch.mods.resid = richness.ch.mods.resid %>% 
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

ggplot() +
  geom_point(data = df, aes(x = canopy_height, y = rich, color = color), size = 0.2, alpha = 0.05, pch = 20) +
  geom_line(data = richness.ch.mods, aes(x = canopy_height, y = richness.pred, color = color)) +
  facet_wrap(~taxa, scales = "free_y", ncol = 1) +
  scale_color_identity(guide = "legend",
                       labels = levels(richness.ch.mods$biome2),
                       breaks = levels(richness.ch.mods$color)) +
  scale_x_continuous("Canopy height (m)") +
  scale_y_continuous("Species richness") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank())

ggsave("figures/main_figs/richness~canopy_height_biome2.png", width = 180, height = 120, units = "mm", dpi = 300)

ggplot(richness.ch.mods.resid, aes(log10(precip_dry+1), resid.deviance, color = color)) +
  geom_point(size = 0.2, alpha = 0.05, pch = 20) +
  geom_smooth(method = "lm") +
  facet_wrap(~taxa, scales = "free_y", ncol = 1) +
  scale_color_identity(guide = "legend",
                       labels = levels(richness.ch.mods$biome2),
                       breaks = levels(richness.ch.mods$color)) +
  scale_x_continuous("log10(Dry season precipitation)") +
  scale_y_continuous("Biome residuals") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank())



# Biome level verticality models ------------------------------------------------------

vert.ch.mods = data.frame()
for(i in unique(df$taxa)) {
  for(b in unique(df$biome2)) {
    d = df %>% filter(taxa == i & biome2 == b) %>% 
      drop_na(vert.mean.ses, rich, canopy_height) %>% 
      filter(rich > 5 & canopy_height > 0) %>% 
      mutate(log_ch = log(canopy_height))
    
    # model richness by canopy height
    m1 = glm(vert.mean.ses ~ canopy_height,
             data = d,
             family = gaussian())
    
    nd = data.frame(canopy_height = seq(min(d$canopy_height), max(d$canopy_height), 0.1))
    nd$log_ch = log(nd$canopy_height)
    p = predict(m1, newdata = nd, type = "response")
    nd$sesvert.pred = p
    nd$biome2 = b
    nd$taxa = i
    
    vert.ch.mods = rbind(vert.ch.mods, nd)
  }
}

vert.ch.mods = vert.ch.mods %>% 
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


ggplot() +
  geom_point(data = df, aes(x = canopy_height, y = vert.mean.ses, color = color), size = 0.2, alpha = 0.05, pch = 20) +
  geom_line(data = vert.ch.mods, aes(x = canopy_height, y = sesvert.pred, color = color)) +
  facet_wrap(~taxa, scales = "free_y") +
  scale_color_identity(guide = "legend",
                       labels = levels(vert.ch.mods$biome2),
                       breaks = levels(vert.ch.mods$color)) +
  scale_x_continuous("Canopy height (m)") +
  scale_y_continuous("Species richness") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank())

ggsave("figures/main_figs/vert~canopy_height_biome2.png", width = 180, height = 120, units = "mm", dpi = 300)

mat = rast("data/original/env_data/chelsa/CHELSA_bio1_1981-2010_V.2.1.tif")
plot(mat)
