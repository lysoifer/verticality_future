library(tidyverse)
library(Data4Ecologists)
library(terra)
library(foreach)
library(colorspace)
library(tidyterra)
library(MASS)
source("scripts/04_plots/fig3_coef_plot.R")


# Biome level richness models ------------------------------------------------------

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


