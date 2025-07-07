# GOAL: examine correlations of verticality with both canopy height and species richness
# to evaluate expansion and packing across vertical strata within biomes

library(tidyverse)
library(sdmTMB)
library(lme4)
library(terra)
library(car)
library(data.table)
library(DHARMa)
library(colorspace)
source("scripts/00_functions/manuscript_functions.R")


amph = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/amph_comdat.csv")
birds = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/birds_comdat.csv")
rept = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/reptiles_comdat.csv")
mammals = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/mammals_comdat.csv")

amph$taxa = "Amphibians"
birds$taxa = "Birds"
rept$taxa = "Reptiles"
mammals$taxa = "Mammals"

amph = amph %>% 
  # subset to richness >= 5
  filter(rich >= 5) %>% 
  filter(!is.na(vert.mean.ses)) %>% 
  drop_na() %>% 
  dplyr::select(x,y,rich, vert.mean, vert.mean.ses, ecoregion, biome, realm, canopy_height) %>% 
  mutate(x = x/1e5,
         y = y/1e5,
         ecoregion = factor(ecoregion),
         biome = factor(biome),
         biorealm = factor(paste(biome, realm, sep = "_")),
         realm = factor(realm))

birds = birds %>% 
  # subset to richness >= 5
  filter(rich >= 5) %>% 
  filter(!is.na(vert.mean.ses)) %>% 
  drop_na() %>% 
  dplyr::select(x,y,rich, vert.mean, vert.mean.ses, ecoregion, biome, realm, canopy_height) %>% 
  mutate(x = x/1e5,
         y = y/1e5,
         ecoregion = factor(ecoregion),
         biome = factor(biome),
         biorealm = factor(paste(biome, realm, sep = "_")),
         realm = factor(realm))

rept = rept %>% 
  # subset to richness >= 5
  filter(rich >= 5) %>% 
  filter(!is.na(vert.mean.ses)) %>% 
  drop_na() %>% 
  dplyr::select(x,y,rich, vert.mean, vert.mean.ses, ecoregion, biome, realm, canopy_height) %>% 
  mutate(x = x/1e5,
         y = y/1e5,
         ecoregion = factor(ecoregion),
         biome = factor(biome),
         biorealm = factor(paste(biome, realm, sep = "_")),
         realm = factor(realm))

mammals = mammals %>% 
  # subset to richness >= 5
  filter(rich >= 5) %>% 
  filter(!is.na(vert.mean.ses)) %>% 
  drop_na() %>% 
  dplyr::select(x,y,rich, vert.mean, vert.mean.ses, ecoregion, biome, realm, canopy_height) %>% 
  mutate(x = x/1e5,
         y = y/1e5,
         ecoregion = factor(ecoregion),
         biome = factor(biome),
         biorealm = factor(paste(biome, realm, sep = "_")),
         realm = factor(realm))


ggplot(mammals, aes(x = canopy_height, y = vert.mean.ses, color = biome)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw()

amph %>% 
  group_by(biome) %>% 
  summarise(count = n())
rept %>% 
  group_by(biome) %>% 
  summarise(count = n())
birds %>% 
  group_by(biome) %>% 
  summarise(count = n())
mammals %>% 
  group_by(biome) %>% 
  summarise(count = n())


dat = list(amph, birds, rept, mammals)
taxa = c("amph", "birds", "rept", "mammals")
biomes = unique(birds$biome)
#cv = matrix(ncol = 4, nrow = 13)
#colnames(cv) = c("amph", "birds", "rept", "mammals")
#rownames(cv) = unique(birds$biome)
res = list()
count = 1
for(t in 1:4) {
  for(b in 1:length(biomes)) {
    if(biomes[b] %in% dat[[t]]$biome) {
      d = dat[[t]] %>% 
        filter(biome == as.character(biomes[b]))
      if(nrow(d)>5) {
        cor_cvses = cor.test(d$canopy_height, d$vert.mean.ses)
        cor_cvmean = cor.test(d$canopy_height, d$vert.mean)
        cor_rvses = cor.test(d$rich, d$vert.mean.ses)
        cor_rvmean = cor.test(d$rich, d$vert.mean)
        #cv[as.character(biomes[b]), a[i]] = cor$estimate
        res[[count]] = c(taxa[t], as.character(biomes[b]), 
                         cor_cvses$estimate, cor_cvses$p.value,
                         cor_cvmean$estimate, cor_cvmean$p.value,
                         cor_rvses$estimate, cor_rvses$p.value,
                         cor_rvmean$estimate, cor_rvmean$p.value)
        count = count + 1
      }
    }
  }
}

df = do.call("rbind", res)
colnames(df) = c("taxa", "biome", "corcv_ses", "pcv_ses",
                 "corcv_mean", "pcv_mean",
                 "corrv_ses", "prv_ses",
                 "corrv_mean", "prv_mean")
df = df %>% 
  as.data.frame() %>% 
  mutate_at(.vars = 3:10, .funs = function(x){as.numeric(x)}) %>% 
  pivot_longer(cols = 3:10, names_to = c(".value", "var"), names_sep = "_") %>% 
  mutate(region = case_when(grepl("Tropical", biome) ~ "Tropical",
                            grepl("Temperate", biome) ~ "Temperate",
                            grepl("Boreal", biome) ~ "Boreal",
                            biome == "Tundra" ~ "Boreal",
                            .default = "Other"),
         region = factor(region, levels = c("Boreal", "Temperate", "Tropical", "Other")),
         fillcol = case_when(biome == "Boreal Forests/Taiga" ~ "skyblue2",
                             biome == "Tundra" ~ "gray",
                             biome == "Temperate Conifer Forests" ~ "purple",
                             biome == "Temperate Broadleaf & Mixed Forests" ~ "orange3",
                             biome == "Temperate Grasslands, Savannas & Shrublands" ~ "goldenrod1",
                             biome == "Deserts & Xeric Shrublands" ~ "pink",
                             biome == "Flooded Grasslands & Savannas" ~ "lightblue1",
                             biome == "Mediterranean Forests, Woodlands & Scrub" ~ "red",
                             biome == "Montane Grasslands & Shrublands" ~ "tan",
                             biome == "Tropical & Subtropical Moist Broadleaf Forests" ~ "forestgreen",
                             biome == "Tropical & Subtropical Coniferous Forests" ~ "maroon",
                             biome == "Tropical & Subtropical Dry Broadleaf Forests"  ~ "yellowgreen",
                             biome == "Tropical & Subtropical Grasslands, Savannas & Shrublands" ~ "yellow"),
         biome = case_when(biome == "Tropical & Subtropical Grasslands, Savannas & Shrublands" ~ "Tropical & Subtropical Grasslands,\nSavannas & Shrublands",
                           biome == "Tropical & Subtropical Dry Broadleaf Forests" ~ "Tropical & Subtropical Dry\nBroadleaf Forests",
                           biome == "Tropical & Subtropical Coniferous Forests" ~ "Tropical & Subtropical\nConiferous Forests",
                           biome == "Tropical & Subtropical Moist Broadleaf Forests" ~ "Tropical & Subtropical\nMoist Broadleaf Forests",
                           biome == "Mediterranean Forests, Woodlands & Scrub" ~ "Mediterranean Forests,\nWoodlands & Scrub",
                           biome == "Temperate Grasslands, Savannas & Shrublands" ~ "Temperate Grasslands,\nSavannas & Shrublands",
                           biome == "Temperate Broadleaf & Mixed Forests" ~ "Temperate Broadleaf &\nMixed Forests",
                           .default = biome)) %>% 
  drop_na() %>% 
  filter(var == "ses")

# apply bonferroni correction
nbiome = df %>% 
  group_by(taxa, biome) %>% 
  distinct() %>% 
  group_by(taxa) %>% 
  summarize(nbiome = n())

df = left_join(df, nbiome, by = "taxa") %>% 
  mutate(sigval = 0.05/nbiome)

plt = df %>% 
  mutate(taxa = case_when(taxa == "birds" ~ "Birds",
                          taxa == "mammals" ~ "Mammals",
                          taxa == "rept" ~ "Reptiles",
                          taxa == "amph" ~ "Amphibians"),
         taxa = factor(taxa, levels = c("Birds", "Mammals", "Reptiles", "Amphibians"))) %>% 
  ggplot(aes(x = corrv, y = corcv)) +
  geom_point(aes(color = pcv<sigval), pch = 21, size = 5, stroke = 3) +
  geom_point(aes(color = prv<sigval, fill = fillcol), pch = 21, size = 4, stroke = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("white", "black")) +
  scale_fill_identity(guide = "legend",
                      labels = df$biome,
                      breaks = df$fillcol) +
  scale_x_continuous("r (richness ~ verticality)") +
  scale_y_continuous("r (verticality ~ canopy height)") +
  facet_wrap(vars(taxa)) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_blank())

svg("figures/main_figs/fig4/fig4_correlation_coefs.svg", width = 12, height = 6)
plt
dev.off()


# PLOT RICH~VERT CORRELATION COEF BY DRY SEASON PRECIP --------------------------------
amph = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/amph_comdat.csv")
birds = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/birds_comdat.csv")
rept = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/reptiles_comdat.csv")
mammals = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/mammals_comdat.csv")

amph$taxa = "Amphibians"
birds$taxa = "Birds"
rept$taxa = "Reptiles"
mammals$taxa = "Mammals"

clim = rbind(amph, birds, rept, mammals)
clim = clim %>% 
  group_by(taxa, biome) %>% 
  summarise(precipdry_mean = mean(precip_dry, na.rm = T),
            precipdry_max = max(precip_dry, na.rm = T),
            precipdry_min = min(precip_dry, na.rm = T),
            precipdry_range = precipdry_max - precipdry_min,
            tmin_mean = mean(tmin_cold, na.rm = T)) %>% 
  mutate(biome = case_when(biome == "Tropical & Subtropical Grasslands, Savannas & Shrublands" ~ "Tropical & Subtropical Grasslands,\nSavannas & Shrublands",
                                   biome == "Tropical & Subtropical Dry Broadleaf Forests" ~ "Tropical & Subtropical Dry\nBroadleaf Forests",
                                   biome == "Tropical & Subtropical Coniferous Forests" ~ "Tropical & Subtropical\nConiferous Forests",
                                   biome == "Tropical & Subtropical Moist Broadleaf Forests" ~ "Tropical & Subtropical\nMoist Broadleaf Forests",
                                   biome == "Mediterranean Forests, Woodlands & Scrub" ~ "Mediterranean Forests,\nWoodlands & Scrub",
                                   biome == "Temperate Grasslands, Savannas & Shrublands" ~ "Temperate Grasslands,\nSavannas & Shrublands",
                                   biome == "Temperate Broadleaf & Mixed Forests" ~ "Temperate Broadleaf &\nMixed Forests",
                                   .default = biome))

dat = df %>% 
  mutate(taxa = case_when(taxa == "birds" ~ "Birds",
                          taxa == "mammals" ~ "Mammals",
                          taxa == "rept" ~ "Reptiles",
                          taxa == "amph" ~ "Amphibians")) %>% 
  left_join(clim, by = c("taxa", "biome"))

dat.cors = list()
for(i in unique(dat$taxa)) {
  d = dat %>% filter(taxa == i)
  dat.cors[[i]] = cor.test(d$precipdry_mean, d$corcv)
}

ggplot(dat, aes(x = precipdry_mean, y = corrv, color = taxa, fill = taxa)) +
  geom_smooth(method = "lm", alpha = 0.2, linewidth = 0.5) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_continuous("Mean precipitation of the\ndriest month (kg m\u207B\u00B2 month\u207B\u00B9)") +
  scale_y_continuous("Correlation coefficient between\nverticality and richness") +
  scale_color_discrete_qualitative(palette = "Dark 3") +
  colorspace::scale_fill_discrete_qualitative (palette = "Dark 3") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.title = element_text(size = 6),
        axis.text = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(3, "mm"),
        legend.position = "inside",
        legend.position.inside = c(0.125,0.89),
        legend.background = element_blank())

ggsave("figures/main_figs/fig4-2.png", height = 90, width = 90, units = "mm", dpi = 300)

dat %>% 
  mutate(fillcol = case_when(biome == "Boreal Forests/Taiga" ~ "skyblue2",
                             biome == "Tundra" ~ "gray",
                             biome == "Temperate Conifer Forests" ~ "purple",
                             biome == "Temperate Broadleaf &\nMixed Forests" ~ "orange3",
                             biome == "Temperate Grasslands,\nSavannas & Shrublands" ~ "goldenrod1",
                             biome == "Deserts & Xeric Shrublands" ~ "pink",
                             biome == "Flooded Grasslands & Savannas" ~ "lightblue1",
                             biome == "Mediterranean Forests,\nWoodlands & Scrub" ~ "red",
                             biome == "Montane Grasslands & Shrublands" ~ "tan",
                             biome == "Tropical & Subtropical\nMoist Broadleaf Forests" ~ "forestgreen",
                             biome == "Tropical & Subtropical\nConiferous Forests" ~ "maroon",
                             biome == "Tropical & Subtropical Dry\nBroadleaf Forests"  ~ "yellowgreen",
                             biome == "Tropical & Subtropical Grasslands,\nSavannas & Shrublands" ~ "yellow")) %>% 
  ggplot(aes(x = precipdry_mean, y = corrv)) +
  geom_smooth(method = "lm", alpha = 0.2, linewidth = 0.5, color = "black") +
  geom_point(size = 1.25, aes(color = fillcol)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_continuous("Mean precipitation of the\ndriest month (kg m\u207B\u00B2 month\u207B\u00B9)") +
  scale_y_continuous("Correlation coefficient between\nverticality and richness") +
  scale_color_identity(guide = "legend",
                       labels = df$biome,
                       breaks = df$fillcol) +
  facet_wrap(~taxa) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(3, "mm"),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "gray90", color = NA))

ggsave("figures/main_figs/fig4-3.png", height = 160, width = 180, units = "mm", dpi = 300)


dat %>% 
  ggplot(aes(x = precipdry_mean, y = tmin_mean, color = corrv)) +
  geom_point() +
  facet_wrap(~taxa) +
  scale_color_continuous_divergingx(palette = "spectral", guide = guide_colorbar(title = "richness~verticality\ncorrelation coef")) +
  scale_x_continuous("Precipitation of the driest month") +
  scale_y_continuous("Minimum temperature of the coldest month") +
  theme(legend.title = element_text(angle = 270),
        legend.title.position = "right")

ggsave("figures/main_figs/fig4-4.png", height = 140, width = 160, units = "mm", dpi = 300)


# SUPPLEMENTARY PLOT CORRELATIONS -----------------------------------------

mammals$taxa = "Mammals"
birds$taxa = "Birds"
amph$taxa = "Amphibians"
rept$taxa = "Reptiles"

df = bind_rows(mammals, birds, amph, rept)
df = df %>% 
  mutate(fillcol = case_when(biome == "Boreal Forests/Taiga" ~ "skyblue2",
                              biome == "Tundra" ~ "gray",
                              biome == "Temperate Conifer Forests" ~ "purple",
                              biome == "Temperate Broadleaf & Mixed Forests" ~ "orange3",
                              biome == "Temperate Grasslands, Savannas & Shrublands" ~ "goldenrod1",
                              biome == "Deserts & Xeric Shrublands" ~ "pink",
                              biome == "Flooded Grasslands & Savannas" ~ "lightblue1",
                              biome == "Mediterranean Forests, Woodlands & Scrub" ~ "red",
                              biome == "Montane Grasslands & Shrublands" ~ "tan",
                              biome == "Tropical & Subtropical Moist Broadleaf Forests" ~ "forestgreen",
                              biome == "Tropical & Subtropical Coniferous Forests" ~ "maroon",
                              biome == "Tropical & Subtropical Dry Broadleaf Forests"  ~ "yellowgreen",
                              biome == "Tropical & Subtropical Grasslands, Savannas & Shrublands" ~ "yellow"),
            biome = case_when(biome == "Tropical & Subtropical Grasslands, Savannas & Shrublands" ~ "Tropical & Subtropical Grasslands,\nSavannas & Shrublands",
                              biome == "Tropical & Subtropical Dry Broadleaf Forests" ~ "Tropical & Subtropical Dry\nBroadleaf Forests",
                              biome == "Tropical & Subtropical Coniferous Forests" ~ "Tropical & Subtropical\nConiferous Forests",
                              biome == "Tropical & Subtropical Moist Broadleaf Forests" ~ "Tropical & Subtropical\nMoist Broadleaf Forests",
                              biome == "Mediterranean Forests, Woodlands & Scrub" ~ "Mediterranean Forests,\nWoodlands & Scrub",
                              biome == "Temperate Grasslands, Savannas & Shrublands" ~ "Temperate Grasslands,\nSavannas & Shrublands",
                              biome == "Temperate Broadleaf & Mixed Forests" ~ "Temperate Broadleaf &\nMixed Forests",
                              .default = biome)) %>% 
  dplyr::select(rich, vert.mean.ses, biome, canopy_height, taxa, fillcol) %>% 
  pivot_longer(cols = c(rich, canopy_height), names_to = "var", values_to = "val") %>% 
  mutate(var = ifelse(var == "rich", "Richness", "Canopy height"),
         taxa = factor(taxa, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")))

df %>% 
  ggplot(aes(val, vert.mean.ses, color = fillcol, fill = fillcol)) +
  geom_point(pch = ".", alpha = 0.2) +
  geom_smooth(method = "lm", size = 0.5) +
  facet_grid2(taxa~var, axes = "all", remove_labels = "y", scales = "free", independent = "x", switch = "x") +
  scale_color_identity(guide = "legend", labels = df$biome, breaks = df$fillcol) +
  scale_fill_identity(guide = "legend", labels = df$biome, breaks = df$fillcol) +
  scale_y_continuous("SES Verticality") +
  theme_classic() +
  theme(strip.placement.x = "outside",
        strip.background.x = element_blank(),
        strip.background.y = element_rect(fill = "gray", color = "black"),
        panel.background = element_rect(fill = NA, color = "black"),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        strip.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.text = element_text(size = 6),
        axis.text = element_text(size = 6))

ggsave("figures/supp_figs/vert_richness_correlations.png", width = 130, height = 140, units = "mm", dpi = 600)





