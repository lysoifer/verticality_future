library(terra)
library(tidyverse)
library(data.table)
library(colorspace)
library(foreach)

env = read.csv("data/derivative_data/env_data.csv")

amph = fread("data/derivative_data/gridcell_data/env_forest/50_km/amph_comdat.csv") %>% 
  dplyr::select(vert.mean, vert.mean.ses, biome) %>% 
  mutate(class = "Amphibians")

bird = fread("data/derivative_data/gridcell_data/env_forest/50_km/birds_comdat.csv") %>% 
  dplyr::select(vert.mean, vert.mean.ses, biome) %>% 
  mutate(class = "Birds")

mammals = fread("data/derivative_data/gridcell_data/env_forest/50_km/mammals_comdat.csv") %>% 
  dplyr::select(vert.mean, vert.mean.ses, biome) %>% 
  mutate(class = "Mammals")

rept = fread("data/derivative_data/gridcell_data/env_forest/50_km/reptiles_comdat.csv") %>% 
  dplyr::select(vert.mean, vert.mean.ses, biome) %>% 
  mutate(class = "Reptiles")

df = bind_rows(amph, bird, mammals, rept) %>% 
  mutate(class = factor(class, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")))

df.mean = df %>% 
  group_by(biome, class) %>% 
  summarise(vert.mean.ses.mean = mean(vert.mean.ses, na.rm = T),
            vert.mean.mean = mean(vert.mean, na.rm = T)) %>% 
  mutate(class = factor(class, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")))

plt_vertses_biome = ggplot(df, aes(x = vert.mean.ses)) +
  geom_histogram() +
  geom_vline(data = df.mean, aes(xintercept = vert.mean.ses.mean), linetype = "dashed", color = "red3", linewidth = 1) +
  facet_grid(biome~class, scales = "free", labeller = label_wrap_gen(width = 5)) +
  scale_x_continuous("SES Mean Verticality") +
  theme(panel.grid.major = element_line(color = "gray80"),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color = "black", fill = NA),
        strip.background = element_rect(color = "black"))

png("figures/supp_figs/hist_vertses_biome.png", height = 300, width = 300, res = 300, unit = "mm")
plt_vertses_biome
dev.off()

ggplot(df, aes(x = vert.mean.ses)) +
  geom_histogram() +
  geom_vline(data = df.mean, aes(xintercept = vert.mean.ses.mean), linetype = "dashed", color = "red3", linewidth = 1) +
  facet_grid(biome~class, scales = "free", labeller = label_wrap_gen(width = 5)) +
  scale_x_continuous("SES Mean Verticality") +
  theme(panel.grid.major = element_line(color = "gray80"),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color = "black", fill = NA),
        strip.background = element_rect(color = "black"))


# histograms stacked ------------------------------------------------------


biomes = c(
  "Tropical & Subtropical Moist Broadleaf Forests",
  "Tropical & Subtropical Dry Broadleaf Forests",
  "Tropical & Subtropical Coniferous Forests",
  "Tropical & Subtropical Grasslands, Savannas & Shrublands",
  "Temperate Broadleaf & Mixed Forests",
  "Mediterranean Forests, Woodlands & Scrub",
  "Temperate Conifer Forests",
  "Temperate Grasslands, Savannas & Shrublands",
  "Boreal Forests/Taiga",
  "Montane Grasslands & Shrublands",
  "Flooded Grasslands & Savannas",
  "Deserts & Xeric Shrublands",
  "Tundra"
)

library(RColorBrewer)
library(cowplot)
dummy = data.frame(vert.mean = c(0,1), vert.mean.ses = range(df$vert.mean.ses, na.rm = T)) %>% 
  rename("Mean Verticality" = vert.mean, "SES Mean Verticality" = vert.mean.ses) %>% 
  pivot_longer(cols = 1:2, names_to = "vert", values_to = "vals")

p = foreach(i = unique(biomes)) %do% {
  df %>% 
    rename("Mean Verticality" = vert.mean, "SES Mean Verticality" = vert.mean.ses) %>% 
    pivot_longer(cols = 1:2, names_to = "vert", values_to = "vals") %>% 
    mutate(class = factor(class, levels = c("Birds", "Mammals", "Reptiles", "Amphibians"))) %>% 
    filter(biome == i) %>% 
    ggplot() +
    geom_histogram(aes(x = vals, fill = class)) +
    geom_blank(data = dummy, aes(x = vals)) +
    #geom_vline(data = sesvert.summ %>% filter(biome == i), aes(xintercept = sesvert.mean, color = class)) +
    scale_fill_manual(values = brewer.pal(4, "Set2")) +
    scale_color_manual(values = brewer.pal(4, "Set2")) +
    scale_x_continuous("") +
    scale_y_continuous(expand = c(0,0), n.breaks = 3) +
    facet_grid(rows = vars(class), cols = vars(vert), scales = "free") +
    ggtitle(i) +
    theme_classic() +
    theme(panel.background = element_rect(color = "black", fill = NA, linewidth = 1),
          panel.grid.major = element_line(color = "gray80"),
          plot.title = element_text(hjust = 0.5),
          strip.text.y = element_blank(),
          panel.spacing = unit(0, units = "mm"),
          strip.background.x = element_rect(fill = "grey80", color = "black"),
          legend.key.size = unit(10, units = "mm"),
          legend.title = element_blank(),
          legend.text = element_text(size = 12))
}

legend = get_legend(p[[1]])

p = foreach(i = 1:length(p)) %do% {
  p[[i]] +
    theme(legend.position = "none")
}

library(patchwork)
p = wrap_plots(p) + legend + plot_annotation(tag_levels = "A")

png("figures/supp_figs/hist_biome.png", width = 550, height = 400, res = 300, unit = "mm")
p
dev.off()



# average vert per biome --------------------------------------------------
dat = bind_rows(bird, mammals, rept, amph) %>% 
  dplyr::select(biome, class, vert.mean, vert.mean.ses) %>%
  pivot_longer(cols = 3:4, names_to = "vertvar", values_to = "vals") %>% 
  group_by(biome, vertvar) %>% 
  mutate(biome.meanvert = mean(vals, na.rm = T))

dummy = data.frame(vert.mean = c(0,1), vert.mean.ses = c(-7.5, 6)) %>% 
  pivot_longer(cols = 1:2, names_to = "vertvar", values_to = "vert") %>% 
  mutate(vertvar = ifelse(vertvar == "vert.mean", "Mean Verticality", "SES Mean Verticality"))

library(ggh4x)
# biomes are ordered from most positive sesvert.dif to most negative sesvert.dif (where order is by average sesvert.dif for all taxa in all realms within biomes)
vert.biome = dat %>% 
  mutate(class = factor(class, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")),
         vertvar = ifelse(vertvar == "vert.mean", "Mean Verticality", "SES Mean Verticality")) %>% 
  group_by(biome, class, biome.meanvert, vertvar) %>% 
  summarise(vert = mean(vals, na.rm = T),
            lci = mean(vals, na.rm = T) - sd(vals, na.rm = T),
            uci = mean(vals, na.rm = T) + sd(vals, na.rm = T)) %>% 
  filter(vertvar == "SES Mean Verticality") %>% 
  ggplot() +
  geom_pointrange(aes(y = reorder(biome, biome.meanvert), x = vert, xmin = lci, xmax = uci), size = 0.15) +
  geom_blank(data = dummy, aes(x = vert)) +
  #geom_vline(data = data.frame(x = c(0,0.5), vertvar = c("SES Mean Verticality", "Mean Verticality")), aes(xintercept = x), linetype = "dashed", color = "red4", linewidth = 1) + 
  geom_vline(data = data.frame(x = c(0), vertvar = c("SES Mean Verticality")), aes(xintercept = x), linetype = "dashed", color = "red4", linewidth = 1) + 
  scale_y_discrete("Biome") +
  scale_x_continuous("SES Verticality", limits = c(-6,6), breaks = seq(-6,6,2)) +
  #scale_color_manual(values = c("tan4", "forestgreen")) +
  facet_wrap(~class, scales = "free_x", nrow = 1) +
  #facet_grid2(vars(vertvar), vars(class), scales = "free_x", independent = "x") +
  theme(panel.grid = element_line(color = "gray"),
        panel.background = element_rect(color = "black", fill = NA),
        strip.background = element_rect(color = "black"),
        legend.position = "none",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 8))

png("figures/supp_figs/vert_biome_forestOnly_forestSES.png", width = 180, height = 70, res = 300, unit = "mm")
vert.biome
dev.off()

# histograms present vs future --------------------------------------------

load("results/sdmTMB_models/amphibians_sesvert.RData")
amph.sesvert = pred
amph.sesvert$est.f = pred.f$est
amph.sesvert$class = "Amphibians"
amph.sesvert$biome = pred$biome

load("results/sdmTMB_models/mammals_sesvert.RData")
mammals.sesvert = pred
mammals.sesvert$est.f = pred.f$est
mammals.sesvert$class = "Mammals"
mammals.sesvert$biome = pred$biome

load("results/sdmTMB_models/birds_sesvert.RData")
birds.sesvert = pred
birds.sesvert$est.f = pred.f$est
birds.sesvert$class = "Birds"
birds.sesvert$biome = pred$biome


load("results/sdmTMB_models/reptiles_sesvert.RData")
repts.sesvert = pred
repts.sesvert$est.f = pred.f$est
repts.sesvert$class = "Reptiles"
repts.sesvert$biome = pred$biome


sesvert.df = bind_rows(amph.sesvert, mammals.sesvert, birds.sesvert, repts.sesvert)

sesvert.summ = sesvert.df %>% 
  select(est, est.f, biome, class) %>% 
  pivot_longer(cols = 1:2, names_to = "time", values_to = "vert") %>% 
  mutate(time = ifelse("est", "Predicted Present", "Predicted Future")) %>% 
  #mutate(time = factor(time, levels = c("Predicted Present", "Predicted Future"))) %>% 
  group_by(time, class, biome) %>% 
  summarise(sesvert.mean = mean(vert))

biomes = c(
  "Tropical & Subtropical Moist Broadleaf Forests",
  "Tropical & Subtropical Dry Broadleaf Forests",
  "Tropical & Subtropical Coniferous Forests",
  "Tropical & Subtropical Grasslands, Savannas & Shrublands",
  "Temperate Broadleaf & Mixed Forests",
  "Mediterranean Forests, Woodlands & Scrub",
  "Temperate Conifer Forests",
  "Temperate Grasslands, Savannas & Shrublands",
  "Boreal Forests/Taiga",
  "Montane Grasslands & Shrublands",
  "Flooded Grasslands & Savannas",
  "Deserts & Xeric Shrublands",
  "Tundra"
)

library(RColorBrewer)
p = foreach(i = unique(biomes)) %do% {
  sesvert.df %>% 
  select(est, est.f, biome, class) %>% 
  pivot_longer(cols = 1:2, names_to = "time", values_to = "vert") %>% 
  mutate(time = ifelse(time == "est", "Predicted Present", "Predicted Future")) %>% 
  mutate(time = factor(time, levels = c("Predicted Present", "Predicted Future")),
         class = factor(class, levels = c("Birds", "Mammals", "Reptiles", "Amphibians"))) %>% 
  filter(biome == i) %>% 
  ggplot(aes(x = vert, fill = class)) +
  geom_histogram() +
  #geom_vline(data = sesvert.summ %>% filter(biome == i), aes(xintercept = sesvert.mean, color = class)) +
  scale_fill_manual(values = brewer.pal(4, "Set2")) +
  scale_color_manual(values = brewer.pal(4, "Set2")) +
  scale_x_continuous("SES Mean Verticality") +
  scale_y_continuous(expand = c(0,0), n.breaks = 3) +
  facet_grid(rows = vars(class), cols = vars(time), scales = "free_y") +
  ggtitle(i) +
  theme_classic() +
  theme(panel.background = element_rect(color = "black", fill = NA),
        panel.grid.major = element_line(color = "gray80"),
        plot.title = element_text(hjust = 0.5),
        strip.text.y = element_blank(),
        panel.spacing = unit(0, units = "mm"),
        strip.background.x = element_rect(fill = "grey80", color = "black"))
}

legend = get_legend(p[[1]])

p = foreach(i = 1:length(p)) %do% {
  p[[i]] +
    theme(legend.position = "none")
}

library(patchwork)
p = wrap_plots(p) + legend + plot_annotation(tag_levels = "A")

png("figures/supp_figs/vert_difs/sesvert_difs_histograms.png", width = 500, height = 400, res = 300, unit = "mm")
p
dev.off()





















