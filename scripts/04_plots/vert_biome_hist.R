library(terra)
library(tidyverse)
library(data.table)
library(colorspace)
library(foreach)

env = read.csv("data/derivative_data/env_data.csv")

amph = fread("data/derivative_data/gridcell_data/amphibians_comdat/amph_comdat_parallel_forestsOnly.csv") %>% 
  dplyr::select(vert.mean, vert.mean.ses, biome) %>% 
  mutate(class = "Amphibians")

bird = fread("data/derivative_data/gridcell_data/birds_comdat/birds_comdat_parallel_elton_forestsOnly.csv") %>% 
  dplyr::select(vert.mean, vert.mean.ses, biome) %>% 
  mutate(class = "Birds")

mammals = fread("data/derivative_data/gridcell_data/mammals_comdat/mammals_comdat_parallel_forestsOnly.csv") %>% 
  dplyr::select(vert.mean, vert.mean.ses, biome) %>% 
  mutate(class = "Mammals")

rept = fread("data/derivative_data/gridcell_data/reptiles_comdat/rept_comdat_parallel_forestsOnly.csv") %>% 
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
dummy = data.frame(vert.mean = c(0,1), vert.mean.ses = range(df$vert.mean.ses)) %>% 
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
    scale_y_continuous(expand = c(0,0)) +
    facet_wrap(~vert, nrow = 1, scales = "free_x") +
    ggtitle(i) +
    theme_classic() +
    theme(panel.background = element_rect(color = "black", fill = NA),
          panel.grid.major = element_line(color = "gray80"),
          plot.title = element_text(hjust = 0.5))
}

legend = get_legend(p[[1]])

p = foreach(i = 1:length(p)) %do% {
  p[[i]] +
    theme(legend.position = "none")
}

library(patchwork)
p = wrap_plots(p) + legend + plot_annotation(tag_levels = "A")

png("figures/supp_figs/hist_biome.png", width = 550, height = 300, res = 300, unit = "mm")
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
  ggplot() +
  geom_pointrange(aes(y = reorder(biome, biome.meanvert), x = vert, xmin = lci, xmax = uci), size = 0.25) +
  geom_blank(data = dummy, aes(x = vert)) +
  #geom_vline(xintercept = 0, linetype = "dashed", color = "red4", linewidth = 1) + 
  scale_y_discrete("Biome") +
  scale_x_continuous("") +
  #scale_color_manual(values = c("tan4", "forestgreen")) +
  facet_grid2(vars(vertvar), vars(class), scales = "free_x", independent = "x") +
  theme(panel.grid = element_line(color = "gray"),
        panel.background = element_rect(color = "black", fill = NA),
        strip.background = element_rect(color = "black"),
        legend.position = "none")

png("figures/supp_figs/vert_biome_forestOnly_forestSES.png", width = 300, height = 150, res = 300, unit = "mm")
vert.biome
dev.off()

# histograms present vs future --------------------------------------------

load("results/sar_mods_forestOnly_forestSES/amphibians/amphibians_sar_sesvert.RData")
amph.sesvert = pred.df
amph.sesvert$class = "Amphibians"

load("results/sar_mods_forestOnly_forestSES/mammals/mammals_sar_sesvert.RData")
mammals.sesvert = pred.df
mammals.sesvert$class = "Mammals"

load("results/sar_mods_forestOnly_forestSES/birds/birds_sar_sesvert.RData")
birds.sesvert = pred.df
birds.sesvert$class = "Birds"

load("results/sar_mods_forestOnly_forestSES/reptiles/reptiles_sar_sesvert.RData")
repts.sesvert = pred.df
repts.sesvert$class = "Reptiles"



sesvert.df = bind_rows(amph.sesvert, mammals.sesvert, birds.sesvert, repts.sesvert)
sesvert.df = sesvert.df %>% 
  left_join(env, by = c("x", "y")) %>% 
  dplyr::select(pred.pres.trend, pred.future, biome, class) %>% 
  rename("Predicted Present Trend" = pred.pres.trend, "Predicted Future Trend" = pred.future) %>% 
  pivot_longer(1:2, names_to = "time", values_to = "vert") %>% 
  mutate(vert.type = "SES Mean Verticality")

sesvert.df %>% 
  group_by(biome) %>% 
  summarise(n = n())

sesvert.df %>% 
  mutate(time = factor(time, levels = c("Predicted Present Trend", "Predicted Future Trend"))) %>% 
  #filter(class == "Birds") %>% 
  ggplot(aes(x = vert, fill = class, group = class)) +
  geom_histogram(alpha = 1, position = "stack") +
  facet_wrap(~time, nrow = 1) +
  scale_fill_discrete_qualitative(palette = "Set2") +
  theme_classic()
  
sesvert.df %>% 
  mutate(time = factor(time, levels = c("Predicted Present Trend", "Predicted Future Trend"))) %>% 
  ggplot(aes(x = vert, fill = class, group = class)) +
  geom_histogram() +
  facet_wrap(biome~time) +
  scale_fill_discrete_qualitative(palette = "Set2") +
  theme_classic()

sesvert.summ = sesvert.df %>% 
  mutate(time = factor(time, levels = c("Predicted Present Trend", "Predicted Future Trend"))) %>% 
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
  mutate(time = factor(time, levels = c("Predicted Present Trend", "Predicted Future Trend")),
         class = factor(class, levels = c("Birds", "Mammals", "Reptiles", "Amphibians"))) %>% 
  filter(biome == i) %>% 
  ggplot(aes(x = vert, fill = class)) +
  geom_histogram() +
  #geom_vline(data = sesvert.summ %>% filter(biome == i), aes(xintercept = sesvert.mean, color = class)) +
  scale_fill_manual(values = brewer.pal(4, "Set2")) +
  scale_color_manual(values = brewer.pal(4, "Set2")) +
  scale_x_continuous("SES Mean Verticality", limits = c(-5.5, 2)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~time, nrow = 1) +
  ggtitle(i) +
  theme_classic() +
  theme(panel.background = element_rect(color = "black", fill = NA),
        panel.grid.major = element_line(color = "gray80"),
        plot.title = element_text(hjust = 0.5))
}

legend = get_legend(p)

p = foreach(i = 1:length(p)) %do% {
  p[[i]] +
    theme(legend.position = "none")
}

library(patchwork)
p = wrap_plots(p) + legend + plot_annotation(tag_levels = "A")

png("figures/vert_difs/sesvert_difs_histograms.png", width = 500, height = 400, res = 300, unit = "mm")
p
dev.off()

load("results/amphibians_sar_vertmean2.RData")
amph.meanvert = vert.scale.sub
amph.meanvert$class = "Amphibians"

load("results/mammals_sar_vertmean2.RData")
mammals.meanvert = vert.scale.sub
mammals.meanvert$class = "Mammals"

load("results/birdselton_sar_vertmean.RData")
birds.meanvert = vert.scale.sub
birds.meanvert$class = "Birds"

load("results/reptiles_sar_vertmean2.RData")
repts.meanvert = vert.scale.sub
repts.meanvert$class = "Reptiles"

meanvert.df = bind_rows(amph.meanvert, mammals.meanvert, birds.meanvert, repts.meanvert)
meanvert.df = meanvert.df %>% 
  dplyr::select(vert.mean, vert.mean.future)



















