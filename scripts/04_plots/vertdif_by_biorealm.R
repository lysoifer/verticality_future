library(terra)
library(tidyterra)
library(tidyverse)


# Load SES vert and Mean vert models ----------------------------------------------------
load("results/sdmTMB_models/amphibians_sesvert.RData")
amph.sesvert = pred.f
amph.sesvert$vertvar = "SES Mean Verticality"

load("results/sdmTMB_models/amphibians_meanvert.RData")
amph.meanvert = pred.f
amph.meanvert$vertvar = "Mean Verticality"

amph.pdat = amph.sesvert %>% 
  bind_rows(amph.meanvert) %>%
  dplyr::select(x,y,est.dif, vertvar) %>% 
  mutate(class = "Amphibians")


load("results/sdmTMB_models/mammals_sesvert.RData")
mammals.sesvert = pred.f
mammals.sesvert$vertvar = "SES Mean Verticality"

load("results/sdmTMB_models/mammals_meanvert.RData")
mammals.meanvert = pred.f
mammals.meanvert$vertvar = "Mean Verticality"

mammals.pdat = mammals.sesvert %>% 
  bind_rows(mammals.meanvert) %>%
  dplyr::select(x,y,est.dif, vertvar) %>% 
  mutate(class = "Mammals")


load("results/sdmTMB_models/birds_sesvert.RData")
birds.sesvert = pred.f
birds.sesvert$vertvar = "SES Mean Verticality"

load("results/sdmTMB_models/birds_meanvert.RData")
birds.meanvert = pred.f
birds.meanvert$vertvar = "Mean Verticality"

birds.pdat = birds.sesvert %>% 
  bind_rows(birds.meanvert) %>%
  dplyr::select(x,y,est.dif, vertvar) %>% 
  mutate(class = "Birds")


load("results/sdmTMB_models/reptiles_sesvert.RData")
reptiles.sesvert = pred.f
reptiles.sesvert$vertvar = "SES Mean Verticality"

load("results/sdmTMB_models/reptiles_meanvert.RData")
reptiles.meanvert = pred.f
reptiles.meanvert$vertvar = "Mean Verticality"

repts.pdat = reptiles.sesvert %>% 
  bind_rows(reptiles.meanvert) %>%
  dplyr::select(x,y,est.dif, vertvar) %>% 
  mutate(class = "Reptiles")

biome = read.csv("data/derivative_data/env_data.csv") %>% 
  dplyr::select(x,y,biome) %>% 
  mutate(biome = factor(biome)) %>% 
  drop_na()

headbiomekey = biome %>% 
  distinct(biome) %>% 
  drop_na() %>% 
  mutate(id = 1:14)

biome = left_join(biome, biomekey, by = "biome") 

biome = rast(biome %>% dplyr::select(x,y,id))

dat = bind_rows(birds.pdat, mammals.pdat, repts.pdat, amph.pdat) %>% 
  mutate(x = x*1e5, y = y*1e5)

dat$id = terra::extract(biome, dat[, c("x", "y")])[,2]
dat = dat %>% 
  left_join(biomekey, by = "id") %>% 
  dplyr::select(x,y, est.dif, biome, class, vertvar) %>% 
  group_by(biome, vertvar) %>% 
  mutate(biome.meandif = mean(est.dif, na.rm = T))

# dat.sesvert = full_join(birds.pdat, mammals.pdat, by = c("x", "y")) %>% 
#   full_join(amph.pdat, by = c("x", "y")) %>% 
#   full_join(repts.pdat, by = c("x", "y")) %>% 
#   left_join(env, by = c("x", "y")) %>% 
#   dplyr::select(x, y, sesvert.dif.birds, sesvert.dif.mammals, sesvert.dif.amph, sesvert.dif.repts, biome, realm) %>% 
#   mutate(biorealm = paste(biome, realm, sep = "_")) %>% 
#   pivot_longer(cols = 3:6, names_to = "class", values_to = "sesvert.dif") %>% 
#   group_by(biome) %>% 
#   mutate(biome.meandif = mean(sesvert.dif, na.rm = T))


library(ggh4x)
# biomes are ordered from most positive sesvert.dif to most negative sesvert.dif (where order is by average sesvert.dif for all taxa in all realms within biomes)
dif.biome = dat %>% 
  mutate(class = factor(class, levels = c("Birds", "Mammals", "Reptiles", "Amphibians"))) %>% 
  drop_na(biome) %>% 
  # group_by(biome, class, biome.meandif, vertvar) %>% 
  # summarise(dif.mean = mean(est.dif, na.rm = T),
  #           dif.lci = mean(est.dif, na.rm = T) - sd(est.dif, na.rm = T),
  #           dif.uci = mean(est.dif, na.rm = T) + sd(est.dif, na.rm = T)) %>% 
  ggplot() +
  # geom_pointrange(aes(y = reorder(biome, biome.meandif), x = dif.mean, xmin = dif.lci, xmax = dif.uci,
  #                     color = dif.mean > 0), size = 0.5, linewidth = 1) +
  geom_boxplot(aes(x = est.dif, y = biome)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 1) + 
  scale_y_discrete("Biome") +
  scale_x_continuous("Verticality difference") +
  scale_fill_manual(values = c("tan4", "forestgreen")) +
  facet_grid2(vars(vertvar), vars(class), scales = "free_x", independent = "x") +
  scale_x_facet(ROW == 1, limits = c(-0.0255, 0.0255)) +
  scale_x_facet(ROW == 2, limits = c(-0.9, 0.9)) +
  theme(panel.grid.major = element_line(color = "gray"),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color = "black", fill = NA),
        strip.background = element_rect(color = "black"),
        legend.position = "none")

png("figures/supp_figs/dif_biome_forestOnly_forestSES.png", width = 4000, height = 2500, res = 300)
dif.biome
dev.off()





  
  
  
  
  
  
  
  
  








