library(terra)
library(tidyterra)
library(tidyverse)

# Supplementary figure: Plot the percent relative change ((vert.future - vert.present)/(vertmax - vertmin)*100).
# biomes are ordered along the y-axis from the greatest average positive change at the top to the greatest
# average negative change at the bottom. Averages for ordering the y-axis are calculate for mean or ses verticality across all taxa.


# Load SES vert and Mean vert models ----------------------------------------------------
load("results/sdmTMB_models/amphibians_sesvert.RData")
amph.sesvert = pred.f
amph.sesvert$vertvar = "SES Mean Verticality"
amph.sesvert$est.reldif = pred.f$est.dif/(max(pred$vert.mean.ses) - min(pred$vert.mean.ses))*100
amph.sesvert$biome = pred$biome

load("results/sdmTMB_models/amphibians_meanvert.RData")
amph.meanvert = pred.f
amph.meanvert$vertvar = "Mean Verticality"
amph.meanvert$est.reldif = pred.f$est.dif/(max(pred$vert.mean) - min(pred$vert.mean))*100
amph.meanvert$biome = pred$biome

amph.pdat = amph.sesvert %>% 
  bind_rows(amph.meanvert) %>%
  dplyr::select(x,y,est.dif, est.reldif, vertvar, biome) %>% 
  mutate(class = "Amphibians")


load("results/sdmTMB_models/mammals_sesvert.RData")
mammals.sesvert = pred.f
mammals.sesvert$vertvar = "SES Mean Verticality"
mammals.sesvert$est.reldif = pred.f$est.dif/(max(pred$vert.mean.ses) - min(pred$vert.mean.ses))*100
mammals.sesvert$biome = pred$biome


load("results/sdmTMB_models/mammals_meanvert.RData")
mammals.meanvert = pred.f
mammals.meanvert$vertvar = "Mean Verticality"
mammals.meanvert$est.reldif = pred.f$est.dif/(max(pred$vert.mean) - min(pred$vert.mean))*100
mammals.meanvert$biome = pred$biome

mammals.pdat = mammals.sesvert %>% 
  bind_rows(mammals.meanvert) %>%
  dplyr::select(x,y,est.dif, est.reldif, vertvar, biome) %>% 
  mutate(class = "Mammals")


load("results/sdmTMB_models/birds_sesvert.RData")
birds.sesvert = pred.f
birds.sesvert$vertvar = "SES Mean Verticality"
birds.sesvert$est.reldif = pred.f$est.dif/(max(pred$vert.mean.ses) - min(pred$vert.mean.ses))*100
birds.sesvert$biome = pred$biome

load("results/sdmTMB_models/birds_meanvert.RData")
birds.meanvert = pred.f
birds.meanvert$vertvar = "Mean Verticality"
birds.meanvert$est.reldif = pred.f$est.dif/(max(pred$vert.mean) - min(pred$vert.mean))*100
birds.meanvert$biome = pred$biome

birds.pdat = birds.sesvert %>% 
  bind_rows(birds.meanvert) %>%
  dplyr::select(x,y,est.dif, est.reldif, vertvar, biome) %>% 
  mutate(class = "Birds")


load("results/sdmTMB_models/reptiles_sesvert.RData")
reptiles.sesvert = pred.f
reptiles.sesvert$vertvar = "SES Mean Verticality"
reptiles.sesvert$est.reldif = pred.f$est.dif/(max(pred$vert.mean.ses) - min(pred$vert.mean.ses))*100
reptiles.sesvert$biome = pred$biome

load("results/sdmTMB_models/reptiles_meanvert.RData")
reptiles.meanvert = pred.f
reptiles.meanvert$vertvar = "Mean Verticality"
reptiles.meanvert$est.reldif = pred.f$est.dif/(max(pred$vert.mean) - min(pred$vert.mean))*100
reptiles.meanvert$biome = pred$biome

repts.pdat = reptiles.sesvert %>% 
  bind_rows(reptiles.meanvert) %>%
  dplyr::select(x,y,est.dif, est.reldif, vertvar, biome) %>% 
  mutate(class = "Reptiles")

# biome = read.csv("data/derivative_data/env_data.csv") %>% 
#   dplyr::select(x,y,biome) %>% 
#   mutate(biome = factor(biome)) %>% 
#   drop_na()
# 
# biomekey = biome %>% 
#   distinct(biome) %>% 
#   drop_na() %>% 
#   mutate(id = 1:14)
# 
# biome = left_join(biome, biomekey, by = "biome") 
# 
# biome = rast(biome %>% dplyr::select(x,y,id))

dat = bind_rows(birds.pdat, mammals.pdat, repts.pdat, amph.pdat) %>% 
  mutate(x = x*1e5, y = y*1e5)

#dat$id = terra::extract(biome, dat[, c("x", "y")])[,2]
dat = dat %>% 
  #left_join(biomekey, by = "id") %>% 
  dplyr::select(x,y, est.dif, est.reldif, biome, class, vertvar) %>% 
  group_by(biome, vertvar) %>% 
  mutate(biome.meandif = mean(est.dif, na.rm = T),
         biome.meanReldif = mean(est.reldif, na.rm = T))

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
# biomes are ordered from most positive sesvert.dif to most negative relative difference (where order is by average est.reldif for all taxa in all realms within biomes)
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
  geom_boxplot(aes(x = est.reldif, y = reorder(biome, biome.meanReldif)), fill = "forestgreen", alpha = 0.5, width = 0.25) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) + 
  scale_y_discrete("Biome") +
  scale_x_continuous("% relative change") +
  facet_grid(rows = vars(vertvar), cols = vars(class)) +
  #scale_fill_manual(values = c("tan4", "forestgreen")) +
  #facet_grid2(vars(vertvar), vars(class), scales = "free_x", independent = "x") +
  #scale_x_facet(ROW == 1, limits = c(-5.2, 5.2)) +
  #scale_x_facet(ROW == 2, limits = c(-5.2, 5.2)) +
  theme(panel.grid.major = element_line(color = "gray80"),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color = "black", fill = NA),
        strip.background = element_rect(color = "black"),
        legend.position = "none")

png("figures/supp_figs/dif_biome_forestOnly_forestSES.png", width = 4000, height = 2000, res = 300)
dif.biome
dev.off()





  
  
  
  
  
  
  
  
  








