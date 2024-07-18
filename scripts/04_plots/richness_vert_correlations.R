# plot verticality x richness for regions with canopy_coef > 1 and canopy_coef < 1

library(ggplot2)
library(tidyverse)
library(patchwork)

amph = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/amph_comdat.csv")
birds = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/birds_comdat.csv")
rept = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/reptiles_comdat.csv")
mammals = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/mammals_comdat.csv")

amph$taxa = "Amphibians"
birds$taxa = "Birds"
rept$taxa = "Reptiles"
mammals$taxa = "Mammals"

dat = rbind(amph, birds, rept, mammals)

dat %>% 
  filter(rich >= 5) %>% 
  ggplot(aes(x = rich, y = vert.mean.ses, color = biome)) +
    geom_point(pch = ".") +
    geom_smooth(method = "lm") +
    facet_grid(rows = vars(biome), cols = vars(taxa), scales = "free_x") +
    theme_classic()

dat %>% 
  filter(rich >= 5) %>% 
  ggplot(aes(x = rich, y = vert.mean, color = biome)) +
  geom_point(pch = ".") +
  geom_smooth(method = "lm") +
  facet_wrap(~taxa, scales = "free_x") +
  theme_classic()


plot_cor = function(mod, pred, taxa) {
  
  coefs = tidy(mod)
  coef.var = as.numeric(coefs[which(coefs$term == "canopy_height"), "estimate"])
  
  pred.canopy.pos = pred %>% 
    filter(zeta_s_canopy_height + coef.var > 1)
  
  pred.canopy.n = pred %>% 
    filter(zeta_s_canopy_height + coef.var <= 1)
  
  cor.p = cor.test(pred.canopy.pos$vert.mean.ses, pred.canopy.pos$rich)
  cor.n = cor.test(pred.canopy.n$vert.mean.ses, pred.canopy.n$rich)

  sesvert = pred %>% 
    ggplot(aes(x = vert.mean.ses, y = rich, color = (zeta_s_canopy_height + coef.var) > 1, 
               fill = (zeta_s_canopy_height + coef.var) > 1,
               alpha = (zeta_s_canopy_height + coef.var) > 1)) +
    geom_point(size = 0.05) +
    geom_smooth(method = "lm", alpha = 0.3) +
    scale_fill_manual(values = c("gray80", "blue")) +
    scale_color_manual(values = c("gray80", "blue")) +
    scale_alpha_manual(values = c(0.2, 1)) +
    annotate(geom = "text", label = paste0("R\u00b2 = ", format(round(cor.n$estimate^2, 2), nsmall = 2)),  
             x = -Inf, y = Inf, hjust = -0.5, vjust = 3, color = "gray50") +
    annotate(geom = "text", label = paste0("R\u00b2 = ", format(round(cor.p$estimate^2, 2), nsmall = 2)),  
             x = -Inf, y = Inf, hjust = -0.5, vjust = 5, color = "blue") +
    scale_x_continuous("SES Mean Verticality") +
    scale_y_continuous(paste0(taxa, " Richness")) +
    theme_classic()
  
  cor.p = cor.test(pred.canopy.pos$vert.mean, pred.canopy.pos$rich)
  cor.n = cor.test(pred.canopy.n$vert.mean, pred.canopy.n$rich)
  
  meanvert = pred %>% 
    ggplot(aes(x = vert.mean.ses, y = rich, color = (zeta_s_canopy_height + coef.var) > 1, 
               fill = (zeta_s_canopy_height + coef.var) > 1,
               alpha = (zeta_s_canopy_height + coef.var) > 1)) +
    geom_point(size = 0.05) +
    geom_smooth(method = "lm", alpha = 0.3) +
    scale_fill_manual(values = c("gray80", "blue")) +
    scale_color_manual(values = c("gray80", "blue")) +
    scale_alpha_manual(values = c(0.2, 1)) +
    annotate(geom = "text", label = paste0("R\u00b2 = ", format(round(cor.n$estimate^2, 2), nsmall = 2)),  
             x = -Inf, y = Inf, hjust = -0.5, vjust = 3, color = "gray50") +
    annotate(geom = "text", label = paste0("R\u00b2 = ", format(round(cor.p$estimate^2, 2), nsmall = 2)),  
             x = -Inf, y = Inf, hjust = -0.5, vjust = 5, color = "blue") +
    scale_x_continuous("SES Mean Verticality") +
    scale_y_continuous(paste0(taxa, " Richness")) +
    theme_classic()

  return(list(sesvert, meanvert))
}

plot_cor_biome = function(mod, pred, taxa) {
  
  #coefs = tidy(mod)
  #coef.var = as.numeric(coefs[which(coefs$term == "canopy_height"), "estimate"])
  
  pred.filter = pred %>% 
    filter(biome == "Tropical & Subtropical Moist Broadleaf Forests")
  
  cor = cor.test(pred$vert.mean.ses, pred$rich)
  cor.filter = cor.test(pred.filter$vert.mean.ses, pred.filter$rich)
  
  # set up legend
  cols = c("all points" = "gray80", "Tropical & subtropical\nmoist broadleaf forest" = "forestgreen")
  
  sesvert = pred %>% 
    ggplot(aes(x = vert.mean.ses, y = rich)) +
    geom_point(size = 0.05, alpha = 0.4, color = "gray80") +
    geom_point(data = pred.filter, color = "forestgreen", size = 0.05) +
    geom_smooth(method = "lm", alpha = 0.3, color = "gray60", aes(fill = "all points")) +
    geom_smooth(data = pred.filter, method = "lm", color = "forestgreen", aes(fill = "Tropical & subtropical\nmoist broadleaf forest"), alpha = 0.3) +
    scale_fill_manual("", values = cols) +
    scale_alpha_manual(values = c(0.2, 1)) +
    annotate(geom = "text", label = paste0("R\u00b2 = ", format(round(cor$estimate^2, 2), nsmall = 2)),  
             x = -Inf, y = Inf, hjust = -0.5, vjust = 3, color = "gray50") +
    annotate(geom = "text", label = paste0("R\u00b2 = ", format(round(cor.filter$estimate^2, 2), nsmall = 2)),  
             x = -Inf, y = Inf, hjust = -0.5, vjust = 5, color = "forestgreen") +
    scale_x_continuous("SES Mean Verticality") +
    scale_y_continuous(paste0(taxa, " Richness")) +
    theme_classic() + theme(panel.grid = element_blank())
  
  cor = cor.test(pred$vert.mean, pred$rich)
  cor.filter = cor.test(pred.filter$vert.mean, pred.filter$rich)
  
  meanvert = pred %>% 
    ggplot(aes(x = vert.mean, y = rich)) +
    geom_point(size = 0.05, alpha = 0.4, color = "gray80") +
    geom_point(data = pred.filter, color = "forestgreen", size = 0.05) +
    geom_smooth(method = "lm", alpha = 0.3, color = "gray60", aes(fill = "all points")) +
    geom_smooth(data = pred.filter, method = "lm", color = "forestgreen", aes(fill = "Tropical & subtropical\nmoist broadleaf forest"), alpha = 0.3) +
    scale_fill_manual("", values = cols) +
    scale_alpha_manual(values = c(0.2, 1)) +
    annotate(geom = "text", label = paste0("R\u00b2 = ", format(round(cor$estimate^2, 2), nsmall = 2)),  
             x = -Inf, y = Inf, hjust = -0.5, vjust = 3, color = "gray50") +
    annotate(geom = "text", label = paste0("R\u00b2 = ", format(round(cor.filter$estimate^2, 2), nsmall = 2)),  
             x = -Inf, y = Inf, hjust = -0.5, vjust = 5, color = "forestgreen") +
    scale_x_continuous("Mean Verticality") +
    scale_y_continuous(paste0(taxa, " Richness")) +
    theme_classic()
  
  return(list(sesvert, meanvert))
}

load("results/sdmTMB_models/amphibians_sesvert.RData")
amph = plot_cor(m.final, pred, "Amphibians")
amph_biome = plot_cor_biome(m.final, pred, "Amphibians")

load("results/sdmTMB_models/reptiles_sesvert.RData")
rept = plot_cor(m.final, pred, "Reptiles")
rept_biome = plot_cor_biome(m.final, pred, "Reptiles")

load("results/sdmTMB_models/birds_sesvert.RData")
birds = plot_cor(m.final, pred, "Birds")
birds_biome = plot_cor_biome(m.final, pred, "Birds")


load("results/sdmTMB_models/mammals_sesvert.RData")
mammals = plot_cor(m.final, pred, "Mammals")
mammals_biome = plot_cor_biome(m.final, pred, "Mammals")

design = "
AB
CD
EF
GH"

birds[[1]] + birds[[2]] +
  mammals[[1]] + mammals[[2]] +
  rept[[1]] + rept[[2]] +
  amph[[1]] + amph[[2]] + 
  plot_layout(design = design) &
  theme(legend.position = "none")

birds_biome[[1]] + birds_biome[[2]] +
  mammals_biome[[1]] + mammals_biome[[2]] +
  rept_biome[[1]] + rept_biome[[2]] +
  amph_biome[[1]] + amph_biome[[2]] + 
  plot_layout(design = design) &
  theme(legend.position = "none")











