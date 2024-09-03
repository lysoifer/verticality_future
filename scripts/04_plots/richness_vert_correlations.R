# plot verticality x richness for regions with canopy_coef > 1 and canopy_coef < 1

library(ggplot2)
library(tidyverse)
library(patchwork)
library(terra)
library(sdmTMB)
library(colorspace)
library(lemon)

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
    facet_grid(cols = vars(biome),rows = vars(taxa), scales = "free") +
    theme_classic()

dat %>% 
  filter(rich >= 5) %>% 
  ggplot(aes(x = vert.mean.ses, y = rich, color = biome)) +
  geom_point(pch = ".") +
  geom_smooth(method = "lm") +
  facet_wrap(~taxa, scales = "free") +
  theme_classic()


plot_cor = function(mod, pred, svc, taxa, var, xlab, ylab) {
  # mod = sdmTMB model outpum
  # pred = prediction to current climate conditions
  # svc = spatraster of svc uncertainty from plot_spatial_varying.R
  # taxa = character to label plot
  # var = either "vert.mean" or "vert.mean.ses"
  # xlab = x-axis label
  
  coefs = tidy(mod)
  

  pred.r = pred %>% 
    dplyr::select(x,y,var, rich) %>% 
    mutate(x = x*1e5, y = y*1e5) %>% 
    rast(crs = "+proj=cea +datum=WGS84")
  
  # combine vert, richness, and svc significance in a dataframe
  df = c(pred.r, svc$svc_sig, svc$svc_median_effect_link) %>% 
    as.data.frame()

  # df.sig = df %>% 
  #   filter(svc_sig == 1)
  
  df.sigpos = df %>% 
    filter(svc_sig == 1 & svc_median_effect_link > 0)
  
  cor.all = cor.test(df[,var], df$rich) # correlation between vert and richness with all points
  #cor.sig = cor.test(df.sig[,var], df.sig$rich) # correlation between ver and richness with all points that have a significant canopy height coefficient
  cor.sigpos = cor.test(df.sigpos[,var], df.sigpos$rich)
  
  df$part = "all points"
  #df.sig$part = "CH coef significant"
  df.sigpos$part = "CH coef significant\n& positive"
  
  # df = rbind(df, df.sig, df.sigpos)
  # 
  # df.cor = data.frame(cor = c(round(cor.all$estimate, 2), round(cor.sig$estimate, 2), round(cor.sigpos$estimate, 2)),
  #                     part = c("all", "sig", "sigpos"))
  
  df = rbind(df, df.sigpos)
  
  df.cor = data.frame(cor = c(round(cor.all$estimate, 2), round(cor.sigpos$estimate, 2)),
                      part = c("all", "sigpos"))
  
  
  plot = ggplot(df, aes(x = .data[[var]], y = rich, color = part, fill = part)) +
    geom_point(pch = ".") +
    geom_smooth(method = "lm") +
    scale_color_manual(values = sequential_hcl(5, "Mint")[c(4,2,1)]) +
    scale_fill_manual(values = sequential_hcl(5, "Mint")[c(4,2,1)]) +
    scale_x_continuous(xlab) +
    scale_y_continuous(ylab, limits = c(5,max(df$rich))) +
    annotate(geom = "text", label = paste0("R\u00b2 = ", format(round(cor.all$estimate^2, 2), nsmall = 2)),  
             x = -Inf, y = Inf, hjust = -0.25, vjust = 1, color = sequential_hcl(5, "Mint")[4], size = 2.5, fontface = "bold") +
    # annotate(geom = "text", label = paste0("R\u00b2 = ", format(round(cor.sig$estimate^2, 2), nsmall = 2)),  
    #          x = -Inf, y = Inf, hjust = -0.25, vjust = 3, color = sequential_hcl(5, "Mint")[2], size = 2.5, fontface = "bold") +
    annotate(geom = "text", label = paste0("R\u00b2 = ", format(round(cor.sigpos$estimate^2, 2), nsmall = 2)),  
             x = -Inf, y = Inf, hjust = -0.25, vjust = 3, color = sequential_hcl(5, "Mint")[1], size = 2.5, fontface = "bold") +
    theme_classic() +
    theme(legend.position = "bottom",
          legend.title = element_blank())

  return(list(plot = plot, cor.all = cor.all, cor.sigpos = cor.sigpos))
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
amph.sesvert.svc = rast("results/sdmTMB_models/svc_uncertainty/amph_sesvert.tif")
amph.sesvert = plot_cor(mod, pred, svc = amph.sesvert.svc, taxa = "Amphibians", var = "vert.mean.ses", xlab = "SES Verticality", ylab = "Amphibain\nrichness")
#amph_biome = plot_cor_biome(mod, pred, "Amphibians")

load("results/sdmTMB_models/amphibians_meanvert.RData")
amph.meanvert.svc = rast("results/sdmTMB_models/svc_uncertainty/amph_meanvert.tif")
amph.meanvert = plot_cor(mod, pred, svc = amph.meanvert.svc, taxa = "Amphibians", var = "vert.mean", xlab = "Mean Verticality", ylab = "Amphibain\nrichness")
#amph_biome = plot_cor_biome(mod, pred, "Amphibians")

load("results/sdmTMB_models/reptiles_sesvert.RData")
rept.sesvert.svc = rast("results/sdmTMB_models/svc_uncertainty/rept_sesvert.tif")
rept.sesvert = plot_cor(mod, pred, svc = rept.sesvert.svc, taxa = "Reptiles", var = "vert.mean.ses", xlab = "SES Verticality", ylab = "Reptile\nrichness")

load("results/sdmTMB_models/reptiles_meanvert.RData")
rept.meanvert.svc = rast("results/sdmTMB_models/svc_uncertainty/rept_meanvert.tif")
rept.meanvert = plot_cor(mod, pred, svc = rept.meanvert.svc, taxa = "Reptiles", var = "vert.mean", xlab = "Mean Verticality", ylab = "Reptile\nrichness")

load("results/sdmTMB_models/mammals_sesvert.RData")
mammals.sesvert.svc = rast("results/sdmTMB_models/svc_uncertainty/mammals_sesvert.tif")
mammals.sesvert = plot_cor(mod, pred, svc = mammals.sesvert.svc, taxa = "Mammals", var = "vert.mean.ses", xlab = "SES Verticality", ylab = "Mammal\nrichness")

load("results/sdmTMB_models/mammals_meanvert.RData")
mammals.meanvert.svc = rast("results/sdmTMB_models/svc_uncertainty/mammals_meanvert.tif")
mammals.meanvert = plot_cor(mod, pred, svc = mammals.meanvert.svc, taxa = "Mammals", var = "vert.mean", xlab = "Mean Verticality", ylab = "Mammal\nrichness")

load("results/sdmTMB_models/birds_sesvert.RData")
birds.sesvert.svc = rast("results/sdmTMB_models/svc_uncertainty/birds_sesvert.tif")
birds.sesvert = plot_cor(mod, pred, svc = birds.sesvert.svc, taxa = "Birds", var = "vert.mean.ses", xlab = "SES Verticality", ylab = "Bird\nrichness")

load("results/sdmTMB_models/birds_meanvert.RData")
birds.meanvert.svc = rast("results/sdmTMB_models/svc_uncertainty/birds_meanvert.tif")
birds.meanvert = plot_cor(mod, pred, svc = birds.meanvert.svc, taxa = "Birds", var = "vert.mean", xlab = "Mean Verticality", ylab = "Bird\nrichness")


library(grid)
library(patchwork)
col1 = wrap_elements(panel = textGrob("SES Mean Verticality"))
col2 = wrap_elements(panel = textGrob("Mean Verticality"))

row1 = wrap_elements(panel = textGrob("Birds", rot = 90))
row2 = wrap_elements(panel = textGrob("Mammals", rot = 90))
row3 = wrap_elements(panel = textGrob("Reptiles", rot = 90))
row4 = wrap_elements(panel = textGrob("Amphibians", rot = 90))

design = "
AE
BF
CG
DH"

# plot areas where canopy height coef is significant (non-significant areas are dark gray) - link space
p = 
  birds.sesvert$plot + mammals.sesvert$plot + rept.sesvert$plot + amph.sesvert$plot +
  birds.meanvert$plot + mammals.meanvert$plot + rept.meanvert$plot + amph.meanvert$plot + 
  plot_layout(design = design, heights = c(1,1,1,1), widths = c(1,1), guides = "collect", axis_titles = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag.position = c(0.95, 0.95)) &
  theme(legend.position = "bottom")

png("figures/main_figs/vert_rich_cor.png", width = 150, height = 200, res = 300, units = "mm")
p
dev.off()




# CORRELATIONS BY BIOME ---------------------------------------------------
write.excel <- function(x,row.names=FALSE,col.names=TRUE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)
}

# plot correlations of richness by verticality by biome
cor_biome = function(mod, svc = NULL, pred, var, taxon) {
  p = ggplot(pred, aes(.data[[var]], y = rich, color = biome, fill = biome)) +
    geom_point(pch = ".", alpha = 0.5) +
    geom_smooth(method = "lm", alpha = 0.4) +
    theme_classic()
  
  if(!is.null(svc)) {
    pred2 = pred %>% 
      mutate(x = x*1e5, y = y*1e5) %>% 
      mutate(svc_sig = terra::extract(svc$svc_sig, pred[,c("x", "y")])[,2])
    filter(svc_sig == 1 & svc_median_effect_link > 0)
    
  }
  
  cor = pred %>% 
    group_by(biome) %>% 
    summarise(correlation = cor(.data[[var]], rich)) %>% 
    mutate(taxon = taxon, var = var)
  
  return(list(plot = p, cor = cor))
}


load("results/sdmTMB_models/amphibians_sesvert.RData")
amph.sesvert.svc = rast("results/sdmTMB_models/svc_uncertainty/amph_sesvert.tif")
amph.corbiome.sesvert = cor_biome(mod, pred, "vert.mean.ses", "Amphibians")

load("results/sdmTMB_models/amphibians_meanvert.RData")
amph.corbiome.meanvert = cor_biome(mod, pred, "vert.mean", "Amphibians")

load("results/sdmTMB_models/amphibians_parb.RData")
amph.corbiome.parb = cor_biome(mod, pred, "p.arb", "Amphibians")

load("results/sdmTMB_models/reptiles_sesvert.RData")
rept.corbiome.sesvert = cor_biome(mod, pred, "vert.mean.ses", "Reptiles")

load("results/sdmTMB_models/reptiles_meanvert.RData")
rept.corbiome.meanvert = cor_biome(mod, pred, "vert.mean", "Reptiles")

load("results/sdmTMB_models/reptiles_parb.RData")
rept.corbiome.parb = cor_biome(mod, pred, "p.arb", "Reptiles")

load("results/sdmTMB_models/mammals_sesvert.RData")
mammals.corbiome.sesvert = cor_biome(mod, pred, "vert.mean.ses", "Mammals")

load("results/sdmTMB_models/mammals_meanvert.RData")
mammals.corbiome.meanvert = cor_biome(mod, pred, "vert.mean", "Mammals")

load("results/sdmTMB_models/mammals_parb.RData")
mammals.corbiome.parb = cor_biome(mod, pred, "p.arb", "Mammals")

load("results/sdmTMB_models/birds_sesvert.RData")
birds.corbiome.sesvert = cor_biome(mod, pred, "vert.mean.ses", "Birds")

load("results/sdmTMB_models/birds_meanvert.RData")
birds.corbiome.meanvert = cor_biome(mod, pred, "vert.mean", "Birds")

load("results/sdmTMB_models/birds_parb.RData")
birds.corbiome.parb = cor_biome(mod, pred, "p.arb", "Birds")

cordf = rbind(
  amph.corbiome.sesvert$cor,
  amph.corbiome.meanvert$cor,
  amph.corbiome.parb$cor,
  rept.corbiome.sesvert$cor,
  rept.corbiome.meanvert$cor,
  rept.corbiome.parb$cor,
  mammals.corbiome.sesvert$cor,
  mammals.corbiome.meanvert$cor,
  mammals.corbiome.parb$cor,
  birds.corbiome.sesvert$cor,
  birds.corbiome.meanvert$cor,
  birds.corbiome.parb$cor
)

ggplot(cordf, aes(x = biome, y = correlation, color = var)) +
  geom_point(size = 2) +
  facet_rep_wrap(~taxon) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 305, hjust = 0, vjust = 0.5))



