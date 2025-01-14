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
# run simple linear model to look at correlation between vert and richness
cor_biome = function(mod, svcmap = NULL, pred, var, taxon) {
  p = ggplot(pred, aes(.data[[var]], y = rich, color = biome, fill = biome)) +
    geom_point(pch = ".", alpha = 0.5) +
    geom_smooth(method = "lm", alpha = 0.4) +
    theme_classic()
  
  # correlation with all points
  lm.all = lm(formula(paste0("rich ~ ", var)), data = pred)
  #cor.all = cor.test(pred[,var], pred[,"rich"])
  
  p2 = ggplot(pred, aes(.data[[var]], y = rich)) +
    geom_point(pch = ".", alpha = 0.5, color = "gray") +
    geom_smooth(method = "lm", color = "blue", fill = "blue") +
    scale_y_continuous(paste0(taxon, "\nrichness")) +
    scale_x_continuous(ifelse(var == "vert.mean.ses", "SES Verticality", "Mean Verticality")) +
    theme_classic()
  
  if(!is.null(svcmap)) {
    pred2 = pred %>% 
      mutate(x = x*1e5, y = y*1e5) %>% 
      mutate(svc_sig = terra::extract(svc$svc_sig, pred[,c("x", "y")])[,2])
    filter(svc_sig == 1 & svc_median_effect_link > 0)
    
  }
  
  # cor = pred %>% 
  #   group_by(biome) %>% 
  #   summarise(correlation = cor(.data[[var]], rich),
  #             pval = as.numeric(cor.test(.data[[var]], rich)["p.value"])) %>% 
  #   mutate(taxon = taxon, var = var,
  #          sig = ifelse(pval < 0.05, "sig", "notsig"))
  
  
  biomes = unique(pred$biome)
  
  lms = data.frame(n = numeric(length = length(biomes)),
                   slope = double(length = length(biomes)),
                   p = double(length = length(biomes)),
                   r2 = double(length = length(biomes)),
                   biome = character(length = length(biomes)),
                   taxa = character(length = length(biomes)),
                   var = character(length = length(biomes)))
  
  for(b in 1:length(biomes)) {
    d = pred %>% filter(biome == biomes[b])
    if(nrow(d) >= 30) {
      fit = lm(formula(paste0("rich ~ ", var)), data = d)
      summ = summary(fit)
      lms$n[b] = nrow(d)
      lms$slope[b] = summ$coefficients[2,1]
      lms$p[b] = summ$coefficients[2,4]
      lms$r2[b] = summ$adj.r.squared
      lms$biome[b] = as.character(biomes[b])
    } else {
      lms$n = nrow(d)
      lms$slope = NA
      lms$p = NA
      lms$r2 = NA
      lms$biome = as.character(biomes[b])
    }
  }
  lms$taxa = taxon
  lms$var = var
  
  # add model results including all biomes
  summall = summary(lm.all)
  summall = data.frame(n = nrow(pred),
                       slope = summall$coefficients[2,1],
                       p = summall$coefficients[2,4],
                       r2 = summall$adj.r.squared,
                       biome = "All Biomes", 
                       taxa = taxon,
                       var = var)
  lms = rbind(lms, summall)
  
  # all = data.frame(biome = "All Biomes", correlation = cor.all$estimate,
  #                  pval = cor.all$p.value,
  #                  taxon = taxon, var = var) %>% 
  #   mutate(sig = ifelse(pval < 0.05, "sig", "notsig"))
  # 
  # cor = rbind(cor, all)
  
  return(list(plot = p, plot2 = p2, lms = lms, pred = pred))
}


load("results/sdmTMB_models/amphibians_sesvert.RData")
#amph.sesvert.svc = rast("results/sdmTMB_models/svc_uncertainty/amph_sesvert.tif")
amph.corbiome.sesvert = cor_biome(mod, pred, svcmap = NULL, var = "vert.mean.ses", taxon = "Amphibians")

load("results/sdmTMB_models/amphibians_meanvert.RData")
amph.corbiome.meanvert = cor_biome(mod, pred, svcmap = NULL, "vert.mean", "Amphibians")

load("results/sdmTMB_models/amphibians_parb.RData")
amph.corbiome.parb = cor_biome(mod, pred, svcmap = NULL, "p.arb", "Amphibians")

load("results/sdmTMB_models/reptiles_sesvert.RData")
rept.corbiome.sesvert = cor_biome(mod, pred, svcmap = NULL, "vert.mean.ses", "Reptiles")

load("results/sdmTMB_models/reptiles_meanvert.RData")
rept.corbiome.meanvert = cor_biome(mod, pred, svcmap = NULL, "vert.mean", "Reptiles")

load("results/sdmTMB_models/reptiles_parb.RData")
rept.corbiome.parb = cor_biome(mod, pred, svcmap = NULL, "p.arb", "Reptiles")

load("results/sdmTMB_models/mammals_sesvert.RData")
mammals.corbiome.sesvert = cor_biome(mod, pred, svcmap = NULL, "vert.mean.ses", "Mammals")

load("results/sdmTMB_models/mammals_meanvert.RData")
mammals.corbiome.meanvert = cor_biome(mod, pred, svcmap = NULL, "vert.mean", "Mammals")

load("results/sdmTMB_models/mammals_parb.RData")
mammals.corbiome.parb = cor_biome(mod, pred, svcmap = NULL, "p.arb", "Mammals")

load("results/sdmTMB_models/birds_sesvert.RData")
birds.corbiome.sesvert = cor_biome(mod, pred, svcmap = NULL, "vert.mean.ses", "Birds")

load("results/sdmTMB_models/birds_meanvert.RData")
birds.corbiome.meanvert = cor_biome(mod, pred, svcmap = NULL, "vert.mean", "Birds")

load("results/sdmTMB_models/birds_parb.RData")
birds.corbiome.parb = cor_biome(mod, pred, svcmap = NULL, "p.arb", "Birds")

# cordf = rbind(
#   amph.corbiome.sesvert$cor,
#   amph.corbiome.meanvert$cor,
#   rept.corbiome.sesvert$cor,
#   rept.corbiome.meanvert$cor,
#   mammals.corbiome.sesvert$cor,
#   mammals.corbiome.meanvert$cor,
#   birds.corbiome.sesvert$cor,
#   birds.corbiome.meanvert$cor
# )

cordf = rbind(
  amph.corbiome.sesvert$lms,
  amph.corbiome.meanvert$lms,
  rept.corbiome.sesvert$lms,
  rept.corbiome.meanvert$lms,
  mammals.corbiome.sesvert$lms,
  mammals.corbiome.meanvert$lms,
  birds.corbiome.sesvert$lms,
  birds.corbiome.meanvert$lms
)


cordf %>% 
  mutate(taxa = factor(taxa, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")),
         biome = case_when(biome == "All Biomes" ~ "All Biomes",
                           biome == "Tundra" ~ "Tundra",
                           biome == "Boreal Forests/Taiga" ~ "Boreal Forest",
                           biome == "Deserts & Xeric Shrublands" ~ "Deserts & Xeric Shrublands",
                           biome == "Flooded Grasslands & Savannas" ~ "Flooded GS",
                           biome == "Mediterranean Forests, Woodlands & Scrub" ~ "Mediterranean FWS",
                           biome == "Montane Grasslands & Shrublands" ~ "Montane GS",
                           biome == "Temperate Broadleaf & Mixed Forests" ~ "Temperate BMF",
                           biome == "Temperate Conifer Forests" ~ "Temperate CF",
                           biome == "Temperate Grasslands, Savannas & Shrublands" ~ "Temperate GSS",
                           biome == "Tropical & Subtropical Coniferous Forests" ~ "Tropical CF",
                           biome == "Tropical & Subtropical Dry Broadleaf Forests" ~ "Tropical DBF",
                           biome == "Tropical & Subtropical Grasslands, Savannas & Shrublands" ~ "Tropical GSS",
                           biome == "Tropical & Subtropical Moist Broadleaf Forests" ~ "Tropical MBF"),
         biome = factor(biome, levels = c("All Biomes",
                                          "Tundra",
                                          "Boreal Forest",
                                          "Deserts & Xeric Shrublands",
                                          "Flooded GS",
                                          "Mediterranean FWS",
                                          "Montane GS",
                                          "Temperate BMF",
                                          "Temperate CF",
                                          "Temperate GSS",
                                          "Tropical CF",
                                          "Tropical DBF",
                                          "Tropical GSS",
                                          "Tropical MBF")),
        # var = ifelse(var == "vert.mean.ses", "SES Verticality", "Mean Verticality"),
         sig = ifelse(p > 0.05, "notsig", "sig"),
         fillcol = paste0(var, sig)) %>% 
  drop_na() %>% 
  pivot_wider(names_from = var, values_from = slope) %>% 
  ggplot(aes(fill = fillcol, color = fillcol,size = r2)) +
  annotate(geom = "rect", xmin = 1.5, xmax = 3.5, ymin = -Inf, ymax = Inf, fill = "gray80", alpha = 0.5) +
  annotate(geom = "rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf, fill = "purple", alpha = 0.5) +
  annotate(geom = "rect", xmin = 7.5, xmax = 10.5, ymin = -Inf, ymax = Inf, fill = "orange", alpha = 0.5) +
  annotate(geom = "rect", xmin = 10.5, xmax = 14.5, ymin = -Inf, ymax = Inf, fill = "forestgreen", alpha = 0.5) +
  geom_point(aes(x = biome, y = vert.mean.ses), pch = 21) +
  geom_point(aes(x = biome, y = vert.mean/100), pch = 21) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~taxa, ncol = 2, scales = "free_y") +
  #scale_color_manual(values = c("blue", "skyblue")) +
  scale_fill_manual(values = c("white", "blue", "white", "skyblue")) +
  scale_color_manual(values = c("blue", "blue", "skyblue", "skyblue")) +
  scale_size(range = c(0.5,4)) +
  scale_y_continuous(name = "Slope (SES verticality)", 
                     sec.axis = sec_axis(~ . * 100, name = "Slope (Mean verticality)")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 305, hjust = 0, vjust = 0.5),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 8),
        axis.title.y = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 2))

ggsave("figures/main_figs/vert_rich_lm_slopes.svg", width = 180, height = 130, dpi = 600, units = "mm")
ggsave("figures/main_figs/vert_rich_lm_slopes.png", width = 180, height = 130, dpi = 600, units = "mm")

# cordf %>% 
#   mutate(taxon = factor(taxon, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")),
#          biome = factor(biome, levels = c("All Biomes",
#                                           "Tundra",
#                                           "Boreal Forests/Taiga",
#                                           "Deserts & Xeric Shrublands",
#                                           "Flooded Grasslands & Savannas",
#                                           "Mediterranean Forests, Woodlands & Scrub",
#                                           "Montane Grasslands & Shrublands",
#                                           "Temperate Broadleaf & Mixed Forests",
#                                           "Temperate Conifer Forests",
#                                           "Temperate Grasslands, Savannas & Shrublands",
#                                           "Tropical & Subtropical Coniferous Forests",
#                                           "Tropical & Subtropical Dry Broadleaf Forests",
#                                           "Tropical & Subtropical Grasslands, Savannas & Shrublands",
#                                           "Tropical & Subtropical Moist Broadleaf Forests")),
#          var = ifelse(var == "vert.mean.ses", "SES Verticality", "Mean Verticality"),
#          fillcol = paste0(var, sig)) %>% 
#   drop_na() %>% 
#   ggplot(aes(x = biome, y = correlation, color = var, fill = fillcol)) +
#   geom_point(pch = 21, size = 1) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   facet_rep_wrap(~taxon, ncol = 2) +
#   scale_color_manual(values = c("blue", "skyblue")) +
#   scale_fill_manual(values = c("white", "blue", "white", "skyblue")) +
#   scale_y_continuous("Pearson correlation coefficient (r)") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 305, hjust = 0, vjust = 0.5),
#         axis.title.x = element_blank(),
#         axis.text = element_text(size = 8),
#         axis.title.y = element_text(size = 10),
#         legend.title = element_blank())
# 
# ggsave("figures/main_figs/vert_rich_cor_coefs.png", width = 180, height = 120, dpi = 600, units = "mm")
# ggsave("figures/main_figs/vert_rich_cor_coefs.svg", width = 180, height = 120, dpi = 600, units = "mm")
# 
# 
# 

# scatter plots for supplement


birds.corbiome.sesvert$pred$taxa = "Birds"
birds.corbiome.sesvert$var = "SES Verticality"
birds.corbiome.meanvert$pred$taxa = "Birds"
birds.corbiome.meanvert$var = "Mean Verticality"

mammals.corbiome.sesvert$pred$taxa = "Mammals"
mammals.corbiome.sesvert$var = "SES Verticality"
mammals.corbiome.meanvert$pred$taxa = "Mammals"
mammals.corbiome.meanvert$var = "Mean Verticality"

rept.corbiome.sesvert$pred$taxa = "Reptiles"
rept.corbiome.sesvert$var = "SES Verticality"
rept.corbiome.meanvert$pred$taxa = "Reptiles"
rept.corbiome.meanvert$var = "Mean Verticality"

amph.corbiome.sesvert$pred$taxa = "Amphibians"
amph.corbiome.sesvert$var = "SES Verticality"
amph.corbiome.meanvert$pred$taxa = "Amphibians"
amph.corbiome.meanvert$var = "Mean Verticality"

d = bind_rows(birds.corbiome.meanvert$pred,
          mammals.corbiome.meanvert$pred,
          rept.corbiome.meanvert$pred, 
          amph.corbiome.meanvert$pred) %>% 
  dplyr::select(rich, vert.mean, vert.mean.ses, taxa, biome) %>% 
  pivot_longer(cols = c("vert.mean", "vert.mean.ses"), names_to = "var", values_to = "vert") %>% 
  mutate(taxa = factor(taxa, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")),
         var = ifelse(var == "vert.mean", "Mean Verticality", "SES Verticality"))

plt = d %>% 
  ggplot(aes(x = vert, y = rich, color = biome, fill = biome)) +
  geom_point(pch = ".") +
  geom_smooth(method = "lm", alpha = 0.5, linewidth = 0.5) +
  facet_grid(rows = vars(taxa), cols = vars(var), scales = "free") +
  scale_x_continuous("Verticality") +
  scale_y_continuous("Richness") +
  theme_bw() +
  theme(panel.background = element_rect(color = "black", fill = NA),
        legend.text = element_text(size = 6),
        legend.title = element_blank())

ggsave("figures/supp_figs/vert_richness_correlations.png", width = 180, height = 140, units = "mm", dpi = 600)

