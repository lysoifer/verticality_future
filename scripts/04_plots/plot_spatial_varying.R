library(terra)
library(ggplot2)
library(colorspace)
library(tidyverse)
library(patchwork)
library(cowplot)
library(tidyterra)
library(sdmTMB)

wd = vect("data/original/rnaturalearth_world.shp") %>% 
  project("epsg:4326") %>% 
  crop(ext(-180,180,-60,83.64))

plot_spatial_varying = function(mod, pred, var, v, legend_title) {
  coefs = tidy(mod)
  coef.var = as.numeric(coefs[which(coefs$term == var), "estimate"])
  #pred = predict(mod, type = "response")
  dat = pred %>% 
    mutate(x = x*1e5, y = y*1e5,
           canopy_height_est_slope = zeta_s_canopy_height + coef.var) %>% 
    relocate(x, y, .before = rich) %>% 
    rast(crs = "+proj=cea + datum=WGS84") %>% 
    project("epsg:4326")
    
    
  ggplot() +
    geom_spatvector(data = v, color = "black", fill = "gray70") +
    geom_spatraster(data = dat, aes(fill = canopy_height_est_slope)) +
    geom_spatvector(data = v, color = "black", fill = NA) +
    scale_fill_continuous_divergingx("BrBG", na.value = NA, guide = guide_colorbar(legend_title)) +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(color = "black", fill = NA),
          panel.grid = element_blank())

}

# plot canopy height coefficienby canopy height standard deviation
plot_coef_height = function(pred) {
  height = pred %>% 
    dplyr::select(x,y,canopy_height) %>% 
    mutate(x = x*1e5, y = y*1e5) %>% 
    rast(crs = "+proj=cea +datum=WGS84") %>% 
    focal(w = 5, fun = sd, na.rm = T) %>% 
    as.data.frame(xy = T) %>% 
    mutate(x = round(x/1e5,5), y = round(y/1e5,5)) %>% 
    inner_join(pred, by = c("x", "y")) %>% 
    dplyr::select(x,y,biome, ecoregion, biorealm, zeta_s_canopy_height, focal_sd) %>% 
    rename(canopy_height_sd = focal_sd)
  
  p = height %>% 
    ggplot(aes(x = canopy_height_sd, y = zeta_s_canopy_height)) +
    geom_point(pch = ".") +
    geom_smooth(method = "lm") +
    scale_x_continuous("SD canopy height", breaks = seq(0,1.25,0.25)) +
    scale_y_continuous("Canopy height coefficient") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p)
    
}



# Amphibians

load("results/sdmTMB_models/amphibians_sesvert.RData")
amph.sesvert = plot_spatial_varying(mod = mod, pred = pred, var = "canopy_height", v = wd, legend_title = "coefficient")
coef.a.ses = plot_coef_height(pred)

load("results/sdmTMB_models/amphibians_meanvert.RData")
amph.meanvert = plot_spatial_varying(mod = mod, pred = pred, var = "canopy_height", v = wd, legend_title = "coefficient")
coef.a.mean = plot_coef_height(pred)


# Reptiles

load("results/sdmTMB_models/reptiles_sesvert.RData")
rept.sesvert = plot_spatial_varying(mod = mod, pred = pred, var = "canopy_height", v = wd, legend_title = "coefficient")
coef.r.ses = plot_coef_height(pred)

load("results/sdmTMB_models/reptiles_meanvert.RData")
rept.meanvert.svc = plot_spatial_varying(mod = mod, pred = pred, var = "canopy_height", v = wd, legend_title = "coefficient")
coef.r.mean = plot_coef_height(pred)


# Mammals

load("results/sdmTMB_models/mammals_sesvert.RData")
mammals.sesvert = plot_spatial_varying(mod = mod, pred = pred, var = "canopy_height", v = wd, legend_title = "coefficient")
coef.m.ses = plot_coef_height(pred)

load("results/sdmTMB_models/mammals_meanvert.RData")
mammals.meanvert = plot_spatial_varying(mod = mod, pred = pred, var = "canopy_height", v = wd, legend_title = "coefficient")
coef.m.mean = plot_coef_height(pred)


# Birds

load("results/sdmTMB_models/birds_sesvert.RData")
birds.sesvert = plot_spatial_varying(mod = mod, pred = pred, var = "canopy_height", v = wd, legend_title = "coefficient")
coef.b.ses = plot_coef_height(pred)

load("results/sdmTMB_models/birds_meanvert.RData")
birds.meanvert = plot_spatial_varying(mod = mod, pred = pred, var = "canopy_height", v = wd, legend_title = "coefficient")
coef.b.mean = plot_coef_height(pred)


# Full figure

design = "
AF
BG
CH
DI
EJ"

library(grid)
col1 = wrap_elements(panel = textGrob("SES Mean Verticality"))
col2 = wrap_elements(panel = textGrob("Mean Verticality"))

p = col1 + birds.sesvert+ mammals.sesvert + rept.sesvert + amph.sesvert +
  col2 + birds.meanvert + mammals.meanvert + rept.meanvert + amph.meanvert + 
  plot_layout(design = design, heights = c(0.01, -1, -1, -1, -1)) +
  plot_annotation(tag_levels = list(c("", "A", "B", "C", "D", "", "E", "F", "G", "H")))

png("figures/supp_figs/canopy_height_coefs.png", height = 300, width = 300, res = 300, units = "mm")
p
dev.off()

row1 = wrap_elements(panel = textGrob("Birds", rot = 90))
row2 = wrap_elements(panel = textGrob("Mammals", rot = 90))
row3 = wrap_elements(panel = textGrob("Reptiles", rot = 90))
row4 = wrap_elements(panel = textGrob("Amphibians", rot = 90))

design = "
AFK
BGL
CHM
DIN
EJO"

p = 
  plot_spacer() + row1 + row2 + row3 + row4 +
  col1 + coef.b.ses + coef.m.ses + coef.r.ses + coef.a.ses +
  col2 + coef.b.mean + coef.m.mean + coef.r.mean + coef.a.mean + 
  plot_layout(design = design, heights = c(0.2,1,1,1,1), widths = c(0.2,1,1)) +
  plot_annotation(tag_levels = list(c("","","","","", "A", "B", "C", "D", "", "E", "F", "G", "H")))

p[[12]] <- p[[12]] + theme(axis.title = element_blank())
p[[13]] <- p[[13]] + theme(axis.title = element_blank())
p[[14]] <- p[[14]] + theme(axis.title = element_blank())
p[[15]] <- p[[15]] + theme(axis.title.y = element_blank())
p[[7]] <- p[[7]] + theme(axis.title.x = element_blank())
p[[8]] <- p[[8]] + theme(axis.title.x = element_blank())
p[[9]] <- p[[9]] + theme(axis.title.x = element_blank())

png("figures/supp_figs/canopyHeightCoefs-canopyHeightSD.png", height = 200, width = 200, res = 300, units = "mm")
p
dev.off()



# SPATIAL VARYING UNCERTAINTY ---------------------------------------------

# see snowy owl example from sdmTMB paper
library(foreach)
library(doParallel)
library(ggnewscale)

# Calculate uncertainty of spatially varying canopy height coefficient
# see https://www.biorxiv.org/content/10.1101/2022.03.24.485545v3 (owl example for SVCs)

svc_uncertainty = function(mod, pred,legend_title) {
  cl = makeCluster(6)
  registerDoParallel(cl)
  
  tic()
  zeta_s = foreach(i = 1:6, .combine = "cbind", .packages = c("sdmTMB")) %dopar% {
    sdmTMB:::predict.sdmTMB(mod, nsim = 50, sims_var = "zeta_s")
  }
  stopCluster(cl)
  toc()
  
  sims = spread_sims(mod, nsim = 150)
  combined = sims$canopy_height + t(zeta_s)
  
  # Transformed out of link space
  svc_effect_invlink = mod$family$linkinv(apply(combined, 2, median))
  svc_lwr_invlink = mod$family$linkinv(apply(combined, 2, quantile, probs = 0.025))
  svc_upr_invlink = mod$family$linkinv(apply(combined, 2, quantile, prob = 0.975))
  
  # Here we do NOT apply inverse link function - so coefficients are in link space (for beta models, this is logit link space)
  svc_effect_link = apply(combined, 2, median)
  svc_lwr_link = apply(combined, 2, quantile, probs = 0.025)
  svc_upr_link = apply(combined, 2, quantile, prob = 0.975)
  
  pred$svc_median_effect_link = svc_effect_link
  pred$svc_lwr_link = svc_lwr_link
  pred$svc_upr_link = svc_upr_link
  pred = pred %>% 
    mutate(svc_sig = ifelse((svc_lwr_link < 0 & svc_upr_link < 0) | (svc_lwr_link > 0 & svc_upr_link > 0), 1, 0),
           svc_sig = factor(svc_sig))
  
  pred$svc_median_effect_invlink = svc_effect_invlink
  pred$svc_lwr_invlink = svc_lwr_invlink
  pred$svc_upr_invlink = svc_upr_invlink
  
  dat = pred %>% 
    mutate(x = x*1e5, y = y*1e5) %>% 
    relocate(x, y, .before = rich) %>% 
    dplyr::select(x, y, svc_median_effect_link, svc_lwr_link, svc_upr_link, svc_median_effect_invlink, svc_lwr_invlink, svc_upr_invlink, svc_sig) %>% 
    rast(crs = "+proj=cea + datum=WGS84") %>% 
    project("epsg:4326")
  
  return(dat)
  
}

load("results/sdmTMB_models/amphibians_sesvert.RData")
amph.sesvert = svc_uncertainty(mod = mod, pred = pred, legend_title = "coefficient")
writeRaster(amph.sesvert, file = "results/sdmTMB_models/svc_uncertainty/amph_sesvert.tif")


load("results/sdmTMB_models/amphibians_meanvert.RData")
amph.meanvert = svc_uncertainty(mod = mod, pred = pred, legend_title = "coefficient")
writeRaster(amph.meanvert, file = "results/sdmTMB_models/svc_uncertainty/amph_meanvert.tif")

# Reptiles

load("results/sdmTMB_models/reptiles_sesvert.RData")
rept.sesvert = svc_uncertainty(mod = mod, pred = pred, legend_title = "coefficient")
writeRaster(rept.sesvert, file = "results/sdmTMB_models/svc_uncertainty/rept_sesvert.tif")


load("results/sdmTMB_models/reptiles_meanvert.RData")
rept.meanvert = svc_uncertainty(mod = mod, pred = pred, legend_title = "coefficient")
writeRaster(rept.meanvert, file = "results/sdmTMB_models/svc_uncertainty/rept_meanvert.tif")


# Mammals

load("results/sdmTMB_models/mammals_sesvert.RData")
mammals.sesvert = svc_uncertainty(mod = mod, pred = pred, legend_title = "coefficient")
writeRaster(mammals.sesvert, file = "results/sdmTMB_models/svc_uncertainty/mammals_sesvert.tif")


load("results/sdmTMB_models/mammals_meanvert.RData")
mammals.meanvert = svc_uncertainty(mod = mod, pred = pred, legend_title = "coefficient")
writeRaster(mammals.meanvert, file = "results/sdmTMB_models/svc_uncertainty/mammals_meanvert.tif")


# Birds

load("results/sdmTMB_models/birds_sesvert.RData")
birds.sesvert = svc_uncertainty(mod = mod, pred = pred, legend_title = "coefficient")
writeRaster(birds.sesvert, file = "results/sdmTMB_models/svc_uncertainty/birds_sesvert.tif")


load("results/sdmTMB_models/birds_meanvert.RData")
birds.meanvert = svc_uncertainty(mod = mod, pred = pred, legend_title = "coefficient")
writeRaster(birds.meanvert, file = "results/sdmTMB_models/svc_uncertainty/birds_meanvert.tif")

plot_svc_uncertainty = function(dat, v) {
  
  pred = as.data.frame(dat)
  
  min_link = min(pred$svc_median_effect_link, pred$svc_lwr_link, pred$svc_upr_link)
  max_link = max(pred$svc_median_effect_link, pred$svc_lwr_link, pred$svc_upr_link)
  
  min_invlink = min(pred$svc_median_effect_invlink, pred$svc_lwr_invlink, pred$svc_upr_invlink)
  max_invlink = max(pred$svc_median_effect_invlink, pred$svc_lwr_invlink, pred$svc_upr_invlink)
  
  thm = theme(legend.position = "inside",
              legend.position.inside = c(0.03,0.7),
              legend.justification = c(0,1),
              legend.margin = margin(0,0,0.2,0.3),
              #legend.key = element_blank(),
              legend.background = element_blank(),
              legend.key.height = unit(3, units = "mm"),
              legend.key.width = unit(2, units = "mm"),
              legend.box.margin = margin(0,0,0,0),
              legend.box.spacing = unit(1,"mm"),
              legend.box.background = element_blank(),
              plot.background = element_blank(),
              plot.margin = unit(c(0,0,0,0), units = "mm"),
              panel.spacing = unit(0,"mm"),
              panel.background = element_rect(color = "black", fill = NA),
              #legend.title.position = "top",
              #legend.title = element_text(hjust = 0.5, size = 8),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank(),
              panel.grid = element_blank())
  
  
  median_link = ggplot() +
    geom_spatvector(data = v, color = "black", fill = "gray70") +
    geom_tile(data = dat, aes(x = x, y = y, fill = svc_median_effect_link)) +
    geom_spatvector(data = v, color = "black", fill = NA) +
    coord_sf(crs = 4326) +
    scale_fill_continuous_divergingx("BrBG", na.value = NA, guide = guide_colorbar(""), limits = c(min_link, max_link)) +
    thm
  
    # theme(axis.text = element_blank(),
    #       axis.title = element_blank(),
    #       axis.ticks = element_blank(),
    #       panel.background = element_rect(color = "black", fill = NA),
    #       panel.grid = element_blank(),
    #       legend.position = "inside",
    #       legend.position.inside = c(0.1,0.3),
    #       plot.margin = unit(c(0,0,0,0), units = "mm"))
  
  lwr_link = ggplot() +
    geom_spatvector(data = v, color = "black", fill = "gray70") +
    geom_tile(data = dat, aes(x = x, y = y, fill = svc_lwr_link)) +
    geom_spatvector(data = v, color = "black", fill = NA) +
    coord_sf(crs = 4326) +
    scale_fill_continuous_divergingx("BrBG", na.value = NA, guide = guide_colorbar(""), limits = c(min_link, max_link)) +
    thm
  
  upr_link = ggplot() +
    geom_spatvector(data = v, color = "black", fill = "gray70") +
    geom_tile(data = dat, aes(x = x, y = y, fill = svc_upr_link)) +
    geom_spatvector(data = v, color = "black", fill = NA) +
    coord_sf(crs = 4326) +
    scale_fill_continuous_divergingx("BrBG", na.value = NA, guide = guide_colorbar(""), limits = c(min_link, max_link)) +
    thm
  
  sig_link = ggplot() +
    geom_spatvector(data = v, color = "black", fill = "gray80") +
    geom_tile(data = dat, aes(x = x, y = y, fill = svc_median_effect_link)) +
    scale_fill_continuous_divergingx("BrBG", na.value = NA, guide = guide_colorbar(""), limits = c(min_link, max_link)) +
    ggnewscale::new_scale_fill() +
    #geom_tile(data = dat, aes(x = x, y = y, fill = svc_sig)) +
    geom_spatraster(data = as.factor(dat$svc_sig), aes(fill = svc_sig)) +
    geom_spatvector(data = v, color = "black", fill = NA) +
    scale_fill_manual(values = c("gray50", NA), na.value = NA, guide = guide_colorbar(show = F)) +
    coord_sf(crs = 4326) +
    thm
  
  if(mod$family$link == "identity") {mid = 0}else{mid = 0.5}
  median_invlink = ggplot() +
    geom_spatvector(data = v, color = "black", fill = "gray70") +
    geom_tile(data = dat, aes(x = x, y = y, fill = svc_median_effect_invlink)) +
    geom_spatvector(data = v, color = "black", fill = NA) +
    coord_sf(crs = 4326) +
    scale_fill_continuous_divergingx("BrBG", na.value = NA, guide = guide_colorbar(""), 
                                     limits = c(min_invlink, max_invlink), mid = mid) +
    thm

  lwr_invlink = ggplot() +
    geom_spatvector(data = v, color = "black", fill = "gray70") +
    geom_tile(data = dat, aes(x = x, y = y, fill = svc_lwr_invlink)) +
    geom_spatvector(data = v, color = "black", fill = NA) +
    coord_sf(crs = 4326) +
    scale_fill_continuous_divergingx("BrBG", na.value = NA, guide = guide_colorbar(""), 
                                     limits = c(min_invlink, max_invlink), mid = mid) +
    thm
  
  upr_invlink = ggplot() +
    geom_spatvector(data = v, color = "black", fill = "gray70") +
    geom_tile(data = dat, aes(x = x, y = y, fill = svc_upr_invlink)) +
    geom_spatvector(data = v, color = "black", fill = NA) +
    coord_sf(crs = 4326) +
    scale_fill_continuous_divergingx("BrBG", na.value = NA, guide = guide_colorbar(""), 
                                     limits = c(min_invlink, max_invlink), mid = mid) +
    thm
  
  sig_invlink = ggplot() +
    geom_spatvector(data = v, color = "black", fill = "gray80") +
    geom_tile(data = dat, aes(x = x, y = y, fill = svc_median_effect_invlink)) +
    scale_fill_continuous_divergingx("BrBG", na.value = NA, guide = guide_colorbar(""), 
                                     limits = c(min_invlink, max_invlink), mid = mid,
                                     position = "bottom") +
    ggnewscale::new_scale_fill() +
    #geom_tile(data = dat, aes(x = x, y = y, fill = svc_sig)) +
    geom_spatraster(data = as.factor(dat$svc_sig), aes(fill = svc_sig)) +
    geom_spatvector(data = v, color = "black", fill = NA) +
    scale_fill_manual(values = c("gray50", NA), na.value = NA, guide = guide_colorbar(show = F)) +
    coord_sf(crs = 4326) +
    thm
  
  # link is in linkspace
  # invlink is transformed out of link space
  return(list(median_link = median_link, lwr_link = lwr_link, upr_link = upr_link, sig_link = sig_link,
              median_invlink = median_invlink, lwr_invlink = lwr_invlink, upr_invlink = upr_invlink, sig_invlink = sig_invlink))
  
}

# make plots

amph.sesvert.plot = plot_svc_uncertainty(amph.sesvert, v = wd)
amph.meanvert.plot = plot_svc_uncertainty(amph.meanvert, v = wd)

rept.sesvert.plot = plot_svc_uncertainty(rept.sesvert, v = wd)
rept.meanvert.plot = plot_svc_uncertainty(rept.meanvert, v = wd)

mammals.sesvert.plot = plot_svc_uncertainty(mammals.sesvert, v = wd)
mammals.meanvert.plot = plot_svc_uncertainty(mammals.meanvert, v = wd)

birds.sesvert.plot = plot_svc_uncertainty(birds.sesvert, v = wd)
birds.meanvert.plot = plot_svc_uncertainty(birds.meanvert, v = wd)


library(grid)
library(patchwork)
col1 = wrap_elements(panel = textGrob("SES Mean Verticality"))
col2 = wrap_elements(panel = textGrob("Mean Verticality"))

row1 = wrap_elements(panel = textGrob("Birds", rot = 90))
row2 = wrap_elements(panel = textGrob("Mammals", rot = 90))
row3 = wrap_elements(panel = textGrob("Reptiles", rot = 90))
row4 = wrap_elements(panel = textGrob("Amphibians", rot = 90))

design = "
AFK
BGL
CHM
DIN
EJO"

# plot areas where canopy height coef is significant (non-significant areas are dark gray) - link space
p = 
  plot_spacer() + row1 + row2 + row3 + row4 +
  col1 + birds.sesvert.plot$sig_link + mammals.sesvert.plot$sig_link + rept.sesvert.plot$sig_link + amph.sesvert.plot$sig_link +
  col2 + birds.meanvert.plot$sig_link + mammals.meanvert.plot$sig_link + rept.meanvert.plot$sig_link + amph.meanvert.plot$sig_link + 
  plot_layout(design = design, heights = c(0.2,1,1,1,1), widths = c(0.2,1,1)) +
  plot_annotation(tag_levels = list(c("","","","","", "A", "B", "C", "D", "", "E", "F", "G", "H"))) &
  theme(plot.tag.position = c(0.025, 0.9))

png("figures/supp_figs/svc_canopyHeight_sig_link.png", height = 180, width = 200, res = 300, units = "mm")
p
dev.off()


# plot lower confidence interval of canopy height coef
p = 
  plot_spacer() + row1 + row2 + row3 + row4 +
  col1 + birds.sesvert.plot$lwr_link + mammals.sesvert.plot$lwr_link + rept.sesvert.plot$lwr_link + amph.sesvert.plot$lwr_link +
  col2 + birds.meanvert.plot$lwr_link + mammals.meanvert.plot$lwr_link + rept.meanvert.plot$lwr_link + amph.meanvert.plot$lwr_link + 
  plot_layout(design = design, heights = c(0.2,1,1,1,1), widths = c(0.2,1,1)) +
  plot_annotation(tag_levels = list(c("","","","","", "A", "B", "C", "D", "", "E", "F", "G", "H"))) &
  theme(plot.tag.position = c(0.025, 0.9))

png("figures/supp_figs/svc_canopyHeight_lwr_link.png", height = 180, width = 200, res = 300, units = "mm")
p
dev.off()


# plot upper confidence interval of canopy height coef
p = 
  plot_spacer() + row1 + row2 + row3 + row4 +
  col1 + birds.sesvert.plot$upr_link + mammals.sesvert.plot$upr_link + rept.sesvert.plot$upr_link + amph.sesvert.plot$upr_link +
  col2 + birds.meanvert.plot$upr_link + mammals.meanvert.plot$upr_link + rept.meanvert.plot$upr_link + amph.meanvert.plot$upr_link + 
  plot_layout(design = design, heights = c(0.2,1,1,1,1), widths = c(0.2,1,1)) +
  plot_annotation(tag_levels = list(c("","","","","", "A", "B", "C", "D", "", "E", "F", "G", "H"))) &
  theme(plot.tag.position = c(0.025, 0.9))

png("figures/supp_figs/svc_canopyHeight_upr_link.png", height = 180, width = 200, res = 300, units = "mm")
p
dev.off()










