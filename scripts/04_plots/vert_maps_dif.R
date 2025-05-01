# Author: Lydia Soifer
# plot the difference of future verticality predictions - current verticality for each vertebrate class

library(tidyterra)
library(cowplot)
library(terra)
library(colorspace)
library(patchwork)

plt_rast_margin = function(r, v, plot.title, margin.title, base_text_size = 12, var, rd, scale_limits) {
  # r: spatraster for plotting
  # v: spatvector for outline
  # plot.title = plot title
  # margin.title = y-axis margin title
  
  # get minmax values from vector
  ymin = terra::ext(v)[3]
  if(ymin > -60) {ymin=-60}
  ymax = terra::ext(v)[4]
  xmax = terra::ext(v)[2]
  xmin = terra::ext(v)[1]
  
  # get color scale limits
  col.lims = round(minmax(r[[var]]),3)
  
  # create main plot
  plot_main = ggplot() +
    geom_spatvector(data = wd, color = alpha("black", 0.7), fill = "gray90", linewidth = .1) +
    geom_spatraster(data = r[[var]], maxcell = Inf) +
    # Overlay world
    geom_spatvector(data = wd, color = alpha("black", 0.7), fill = NA, linewidth = .1) +
    #geom_segment(aes(x = xmax, xend = xmax, y = -90, yend = ymax)) +
    geom_segment(aes(x = xmax, xend = xmax, y = ymin, yend = ymax)) +
    geom_segment(aes(x = xmin, xend = xmax, y = ymin, yend = ymin)) +
    #scale_fill_continuous_divergingx("spectral", na.value = NA, limits = col.lims) +
    scale_fill_continuous_divergingx("spectral", na.value = NA, limits = scale_limits) +
    coord_sf(ylim = c(ymin, ymax)) +
    theme_classic() +
    theme(legend.position = "none",
          legend.margin = margin(0,0,0,0),
          legend.key = element_blank(),
          legend.background = element_blank(),
          legend.key.width = unit(10, units = "mm"),
          legend.box.margin = margin(0,0,0,0),
          legend.box.spacing = unit(1,"mm"),
          legend.box.background = element_blank(),
          plot.background = element_blank(),
          panel.spacing = unit(0,"mm"),
          panel.background = element_blank(),
          legend.title = element_blank(),
          axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0,0,0,0), units = "mm"))
  
  # add titles on the secondary axis
  plot_main <- plot_main +
  # titles on secondary axis, for later
    scale_y_continuous(limits = c(ymin-60, ymax), breaks = c(-60,-30,0,30,60), 
                       sec.axis = dup_axis(name = margin.title), expand = expansion(mult = c(0.1,0))) +
    scale_x_continuous(expand = c(0,0)) +
    theme(axis.title.x = element_blank(),
          axis.title.y.left = element_blank(),
          axis.title.y.right = element_text(size = 12))
  
  # Getting averages by x,y
  marg_y <- r[[var]] %>%
    as_tibble(xy = TRUE) %>%
    drop_na() %>%
    group_by(y) %>%
    summarise(avg = mean(.data[[var]], na.rm = T))
  
  # Cowplot would delete axis, we make axis breaks
  br_4marginal <- c(round(min(marg_y$avg),rd), 0, round(max(marg_y$avg),rd))
  
  labs <- data.frame(labs = paste(prettyNum(br_4marginal, big.mark = " ")))
  labs$for_y <- ymin-10
  labs$y <- br_4marginal
  
  # y-axis margin
  plot_y <- axis_canvas(plot_main, axis = "y", coord_flip = T) +
    ggplot() +
    geom_col(data = marg_y, aes(y, avg, fill = avg), color = NA) +
    geom_line(data = marg_y, aes(y, avg), color = "black", show.legend = F) +
    geom_text(data = labs, aes(x = for_y, y = y, label = labs), size = 3) +
    #scale_fill_continuous_divergingx("spectral", limits = col.lims) +
    scale_fill_continuous_divergingx("spectral", limits = scale_limits) +
    #scale_x_continuous(limits = c(ymin-5, ymax), breaks = c(-60,-30,0,30,60)) +
    scale_x_continuous(limits = c(ymin-15, ymax), breaks = c(-60,-30,0,30,60), expand = expansion(mult = c(0.1,0))) +
    scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
    #coord_flip(xlim = c(-90,ymax), clip = "off") +
    coord_flip(xlim = c(ymin,ymax), clip = "off") +
    #geom_vline(xintercept = -90) +
    geom_vline(xintercept = ymin) +
    #geom_segment(aes(x = -90, xend = ymax, y = 0, yend=0), linetype = "dashed") +
    geom_segment(aes(x = ymin, xend = ymax, y = 0, yend=0), linetype = "dashed") +
    #theme_classic() +
    theme(axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(size = 12))
  
  # Combine all plots into one
  sizes_axis <- grid::unit(.3, "null")
  
  plot_final <- insert_yaxis_grob(plot_main, grob = get_panel(plot_y, panel = 2),
                                  position = "right",
                                  width = sizes_axis * 1.25
  )
  
  #gg_final <- ggdraw(plot_final) + draw_label(plot.title, x = 0.5, y = 1, vjust = 1.5)
  gg_final <- ggdraw(plot_final)
  return(gg_final)
}

wd = vect("data/original/rnaturalearth_world.shp") %>% 
  project("epsg:4326") %>% 
  crop(ext(-180,180,-60,83.64))


# 50 km resolution plots from sdmTMB models -------------------------------

# * amphibians SESVERT DIF -------------------------------------------------------------

amph = readRDS("results/sdmTMB_models2/predictions/amphibians_sesvert.rds")
amph.sesvert = amph$pred.f
pred = amph$pred

# load("results/sdmTMB_models/amphibians_sesvert.RData")
# amph.sesvert = pred.f
amph.sesvert$vertvar = "SES Mean Verticality"
# amph.sesvert$est.reldif = pred.f$est.dif/(max(pred$vert.mean.ses) - min(pred$vert.mean.ses))*100
#amph.sesvert$biome = pred$biome

amph.pdat = amph.sesvert %>% 
  mutate(x = x*1e5, y = y*1e5) %>% 
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326")

amph.sesvert.dif.plt = plt_rast_margin(r = amph.pdat, v = wd, margin.title = "Avg. Difference by Latitude",
                                       plot.title = "Amphibians", var = "est.dif", rd = 1, 
                                       scale_limits = c(-1.5,1.5))
amph.sesvert.dif.plt



# * amphibians Meanvert DIF -------------------------------------------------------------

load("results/sdmTMB_models/amphibians_meanvert.RData")
# amph.meanvert = pred.f
# amph.meanvert = pred.f
# amph.meanvert$vertvar = "Mean Verticality"
# amph.meanvert$est.reldif = pred.f$est.dif/(max(pred$vert.mean) - min(pred$vert.mean))*100
# 
# 
# amph.pdat = amph.meanvert %>% 
#   mutate(x = x*1e5, y = y*1e5) %>% 
#   rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
#   project("epsg:4326")
# 
# amph.meanvert.dif.plt = plt_rast_margin(r = amph.pdat, v = wd, margin.title = "avg. relative\ndifference by latitude",
#                                        plot.title = "Amphibians", var = "est.reldif", rd = 2, 
#                                        scale_limits = c(-5.2,5.2))
# amph.meanvert.dif.plt

# * amphibians Parb DIF -------------------------------------------------------------

# load("results/sdmTMB_models/amphibians_parb.RData")
# amph.parb = pred.f
# amph.parb = pred.f
# amph.parb$vertvar = "Proportion arboreal"
# amph.parb$est.reldif = pred.f$est.dif/(max(pred$p.arb) - min(pred$p.arb))*100
# 
# 
# amph.pdat = amph.parb %>% 
#   mutate(x = x*1e5, y = y*1e5) %>% 
#   rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
#   project("epsg:4326")
# 
# amph.parb.dif.plt = plt_rast_margin(r = amph.pdat, v = wd, margin.title = "avg. relative\ndifference by latitude",
#                                         plot.title = "Amphibians", var = "est.reldif", rd = 2, 
#                                         scale_limits = c(-7.5,7.5))
# amph.parb.dif.plt


# * reptiles SESVERT DIF -------------------------------------------------------------

repts = readRDS("results/sdmTMB_models2/predictions/reptiles_sesvert.rds")
reptiles.sesvert = repts$pred.f
# load("results/sdmTMB_models/reptiles_sesvert.RData")
# reptiles.sesvert = pred.f
reptiles.sesvert$vertvar = "SES Mean Verticality"
# reptiles.sesvert$est.reldif = pred.f$est.dif/(max(pred$vert.mean.ses) - min(pred$vert.mean.ses))*100

rept.pdat = reptiles.sesvert %>% 
  mutate(x = x*1e5, y = y*1e5) %>% 
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326")

rept.sesvert.dif.plt = plt_rast_margin(r = rept.pdat, v = wd, margin.title = "Avg. Difference by Latitude",
                                       plot.title = "Reptiles", var = "est.dif", rd = 1, 
                                       scale_limits = c(-1.5,1.5))
rept.sesvert.dif.plt


# * Reptiles MEANVERT DIF -------------------------------------------------------------

# load("results/sdmTMB_models/reptiles_meanvert.RData")
# reptiles.meanvert = pred.f
# reptiles.meanvert$vertvar = "Mean Verticality"
# reptiles.meanvert$est.reldif = pred.f$est.dif/(max(pred$vert.mean) - min(pred$vert.mean))*100
# 
# rept.pdat = reptiles.meanvert %>% 
#   mutate(x = x*1e5, y = y*1e5) %>% 
#   rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
#   project("epsg:4326")
# 
# rept.meanvert.dif.plt = plt_rast_margin(r = rept.pdat, v = wd, margin.title = "avg. relative\ndifference by latitude",
#                                         plot.title = "Reptiles", var = "est.reldif", rd = 3, 
#                                         scale_limits = c(-5.2,5.2))
# rept.meanvert.dif.plt


# * reptiles Parb DIF -------------------------------------------------------------

# load("results/sdmTMB_models/reptiles_parb.RData")
# reptiles.parb = pred.f
# reptiles.parb = pred.f
# reptiles.parb$vertvar = "Proportion arboreal"
# reptiles.parb$est.reldif = pred.f$est.dif/(max(pred$p.arb) - min(pred$p.arb))*100
# 
# 
# reptiles.pdat = reptiles.parb %>% 
#   mutate(x = x*1e5, y = y*1e5) %>% 
#   rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
#   project("epsg:4326")
# 
# reptiles.parb.dif.plt = plt_rast_margin(r = reptiles.pdat, v = wd, margin.title = "avg. relative\ndifference by latitude",
#                                        plot.title = "Reptiles", var = "est.reldif", rd = 2, 
#                                        scale_limits = c(-7.5,7.5))
# reptiles.parb.dif.plt



# * mammals SESVERT DIF -------------------------------------------------------------

mammals = readRDS("results/sdmTMB_models2/predictions/mammals_sesvert.rds")
mammals.sesvert = mammals$pred.f
# load("results/sdmTMB_models/mammals_sesvert.RData")
# mammals.sesvert = pred.f
mammals.sesvert$vertvar = "SES Mean Verticality"
# mammals.sesvert$est.reldif = pred.f$est.dif/(max(pred$vert.mean.ses) - min(pred$vert.mean.ses))*100


mammals.pdat = mammals.sesvert %>% 
  mutate(x = x*1e5, y = y*1e5) %>% 
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326")

mammals.sesvert.dif.plt = plt_rast_margin(r = mammals.pdat, v = wd, margin.title = "Avg. Difference by Latitude",
                                       plot.title = "Mammals", var = "est.dif", rd = 1, 
                                       scale_limits = c(-1.5,1.5))
mammals.sesvert.dif.plt



# * mammals MEANVERT DIF -------------------------------------------------------------

# load("results/sdmTMB_models/mammals_meanvert.RData")
# mammals.meanvert = pred.f
# mammals.meanvert$vertvar = "Mean Verticality"
# mammals.meanvert$est.reldif = pred.f$est.dif/(max(pred$vert.mean) - min(pred$vert.mean))*100
# 
# 
# mammals.pdat = mammals.meanvert %>% 
#   mutate(x = x*1e5, y = y*1e5) %>% 
#   rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
#   project("epsg:4326")
# 
# mammals.meanvert.dif.plt = plt_rast_margin(r = mammals.pdat, v = wd, margin.title = "avg. relative\ndifference by latitude",
#                                         plot.title = "Mammals", var = "est.reldif", rd = 3,
#                                         scale_limits = c(-5.2,5.2))
# mammals.meanvert.dif.plt

# * mammals Parb DIF -------------------------------------------------------------

# load("results/sdmTMB_models/mammals_parb.RData")
# mammals.parb = pred.f
# mammals.parb = pred.f
# mammals.parb$vertvar = "Proportion arboreal"
# mammals.parb$est.reldif = pred.f$est.dif/(max(pred$p.arb) - min(pred$p.arb))*100
# 
# 
# mammals.pdat = mammals.parb %>% 
#   mutate(x = x*1e5, y = y*1e5) %>% 
#   rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
#   project("epsg:4326")
# 
# mammals.parb.dif.plt = plt_rast_margin(r = mammals.pdat, v = wd, margin.title = "avg. relative\ndifference by latitude",
#                                      plot.title = "Mammals", var = "est.reldif", rd = 2, 
#                                      scale_limits = c(-7.5,7.5))
# mammals.parb.dif.plt


# * birds (breeding + resident) SESVERT DIF -------------------------------------------------------------

birds = readRDS("results/sdmTMB_models2/predictions/birds_sesvert.rds")
birds.sesvert = birds$pred.f
# load("results/sdmTMB_models/birds_sesvert.RData")
# birds.sesvert = pred.f
birds.sesvert$vertvar = "SES Mean Verticality"
# birds.sesvert$est.reldif = pred.f$est.dif/(max(pred$vert.mean.ses) - min(pred$vert.mean.ses))*100


birds.pdat = birds.sesvert %>% 
  mutate(x = x*1e5, y = y*1e5) %>% 
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326")

birds.sesvert.dif.plt = plt_rast_margin(r = birds.pdat, v = wd, margin.title = "Avg. Difference by Latitude",
                                       plot.title = "Birds", var = "est.dif", rd = 1,
                                       scale_limits = c(-1.5,1.5))
birds.sesvert.dif.plt




# * Birds MEANVERT DIF -------------------------------------------------------------


# load("results/sdmTMB_models/birds_meanvert.RData")
# birds.meanvert = pred.f
# birds.meanvert$vertvar = "Mean Verticality"
# birds.meanvert$est.reldif = pred.f$est.dif/(max(pred$vert.mean) - min(pred$vert.mean))*100
# 
# 
# birds.pdat = birds.meanvert %>% 
#   mutate(x = x*1e5, y = y*1e5) %>% 
#   rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
#   project("epsg:4326")
# 
# birds.meanvert.dif.plt = plt_rast_margin(r = birds.pdat, v = wd, margin.title = "avg. relative\ndifference by latitude",
#                                         plot.title = "Birds", var = "est.reldif", rd = 3, 
#                                         scale_limits = c(-5.2,5.2))
# birds.meanvert.dif.plt
# 
# 
# # * birds Parb DIF -------------------------------------------------------------
# 
# load("results/sdmTMB_models/birds_parb.RData")
# birds.parb = pred.f
# birds.parb = pred.f
# birds.parb$vertvar = "Proportion arboreal"
# birds.parb$est.reldif = pred.f$est.dif/(max(pred$p.arb) - min(pred$p.arb))*100
# 
# 
# birds.pdat = birds.parb %>% 
#   mutate(x = x*1e5, y = y*1e5) %>% 
#   rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
#   project("epsg:4326")
# 
# birds.parb.dif.plt = plt_rast_margin(r = birds.pdat, v = wd, margin.title = "avg. relative\ndifference by latitude",
#                                     plot.title = "Birds", var = "est.reldif", rd = 2, 
#                                     scale_limits = c(-7.5,7.5))
# birds.parb.dif.plt
# 



# FULL FIGURE -------------------------------------------------------------

design = "
AB
CD"

library(grid)

# make combined legend
legend.sesvert = ggplot(data = data.frame(x = 1, y = 1, z = c(-1.5,1.5))) + 
  geom_point(aes(x,y,fill = z)) + scale_fill_continuous_divergingx("spectral", limits = c(-1.5,1.5), guide = guide_colorbar("")) +
  theme(legend.direction = "horizontal",
        legend.key.width = unit(10, units = "mm"))
legend.sesvert = get_legend(legend.sesvert)


p = birds.sesvert.dif.plt + mammals.sesvert.dif.plt + 
  rept.sesvert.dif.plt + amph.sesvert.dif.plt

svg("figures/main_figs/fig5/fig5_vert_difs.svg", height = 7.5, width = 14)
p
dev.off()

svg("figures/main_figs/fig5/fig5_legend.svg", width = 7, height = 3)
grid.draw(legend.sesvert)
dev.off()


# FULL FIGURE -------------------------------------------------------------

# design = "
# AGM
# BHN
# CIO
# DJP
# EKQ
# FLR"
# 
# library(grid)
# col1 = wrap_elements(panel = textGrob("SES Mean Verticality", gp= gpar(fontface = "bold")))
# col2 = wrap_elements(panel = textGrob("Mean Verticality", gp = gpar(fontface = "bold")))
# 
# row0 = wrap_elements(panel = textGrob("", rot = 90))
# row1 = wrap_elements(panel = textGrob("Birds", rot = 90, gp = gpar(fontface = "bold")))
# row2 = wrap_elements(panel = textGrob("Mammals", rot = 90, gp = gpar(fontface = "bold")))
# row3 = wrap_elements(panel = textGrob("Reptiles", rot = 90, gp = gpar(fontface = "bold")))
# row4 = wrap_elements(panel = textGrob("Amphibians", rot = 90, gp = gpar(fontface = "bold")))
# row5 = wrap_elements(panel = textGrob("", rot = 90))
# 
# # make combined legend
# legend.sesvert = ggplot(data = data.frame(x = 1, y = 1, z = c(-5.2,5.2))) + 
#   geom_point(aes(x,y,fill = z)) + scale_fill_continuous_divergingx("spectral", limits = c(-5.2,5.2), guide = guide_colorbar("")) +
#   theme(legend.direction = "horizontal",
#         legend.key.width = unit(10, units = "mm"))
# legend.sesvert = get_legend(legend.sesvert)
# 
# legend.meanvert = ggplot(data = data.frame(x = 1, y = 1, z = c(-5.2,5.2))) + 
#   geom_point(aes(x,y,fill = z)) + scale_fill_continuous_divergingx("spectral", limits = c(-5.2,5.2), guide = guide_colorbar("")) +
#   theme(legend.direction = "horizontal",
#         legend.key.width = unit(10, units = "mm"))
# legend.meanvert = get_legend(legend.meanvert)
# 
# 
# p = row0 + row1 + row2 + row3 + row4 + row5 +
#   col1 + birds.sesvert.dif.plt + mammals.sesvert.dif.plt + rept.sesvert.dif.plt + amph.sesvert.dif.plt + legend.sesvert +
#   col2 + birds.meanvert.dif.plt + mammals.meanvert.dif.plt + rept.meanvert.dif.plt + amph.meanvert.dif.plt + legend.meanvert +
#   plot_layout(design = design, heights = c(0.1, -1, -1, -1, -1, 0.2), widths = c(0.1, -1, -1, 0.2)) +
#   plot_annotation(tag_levels = list(c("", "", "", "", "", "", "", "a", "b", "c", "d", "", "", "e", "f", "g", "h","")),
#                   theme = theme(plot.margin = margin(0,0,0,0, unit = "mm"),
#                                 panel.spacing = margin(0,0,0,0, unit = "mm"))) &
#   theme(plot.title = element_blank(),
#         panel.spacing = margin(0,0,0,0, unit = "mm"),
#         panel.background = element_blank(),
#         plot.tag.position = c(0.025, 0.95))
# 
# 
# png("figures/main_figs/fig5_vert_reldifs.png", height = 215, width = 300, res = 300, units = "mm")
# p
# dev.off()
# 
# png("figures/main_figs/fig5_vert_reldifs_test.png", height = 150, width = 180, res = 600, units = "mm")
# p
# dev.off()
# 
# ggsave("figures/main_figs/fig5_vert_reldifs_test.svg", height = 150, width = 180, dpi = 600, units = "mm")
# 
# # PROPORTION ARBOREAL FIGURE -------------------------------------------------------------
# 
# design = "
# AG
# BH
# CI
# DJ
# EK
# FL"
# 
# library(grid)
# col1 = wrap_elements(panel = textGrob("Proportion Arboreal", gp= gpar(fontface = "bold")))
# 
# row0 = wrap_elements(panel = textGrob("", rot = 90))
# row1 = wrap_elements(panel = textGrob("Birds", rot = 90, gp = gpar(fontface = "bold")))
# row2 = wrap_elements(panel = textGrob("Mammals", rot = 90, gp = gpar(fontface = "bold")))
# row3 = wrap_elements(panel = textGrob("Reptiles", rot = 90, gp = gpar(fontface = "bold")))
# row4 = wrap_elements(panel = textGrob("Amphibians", rot = 90, gp = gpar(fontface = "bold")))
# row5 = wrap_elements(panel = textGrob("", rot = 90))
# 
# # make combined legend
# legend.parb = ggplot(data = data.frame(x = 1, y = 1, z = c(-7.5,7.5))) + 
#   geom_point(aes(x,y,fill = z)) + scale_fill_continuous_divergingx("spectral", limits = c(-7.5,7.5), guide = guide_colorbar("")) +
#   theme(legend.direction = "horizontal",
#         legend.key.width = unit(10, units = "mm"))
# legend.parb = get_legend(legend.parb)
# 
# 
# 
# p = row0 + row1 + row2 + row3 + row4 + row5 +
#   col1 + birds.parb.dif.plt + mammals.parb.dif.plt + reptiles.parb.dif.plt + amph.parb.dif.plt + legend.parb +
#   plot_layout(design = design, heights = c(0.1, -1, -1, -1, -1, 0.2), widths = c(0.1, -1, 0.2)) +
#   plot_annotation(tag_levels = list(c("", "", "", "", "", "", "", "A", "B", "C", "D", ""))) +
#   theme(plot.title = element_blank())
# 
# 
# png("figures/main_figs/vert_parb_reldifs.png", height = 250, width = 150, res = 300, units = "mm")
# p
# dev.off()

