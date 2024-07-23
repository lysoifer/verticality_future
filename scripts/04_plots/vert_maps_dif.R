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
          axis.ticks.x = element_blank())
  
  # add titles on the secondary axis
  plot_main <- plot_main +
  # titles on secondary axis, for later
    scale_y_continuous(limits = c(ymin-5, ymax), breaks = c(-60,-30,0,30,60), 
                       sec.axis = dup_axis(name = margin.title), expand = expansion(mult = c(0.1,0))) +
    scale_x_continuous(expand = c(0,0)) +
    theme(axis.title.x = element_blank(),
          axis.title.y.left = element_blank())
  
  # Getting averages by x,y
  marg_y <- r[[var]] %>%
    as_tibble(xy = TRUE) %>%
    drop_na() %>%
    group_by(y) %>%
    summarise(avg = mean(.data[[var]], na.rm = T))
  
  # Cowplot would delete axis, we make axis breaks
  br_4marginal <- c(round(min(marg_y$avg),rd), 0, round(max(marg_y$avg),rd))
  
  labs <- data.frame(labs = paste(prettyNum(br_4marginal, big.mark = " ")))
  labs$for_y <- ymin-5
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
          panel.background = element_blank())
  
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

load("results/sdmTMB_models/amphibians_sesvert.RData")
amph = pred.f

amph.pdat = amph %>% 
  mutate(x = x*1e5, y = y*1e5) %>% 
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326")

amph.sesvert.dif.plt = plt_rast_margin(r = amph.pdat, v = wd, margin.title = "avg. difference by latitude",
                                       plot.title = "Amphibians", var = "est.dif", rd = 1, scale_limits = c(-2,2))
amph.sesvert.dif.plt



# * amphibians Meanvert DIF -------------------------------------------------------------

load("results/sdmTMB_models/amphibians_meanvert.RData")
amph = pred.f

amph.pdat = amph %>% 
  mutate(x = x*1e5, y = y*1e5) %>% 
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326")

amph.meanvert.dif.plt = plt_rast_margin(r = amph.pdat, v = wd, margin.title = "avg. difference by latitude",
                                       plot.title = "Amphibians", var = "est.dif", rd = 2, scale_limits = c(-0.5,0.5))
amph.meanvert.dif.plt



# * reptiles SESVERT DIF -------------------------------------------------------------

load("results/sdmTMB_models/reptiles_sesvert.RData")
rept= pred.f

rept.pdat = rept %>% 
  mutate(x = x*1e5, y = y*1e5) %>% 
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326")

rept.sesvert.dif.plt = plt_rast_margin(r = rept.pdat, v = wd, margin.title = "avg. difference by latitude",
                                       plot.title = "Reptiles", var = "est.dif", rd = 1, scale_limits = c(-2,2))
rept.sesvert.dif.plt


# * Reptiles MEANVERT DIF -------------------------------------------------------------

load("results/sdmTMB_models/reptiles_meanvert.RData")
rept = pred.f

rept.pdat = rept %>% 
  mutate(x = x*1e5, y = y*1e5) %>% 
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326")

rept.meanvert.dif.plt = plt_rast_margin(r = rept.pdat, v = wd, margin.title = "avg. difference by latitude",
                                        plot.title = "Reptiles", var = "est.dif", rd = 3, scale_limits = c(-0.5,0.5))
rept.meanvert.dif.plt



# * mammals SESVERT DIF -------------------------------------------------------------

load("results/sdmTMB_models/mammals_sesvert.RData")
mammals = pred.f

mammals.pdat = mammals %>% 
  mutate(x = x*1e5, y = y*1e5) %>% 
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326")

mammals.sesvert.dif.plt = plt_rast_margin(r = mammals.pdat, v = wd, margin.title = "avg. difference by latitude",
                                       plot.title = "Mammals", var = "est.dif", rd = 1, scale_limits = c(-2,2))
mammals.sesvert.dif.plt



# * mammals MEANVERT DIF -------------------------------------------------------------

load("results/sdmTMB_models/mammals_meanvert.RData")
mammals = pred.f

mammals.pdat = mammals %>% 
  mutate(x = x*1e5, y = y*1e5) %>% 
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326")

mammals.meanvert.dif.plt = plt_rast_margin(r = mammals.pdat, v = wd, margin.title = "avg. difference by latitude",
                                        plot.title = "Mammals", var = "est.dif", rd = 3, scale_limits = c(-0.5,0.5))
mammals.meanvert.dif.plt



# * birds (breeding + resident) SESVERT DIF -------------------------------------------------------------

load("results/sdmTMB_models/birds_sesvert.RData")
birds = pred.f

birds.pdat = birds %>% 
  mutate(x = x*1e5, y = y*1e5) %>% 
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326")

birds.sesvert.dif.plt = plt_rast_margin(r = birds.pdat, v = wd, margin.title = "avg. difference by latitude",
                                       plot.title = "Birds", var = "est.dif", rd = 1, scale_limits = c(-2,2))
birds.sesvert.dif.plt




# * Birds MEANVERT DIF -------------------------------------------------------------


load("results/sdmTMB_models/birds_meanvert.RData")
birds = pred.f

birds.pdat = birds %>% 
  mutate(x = x*1e5, y = y*1e5) %>% 
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326")

birds.meanvert.dif.plt = plt_rast_margin(r = birds.pdat, v = wd, margin.title = "avg. difference by latitude",
                                        plot.title = "Birds", var = "est.dif", rd = 3, scale_limits = c(-0.5, 0.5))
birds.meanvert.dif.plt




# FULL FIGURE -------------------------------------------------------------

design = "
AGM
BHN
CIO
DJP
EKQ
FLR"

library(grid)
col1 = wrap_elements(panel = textGrob("SES Mean Verticality", gp= gpar(fontface = "bold")))
col2 = wrap_elements(panel = textGrob("Mean Verticality", gp = gpar(fontface = "bold")))

row0 = wrap_elements(panel = textGrob("", rot = 90))
row1 = wrap_elements(panel = textGrob("Birds", rot = 90, gp = gpar(fontface = "bold")))
row2 = wrap_elements(panel = textGrob("Mammals", rot = 90, gp = gpar(fontface = "bold")))
row3 = wrap_elements(panel = textGrob("Reptiles", rot = 90, gp = gpar(fontface = "bold")))
row4 = wrap_elements(panel = textGrob("Amphibians", rot = 90, gp = gpar(fontface = "bold")))
row5 = wrap_elements(panel = textGrob("", rot = 90))

# make combined legend
legend.sesvert = ggplot(data = data.frame(x = 1, y = 1, z = -2,2)) + 
  geom_point(aes(x,y,fill = z)) + scale_fill_continuous_divergingx("spectral", limits = c(-2,2), guide = guide_colorbar("")) +
  theme(legend.direction = "horizontal",
        legend.key.width = unit(10, units = "mm"))
legend.sesvert = get_legend(legend.sesvert)

legend.meanvert = ggplot(data = data.frame(x = 1, y = 1, z = -0.5,0.5)) + 
  geom_point(aes(x,y,fill = z)) + scale_fill_continuous_divergingx("spectral", limits = c(-0.5,0.5), guide = guide_colorbar("")) +
  theme(legend.direction = "horizontal",
        legend.key.width = unit(10, units = "mm"))
legend.meanvert = get_legend(legend.meanvert)


p = row0 + row1 + row2 + row3 + row4 + row5 +
  col1 + birds.sesvert.dif.plt + mammals.sesvert.dif.plt + rept.sesvert.dif.plt + amph.sesvert.dif.plt + legend.sesvert +
  col2 + birds.meanvert.dif.plt + mammals.meanvert.dif.plt + rept.meanvert.dif.plt + amph.meanvert.dif.plt + legend.meanvert +
  plot_layout(design = design, heights = c(0.1, -1, -1, -1, -1, 0.2), widths = c(0.1, -1, -1, 0.2)) +
  plot_annotation(tag_levels = list(c("", "", "", "", "", "", "", "A", "B", "C", "D", "", "", "E", "F", "G", "H",""))) +
  theme(plot.title = element_blank())


png("figures/main_figs/vert_difs.png", height = 250, width = 300, res = 300, units = "mm")
p
dev.off()
