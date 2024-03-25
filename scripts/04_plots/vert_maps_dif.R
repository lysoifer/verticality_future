# Author: Lydia Soifer
# plot the difference of future verticality predictions - current verticality for each vertebrate class

library(tidyterra)
library(cowplot)
plt_rast_margin = function(r, v, plot.title, margin.title, base_text_size = 12, var) {
  # r: spatraster for plotting
  # v: spatvector for outline
  # plot.title = plot title
  # margin.titl = y-axis margin title
  
  # get minmax values from vector
  ymin = ext(v)[3]
  ymax = ext(v)[4]
  xmax = ext(v)[2]
  
  # get color scale limits
  col.lims = round(minmax(r),1)
  
  # create main plot
  plot_main = ggplot() +
    geom_spatraster(data = r, maxcell = Inf) +
    # Overlay world
    geom_spatvector(data = wd, color = alpha("black", 0.7), fill = NA, linewidth = .1) +
    geom_segment(aes(x = xmax, xend = xmax, y = -90, yend = ymax)) +
    scale_fill_continuous_divergingx("spectral", na.value = NA, limits = col.lims) +
    coord_sf(ylim = c(ymin, ymax)) +
    theme_classic() +
    theme(legend.position = "bottom",
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
          legend.title = element_blank())
  
  # add titles on the secondary axis
  plot_main <- plot_main +
  # titles on secondary axis, for later
    scale_y_continuous(limits = c(ymin-5, ymax), breaks = c(-60,-30,0,30,60), 
                       sec.axis = dup_axis(name = margin.title)) +
    scale_x_continuous(expand = c(0,0)) +
    theme(axis.title.x = element_blank(),
          axis.title.y.left = element_blank())
  
  # Getting averages by x,y
  marg_y <- r %>%
    as_tibble(xy = TRUE) %>%
    drop_na() %>%
    group_by(y) %>%
    summarise(avg = mean(.data[[var]], na.rm = T))
  
  # Cowplot would delete axis, we make axis breaks
  br_4marginal <- c(round(min(marg_y$avg),1), 0, round(max(marg_y$avg),1))
  
  labs <- data.frame(labs = paste(prettyNum(br_4marginal, big.mark = " ")))
  labs$for_y <- ymin-5
  labs$y <- br_4marginal
  
  # y-axis margin
  plot_y <- axis_canvas(plot_main, axis = "y", coord_flip = TRUE) +
    ggplot() +
    geom_col(data = marg_y, aes(y, avg, fill = avg), color = NA) +
    geom_line(data = marg_y, aes(y, avg), color = "black", show.legend = F) +
    geom_text(data = labs, aes(x = for_y, y = y, label = labs), size = 3) +
    scale_fill_continuous_divergingx("spectral", limits = col.lims) +
    scale_x_continuous(limits = c(ymin-5, ymax), breaks = c(-60,-30,0,30,60)) +
    scale_y_continuous(expand = c(0.1,0.1)) +
    coord_flip(xlim = c(-90,ymax), clip = "off") +
    geom_vline(xintercept = -90) +
    geom_segment(aes(x = -90, xend = ymax, y = 0, yend=0), linetype = "dashed") +
    theme_classic() +
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
  
  gg_final <- ggdraw(plot_final) + draw_label(plot.title, x = 0.5, y = 1, vjust = 1.5)
  return(gg_final)
}

wd = vect("data/original/rnaturalearth_world.shp")
wd = project(wd, "epsg:4326")

load("results/amphibians_sar_sesvert2.RData")
amph = sesvert.scale.sub

amph.pdat = amph %>% 
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326") %>% 
  as.data.frame(xy = T) %>% 
  mutate(sesvert.dif = vert.mean.ses.future-vert.mean.ses) %>% 
  dplyr::select(x,y,sesvert.dif) %>% 
  rast(type = "xyz", crs = "epsg:4326")

amph.sesvert.dif.plt = plt_rast_margin(amph.pdat, v = wd, margin.title = "avg. difference by latitude", plot.title = "Amphibians")
amph.sesvert.dif.plt

png("figs_lydia/maps/sesvert_future_dif/amph.png", width = 200, height = 100, res = 300, units = "mm")
amph.sesvert.dif.plt
dev.off()


load("results/mammals_sar_sesvert2.RData")
mammals = sesvert.scale.sub

mammals.pdat = mammals %>% 
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326") %>% 
  as.data.frame(xy = T) %>% 
  mutate(sesvert.dif = vert.mean.ses.future-vert.mean.ses) %>% 
  dplyr::select(x,y,sesvert.dif) %>% 
  rast(type = "xyz", crs = "epsg:4326")

mammals.sesvert.dif.plt = plt_rast_margin(mammals.pdat, v = wd, margin.title = "avg. difference by latitude", plot.title = "Mammals")
mammals.sesvert.dif.plt

png("figs_lydia/maps/sesvert_future_dif/mammals.png", width = 200, height = 100, res = 300, units = "mm")
mammals.sesvert.dif.plt
dev.off()


load("results/birdselton_sar_sesvert2.RData")
birds = sesvert.scale.sub

birds.pdat = birds %>% 
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326") %>% 
  as.data.frame(xy = T) %>% 
  mutate(sesvert.dif = vert.mean.ses.future-vert.mean.ses) %>% 
  dplyr::select(x,y,sesvert.dif) %>% 
  rast(type = "xyz", crs = "epsg:4326")

birds.sesvert.dif.plt = plt_rast_margin(birds.pdat, v = wd, margin.title = "avg. difference by latitude", plot.title = "Birds")
birds.sesvert.dif.plt

png("figs_lydia/maps/sesvert_future_dif/birds.png", width = 200, height = 100, res = 300, units = "mm")
birds.sesvert.dif.plt
dev.off()


load("results/reptiles_sar_sesvert2.RData")
repts = sesvert.scale.sub

repts.pdat = repts %>% 
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326") %>% 
  as.data.frame(xy = T) %>% 
  mutate(sesvert.dif = vert.mean.ses.future-vert.mean.ses) %>% 
  dplyr::select(x,y,sesvert.dif) %>% 
  rast(type = "xyz", crs = "epsg:4326")

repts.sesvert.dif.plt = plt_rast_margin(repts.pdat, v = wd, margin.title = "avg. difference by latitude", plot.title = "Reptiles")
repts.sesvert.dif.plt

png("figs_lydia/maps/sesvert_future_dif/repts.png", width = 200, height = 100, res = 300, units = "mm")
repts.sesvert.dif.plt
dev.off()


# MEAN VERT DIF -----------------------------------------------------------

wd = vect("data/original/rnaturalearth_world.shp")
wd = project(wd, "epsg:4326")

load("results/amphibians_sar_vertmean2.RData")
amph = vert.scale.sub

amph.pdat = amph %>% 
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326") %>% 
  as.data.frame(xy = T) %>% 
  mutate(vert.dif = vert.mean.future-vert.mean) %>% 
  dplyr::select(x,y,vert.dif) %>% 
  rast(type = "xyz", crs = "epsg:4326")

amph.vert.dif.plt = plt_rast_margin(r = amph.pdat, v = wd, margin.title = "avg. difference by latitude",
                                    plot.title = "Amphibians", var = "vert.dif")
amph.vert.dif.plt

png("figs_lydia/maps/vert_future_dif/amph.png", width = 200, height = 100, res = 300, units = "mm")
amph.vert.dif.plt
dev.off()


load("results/mammals_sar_vertmean2.RData")
mammals = vert.scale.sub

mammals.pdat = mammals %>% 
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326") %>% 
  as.data.frame(xy = T) %>% 
  mutate(vert.dif = vert.mean.future-vert.mean) %>% 
  dplyr::select(x,y,vert.dif) %>% 
  rast(type = "xyz", crs = "epsg:4326")

mammals.vert.dif.plt = plt_rast_margin(mammals.pdat, v = wd, margin.title = "avg. difference by latitude", 
                                       plot.title = "Mammals", var = "vert.dif")
mammals.vert.dif.plt

png("figs_lydia/maps/vert_future_dif/mammals.png", width = 200, height = 100, res = 300, units = "mm")
mammals.vert.dif.plt
dev.off()


load("results/birdselton_sar_vertmean.RData")
birds = vert.scale.sub

birds.pdat = birds %>% 
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326") %>% 
  as.data.frame(xy = T) %>% 
  mutate(vert.dif = vert.mean.future-vert.mean) %>% 
  dplyr::select(x,y,vert.dif) %>% 
  rast(type = "xyz", crs = "epsg:4326")

birds.vert.dif.plt = plt_rast_margin(birds.pdat, v = wd, margin.title = "avg. difference by latitude", plot.title = "Birds", var = "vert.dif")
birds.vert.dif.plt

png("figs_lydia/maps/vert_future_dif/birds.png", width = 200, height = 100, res = 300, units = "mm")
birds.vert.dif.plt
dev.off()


load("results/reptiles_sar_vertmean2.RData")
repts = vert.scale.sub

repts.pdat = repts %>% 
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326") %>% 
  as.data.frame(xy = T) %>% 
  mutate(vert.dif = vert.mean.future-vert.mean) %>% 
  dplyr::select(x,y,vert.dif) %>% 
  rast(type = "xyz", crs = "epsg:4326")

repts.vert.dif.plt = plt_rast_margin(repts.pdat, v = wd, margin.title = "avg. difference by latitude", plot.title = "Reptiles", var = "vert.dif")
repts.vert.dif.plt

png("figs_lydia/maps/vert_future_dif/repts.png", width = 200, height = 100, res = 300, units = "mm")
repts.vert.dif.plt
dev.off()


















