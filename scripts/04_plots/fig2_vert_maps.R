library(tidyverse)
library(png)
library(ggplot2)
library(patchwork)
library(terra)
library(data.table)
library(tidyterra)
library(colorspace)
library(grid)
library(tricolore)
library(lemon)


# add proportion plots ----------------------------------------------------

wd = vect("data/original/rnaturalearth_world.shp") %>% 
  project("epsg:4326") %>% 
  crop(ext(-180,180,-60,83.64))
ymin = -60
ymax = 83.64

amph = fread("data/derivative_data/gridcell_data/env_forest/50_km/amph_comdat.csv") %>% 
  dplyr::select(x, y, rich, p.arb, p.fos, p.ter, vert.mean, vert.mean.ses, biome) %>% 
  filter(rich >= 5) %>% 
  filter(!is.na(vert.mean.ses)) %>% 
  drop_na() %>% 
  rast(crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326") %>% 
  as.data.frame(xy = T) %>% 
  mutate(taxa = "Amphibians")

mammals = fread("data/derivative_data/gridcell_data/env_forest/50_km/mammals_comdat.csv") %>% 
  dplyr::select(x, y, rich, p.arb, p.fos, p.ter, vert.mean, vert.mean.ses, biome) %>% 
  filter(rich >= 5) %>% 
  filter(!is.na(vert.mean.ses)) %>% 
  drop_na() %>% 
  rast(crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326") %>% 
  as.data.frame(xy = T) %>% 
  mutate(taxa = "Mammals")

birds = fread("data/derivative_data/gridcell_data/env_forest/50_km/birds_comdat.csv") %>% 
  dplyr::select(x, y, rich, p.arb, p.fos, p.ter, vert.mean, vert.mean.ses, biome) %>% 
  filter(rich >= 5) %>% 
  filter(!is.na(vert.mean.ses)) %>% 
  drop_na()%>% 
  rast(crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326") %>% 
  as.data.frame(xy = T) %>%  
  mutate(taxa = "Birds")

repts = fread("data/derivative_data/gridcell_data/env_forest/50_km/reptiles_comdat.csv") %>% 
  dplyr::select(x, y, rich, p.arb, p.fos, p.ter, vert.mean, vert.mean.ses, biome) %>% 
  filter(rich >= 5) %>% 
  filter(!is.na(vert.mean.ses)) %>% 
  drop_na()%>% 
  rast(crs = "+proj=cea +datum=WGS84") %>% 
  project("epsg:4326") %>% 
  as.data.frame(xy = T) %>%  
  mutate(taxa = "Reptiles")

all = bind_rows(amph, mammals, birds, repts) %>% 
  drop_na(vert.mean.ses) %>% 
  mutate(taxa = factor(taxa, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")))
min = min(all$vert.mean.ses, na.rm = T)
max = max(all$vert.mean.ses, na.rm = T)
col.lims = range(all$vert.mean.ses, na.rm = T)

rgbcols = Tricolore(all,  "p.ter", "p.fos", "p.arb" ,breaks = Inf, spread = 1,
                    show_data = F, hue = 0.2, chroma = 1, lightness = 0.8)
all$rgb = rgbcols$rgb

thm = theme(legend.key.height = unit(3, units = "mm"),
            legend.key.width = unit(10, units = "mm"),
            plot.background = element_blank(),
            plot.margin = unit(c(0,0,0,0), units = "mm"),
            panel.spacing = unit(0,"mm"),
            panel.background = element_rect(fill = "transparent", color = NA),
            strip.background = element_rect(fill = "white", color = "black"),
            panel.grid = element_blank(),
            legend.title = element_text(hjust = 0.5, size = 10),
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.text = element_blank())

library(ggpubr)
p = ggplot(all) +
  geom_sf(data = wd) +
  geom_raster(data = all, aes(x,y,fill = rgb)) +
  geom_sf(data = wd, fill = NA, color = "black") +
  scale_fill_identity() +
  facet_wrap(vars(taxa), ncol = 1, strip.position = "left") +
  coord_sf(crs = "epsg:4326") +
  scale_x_continuous("Vertical Strategies", expand = c(0,0)) +
  thm + 
  theme(strip.background = element_blank(),
        strip.text = element_blank())


p2 = ggplot(all) +
  geom_spatvector(data = wd) +
  geom_raster(aes(x,y,fill = vert.mean.ses)) +
  geom_spatvector(data = wd, fill = NA, color = "black") +
  scale_x_continuous("SES Verticality", expand = c(0,0)) +
  scale_fill_continuous_divergingx("spectral", na.value = NA, limits = col.lims, rev = T, guide = guide_colorbar(title = "SES Verticality")) +
  facet_wrap(vars(taxa), ncol = 1, strip.position = "left") +
  coord_sf(crs = "epsg:4326") +
  thm + theme(
              legend.position = "bottom",
              legend.title.position = "top",
              strip.background = element_blank(),
              strip.text = element_blank())
p2.legend = ggpubr::get_legend(p2)
p2.legend = as_ggplot(p2.legend)
p2 = p2 + theme(legend.position = "none")

cols = rgbcols$key$data
colt = cols[which(cols$p1==1), "rgb"]
colf = cols[which(cols$p2==1), "rgb"]
cola = cols[which(cols$p3==1), "rgb"]


p3 = all %>% 
  dplyr::select(p.arb, p.fos, p.ter, taxa, y) %>% 
  rename(Arboreal = p.arb, Fossorial = p.fos, Terrestrial = p.ter) %>% 
  pivot_longer(1:3, names_to = "micro", values_to = "prop") %>% 
  group_by(y, micro, taxa) %>% 
  summarize(prop.mean = mean(prop), prop.sd = sd(prop)) %>% 
  mutate(micro = factor(micro, levels = c("Arboreal", "Terrestrial", "Fossorial"))) %>% 
  #plot
  ggplot() +
  geom_line(aes(x = y, y = prop.mean, color = micro), size=1) +
  geom_ribbon(aes(x = y, ymin = prop.mean - prop.sd, ymax = prop.mean + prop.sd, fill = micro), alpha = 0.2) +
  scale_y_continuous("Proportion", breaks = seq(0,1,0.5)) +
  scale_x_continuous("Latitude", breaks = seq(-60,60,30)) +
  #scale_color_discrete_divergingx("zissou1", guide = guide_legend(title = ""), rev = T) +
  scale_color_manual(values = c(cola, colt, colf), guide = guide_legend(title = "")) +
  scale_fill_manual(values = c(cola, colt, colf), guide = guide_legend(title = "")) +
  #scale_fill_discrete_divergingx("zissou1", guide = guide_legend(title = ""), rev = T) +
  facet_rep_wrap(vars(taxa), ncol = 1, strip.position = "left") +
  coord_flip(xlim = c(ymin,ymax), clip = "off") +
  theme(axis.title.y = element_blank(),
        legend.position = "none",
        strip.placement = "outside",
        strip.background = element_rect(color = NA, fill = NA),
        axis.line = element_line(color = "black"),
        panel.background = element_blank(),
        strip.text = element_text(size = 10))


p2.grob = ggplotGrob(p2)
p.grob = ggplotGrob(p)
p3.grob = ggplotGrob(p3)

max_height <- unit.pmax(p2.grob$heights, p.grob$heights, p3.grob$heights)
p2.grob$heights <- max_height
p.grob$heights <- max_height
p3.grob$heights <- max_height
p3.grob$widths = p3.grob$widths
p.grob$widths = 2*p.grob$widths
p2.grob$widths = 2*p2.grob$widths

# Use grid.arrange to combine them while preserving aspect ratio for the maps
grid.newpage()
pfinal = cbind(p3.grob, p.grob, p2.grob,  size = "first")
grid.draw(pfinal)

svg("figures/main_figs/fig2_vert_maps/sesvert_legend.svg", width = 6, height = 2)
p2.legend
dev.off()

svg("figures/main_figs/fig2_vert_maps/rgb_legend.svg", width = 3, height = 3)
rgbcols$key +
  theme(axis.text = element_text(size = 20))
dev.off()

svg("figures/main_figs/fig2_vert_maps/fig.svg", width = 7, height = 6)
grid.draw(pfinal)
dev.off()

# figure finished in inkscape

library(rphylopic)
squirrel = get_uuid("Sciurus vulgaris")
squirrel = get_phylopic(squirrel)
save_phylopic(squirrel, path = "figures/main_figs/fig2_vert_maps/squirrel.svg")

lizard = pick_phylopic("Coleonyx fasciatus")
save_phylopic(lizard, path = "figures/main_figs/fig2_vert_maps/lizard.svg")


frog = pick_phylopic("Plectrohyla avia")
save_phylopic(frog, path = "figures/main_figs/fig2_vert_maps/frog.svg")

bird = pick_phylopic("Cuculidae")
save_phylopic(bird, path = "figures/main_figs/fig2_vert_maps/bird.svg")
