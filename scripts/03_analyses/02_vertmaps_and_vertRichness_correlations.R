library(tidyverse)
library(png)
library(ggplot2)
library(patchwork)
library(terra)
library(data.table)
library(tidyterra)
library(colorspace)
library(grid)
#library(tricolore)
library(lemon)
library(cowplot)
library(gridExtra)


# load data ----------------------------------------------------

wd = vect("data/original/rnaturalearth_world.shp") %>% 
  crop(ext(-180,180,-60,83.64)) %>% 
  project("+proj=robin +datum=WGS84") 

ymin = ext(wd)[3]
ymax = ext(wd)[4]

amph = fread("data/derivative_data/gridcell_data/env_forest/50_km/amph_comdat.csv") %>% 
  dplyr::select(x, y, rich, p.arb, p.fos, p.ter, vert.mean, vert.mean.ses, biome) %>% 
  filter(rich >= 5) %>% 
  filter(!is.na(vert.mean.ses)) %>% 
  drop_na() %>% 
  rast(crs = "+proj=cea +datum=WGS84") %>% 
  project("+proj=robin +datum=WGS84") %>% 
  as.data.frame(xy = T) %>% 
  mutate(taxa = "Amphibians")

mammals = fread("data/derivative_data/gridcell_data/env_forest/50_km/mammals_comdat.csv") %>% 
  dplyr::select(x, y, rich, p.arb, p.fos, p.ter, vert.mean, vert.mean.ses, biome) %>% 
  filter(rich >= 5) %>% 
  filter(!is.na(vert.mean.ses)) %>% 
  drop_na() %>% 
  rast(crs = "+proj=cea +datum=WGS84") %>% 
  project("+proj=robin +datum=WGS84") %>% 
  as.data.frame(xy = T) %>% 
  mutate(taxa = "Mammals")

birds = fread("data/derivative_data/gridcell_data/env_forest/50_km/birds_comdat.csv") %>% 
  dplyr::select(x, y, rich, p.arb, p.fos, p.ter, vert.mean, vert.mean.ses, biome) %>% 
  filter(rich >= 5) %>% 
  filter(!is.na(vert.mean.ses)) %>% 
  drop_na()%>% 
  rast(crs = "+proj=cea +datum=WGS84") %>% 
  project("+proj=robin +datum=WGS84") %>% 
  as.data.frame(xy = T) %>%  
  mutate(taxa = "Birds")

repts = fread("data/derivative_data/gridcell_data/env_forest/50_km/reptiles_comdat.csv") %>% 
  dplyr::select(x, y, rich, p.arb, p.fos, p.ter, vert.mean, vert.mean.ses, biome) %>% 
  filter(rich >= 5) %>% 
  filter(!is.na(vert.mean.ses)) %>% 
  drop_na()%>% 
  rast(crs = "+proj=cea +datum=WGS84") %>% 
  project("+proj=robin +datum=WGS84") %>% 
  as.data.frame(xy = T) %>%  
  mutate(taxa = "Reptiles")

all = bind_rows(amph, mammals, birds, repts) %>% 
  drop_na(vert.mean.ses) %>% 
  mutate(taxa = factor(taxa, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")))
min = min(all$vert.mean.ses, na.rm = T)
max = max(all$vert.mean.ses, na.rm = T)
col.lims = range(all$vert.mean.ses, na.rm = T)



# PLOT VERTICALITY --------------------------------------------------------

vertmaps = list()
for(i in levels(all$taxa)) {
  
  d = all %>% filter(taxa == i)
  
  map = ggplot() +
    geom_spatvector(data = wd, fill = "black", color = NA) +
    geom_tile(data = d, aes(x,y,fill = vert.mean.ses)) +
    #geom_spatvector(data = wd, fill = NA, color = "black") +
    scale_x_continuous("SES Verticality", expand = c(0,0)) +
    scale_fill_continuous_divergingx("spectral", na.value = NA, limits = col.lims,
                                     rev = T, guide = guide_colorbar(title = "Verticality")) +
    #facet_wrap(vars(taxa), ncol = 1, strip.position = "left") +
    coord_sf(crs = "+proj=robin +datum=WGS84") +
    theme_void() +
    theme(strip.text = element_blank(),
          legend.position = "bottom",
          legend.title.position = "top",
          legend.title = element_text(hjust = 0.5, size = 10),
          legend.key.height = unit(3, units = "mm"),
          legend.key.width = unit(10, units = "mm"),
          axis.title = element_blank(),
          axis.line = element_line(color = "white"),
          axis.text = element_blank(),
          axis.ticks = element_line(color = "white"),
          panel.background = element_blank())
  
  vert.legend = ggpubr::get_legend(map)
  #vert.legend = as_ggplot(vert.legend)
  vertmaps[[i]] = map + theme(legend.position = "none")

}


# PLOT VERT-RICHNESS CORRELATIONS -------------------------------------------------------

corplot = list()
cor = list()

for(i in levels(all$taxa)) {
  d = all %>% filter(taxa==i)
  
  cor[[i]] = cor.test(d$vert.mean.ses, d$rich)
  r = round(cor[[i]]$estimate,2)
  
  p = ggplot(d, aes(vert.mean.ses, rich)) +
    geom_point(pch = ".", alpha = 0.05) +
    geom_smooth(method = "lm", color = "cyan3", linewidth = 0.5) +
    annotate(geom = "text", x = -Inf, y = Inf, label = paste0("r = ", r),
             hjust = -0.1, vjust = 1.2, size = 3) +
    #facet_rep_wrap(~taxa, scales = "free_y", ncol = 1) +
    scale_y_continuous("Richness", position = "right") +
    scale_x_continuous("Verticality") +
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.text = element_blank())
  
  corplot[[i]] = p
  
}


# PLOT PROPORTIONS ACROSS LATITUDES ---------------------------------------

pal = divergingx_hcl(palette = "zissou", n = 3)
pal2 = divergingx_hcl(palette = "spectral", n =19)

cola = pal2[2]
colt = pal[2]
colf = pal2[17] 


propplots = list()

for(i in levels(all$taxa)) {
  p = all %>% 
    filter(taxa == i) %>% 
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
    scale_x_continuous("Latitude", breaks = c(0)) +
    #scale_color_discrete_divergingx("zissou1", guide = guide_legend(title = ""), rev = T) +
    scale_color_manual(values = c(cola, colt, colf), guide = guide_legend(title = "", direction = "horizontal")) +
    scale_fill_manual(values = c(cola, colt, colf), guide = guide_legend(title = "", direction = "horizontal")) +
    #scale_fill_discrete_divergingx("zissou1", guide = guide_legend(title = ""), rev = T) +
    #facet_rep_wrap(vars(taxa), ncol = 1, strip.position = "left") +
    coord_flip(xlim = c(ymin,ymax), clip = "off") +
    theme(legend.position = "bottom",
          legend.margin = margin(0,0,0,0))
  
  prop_legend = ggpubr::get_legend(p)
  
  p = p +
    theme(legend.position = "none",
          strip.placement = "outside",
          strip.background = element_rect(color = NA, fill = NA),
          axis.line = element_line(color = "black"),
          panel.background = element_blank(),
          strip.text = element_text(size = 10))
  
  propplots[[i]] = p
  
}



# Combine plots -----------------------------------------------------------

library(gridExtra)
library(grid)
library(gtable)

tight_theme <- theme(
  plot.margin = margin(0, 0, 0, 0),
  panel.spacing = unit(0, "pt"),  # for facet_wrap/grid
  #axis.ticks.length = unit(1, "pt"),
  #strip.placement = "outside"
)

# Apply to all plots
propplots <- lapply(propplots, `+`, tight_theme)
for(i in 1:4) {
  if(i!=4){propplots[[i]] = propplots[[i]]+theme(axis.title.x = element_blank())}
  propplots[[i]] = propplots[[i]] + theme(plot.margin = margin(0,0,0,10),
                                          axis.title.y = element_text(margin = margin(r=0)))
  }
vertmaps  <- lapply(vertmaps,  `+`, tight_theme)
# fix aspect ratio for combining plots
vertmaps_fixed <- lapply(vertmaps, function(p) {
  ggdraw(p) + theme(plot.margin = margin(0,0,0,0))
})
#corplot  <- lapply(corplot,  `+`, tight_theme)
for(i in 1:4) {
  p = corplot[[i]]
  p = p + theme(plot.margin = margin(0,20,0,0),
                axis.title.y = element_text(margin = margin(l=25)),
                axis.text.y = element_text(margin = margin(r=25)))
  if(i!=4) {p = p+theme(axis.title.x = element_blank())}
  corplot[[i]]=p
}

rows = list()
icons = list.files("figures/icons/", pattern = ".png", full.names = T)
names(icons) = basename(icons) %>%  gsub(".png", "", .)
names(icons) = c("Birds", "Amphibians", "Reptiles", "Mammals")

for (i in levels(all$taxa)) {
  # Convert each ggplot to grobs
  g_prop = ggplotGrob(propplots[[i]])
  g_map  = ggplotGrob(vertmaps_fixed[[i]])
  g_cor  = ggplotGrob(corplot[[i]])
  
  # Match heights to avoid vertical misalignment
  max_height = grid::unit.pmax(g_prop$heights, g_map$heights, g_cor$heights)
  g_prop$heights <- g_map$heights <- g_cor$heights <- max_height
  g_map$widths <- g_map$widths*2
  g_prop$widths <- g_prop$widths
  g_cor$widths <- g_cor$widths
  
  # Combine into one row (left to right: prop | map | cor)
  row = gtable_cbind(g_prop, g_map, g_cor, size = "first")
  
  # add icon
  img = readPNG(icons[[i]])
  #icon_grob <- rasterGrob(img, width = unit(1, "cm"), height = unit(1, "cm"), just = "center")
  icon_grob <- rasterGrob(
    img,
    x = unit(1.2, "npc"),  # move right within the spanning cell
    y = unit(0.2, "npc"),
    just = "left",
    width = unit(1, "cm"),
    height = unit(1, "cm")
  )

  # Determine layout position (first plot panel = col A, map = col B)
  panel_cols <- which(row$layout$name == "panel")

  # Mid-column index between prop and map
  insert_col <- panel_cols[1] + 1  # place between 1st and 2nd plot
  insert_row <- mean(range(row$layout$t[panel_cols]))  # vertically center

  row_with_icon <- gtable_add_grob(
    row,
    grobs = icon_grob,
    t = insert_row,
    l = insert_col,
    r = insert_col + 1,
    name = "icon",
    z = Inf  # put on top
  )
  row_with_icon$layout$clip <- "off"

  
  rows[[i]] = row_with_icon
}

# Combine all rows vertically
final_gtable = do.call(gtable_rbind, rows)

library(grid)

tags <- list(
  textGrob("a", x = 0, just = "left", gp = gpar(fontsize = 13)),
  textGrob("b", x = 0, just = "left", gp = gpar(fontsize = 13)),
  textGrob("c", x = 0, just = "left", gp = gpar(fontsize = 13))
)

# Create row for tags
tag_row <- gtable(
  widths = final_gtable$widths,
  heights = unit(1.5, "lines")
)

# Assuming your columns are: 1–13 (col 1), 14–26 (col 2), 27–39 (col 3)
tag_row <- gtable_add_grob(tag_row, tags[[1]], t = 1, l = 1,  r = 13)
tag_row <- gtable_add_grob(tag_row, tags[[2]], t = 1, l = 14, r = 26)
tag_row <- gtable_add_grob(tag_row, tags[[3]], t = 1, l = 27, r = 39)

final_gtable <- gtable_rbind(tag_row, final_gtable)

# add legend

# Estimate how many columns legends will span
total_cols <- length(final_gtable$widths)
legend_span <- 30  # columns needed for both legends side by side
start_col <- floor((total_cols - legend_span) / 2) + 1
end_col <- start_col + legend_span - 1

# Create legend container
legend_gtable <- gtable(widths = final_gtable$widths, heights = unit(1, "null"))

# Combine legends into one row
combined_legends <- arrangeGrob(
  grobs = list(prop_legend, vert.legend),
  ncol = 2,
  widths = c(1, 1)
)

legend_gtable <- gtable_add_grob(
  legend_gtable,
  grobs = combined_legends,
  t = 1, l = start_col, r = end_col
)

# Combine plot and centered legend
full_plot <- gtable_rbind(final_gtable, legend_gtable)


# Draw the full aligned plot
png("figures/ms_figures/fig2_vertmaps.png", width = 150, height = 150, units = "mm", res = 300)
grid.newpage()
grid.draw(full_plot)
dev.off()












