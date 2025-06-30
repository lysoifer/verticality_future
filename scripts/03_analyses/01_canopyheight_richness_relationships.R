# Model relationships between canopy height and species richness

library(DHARMa)
library(dplyr)
library(foreach)
library(terra)
library(tidyterra)
library(glmmTMB)
library(ggplot2)
library(patchwork)

amph = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/amph_comdat.csv")
birds = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/birds_comdat.csv")
rept = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/reptiles_comdat.csv")
mammals = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/mammals_comdat.csv")

amph$taxa = "Amphibians"
birds$taxa = "Birds"
rept$taxa = "Reptiles"
mammals$taxa = "Mammals"

df = rbind(amph, birds, rept, mammals)

# identify tropical, temperate, and boreal latitudes
tropical = range(abs(df$y[grepl("Tropical", df$biome)]))
temperate = range(abs(df$y[grepl("Temperate", df$biome)]))
boreal = range(abs(df$y[grepl("Boreal|Tundra", df$biome)]))

# Classify biomes more coarsely
# if not forest and not temperate, then tropical woodland, savanna, shrubland
# if not forest and not tropical or boreal, then temperate woodland, savanna, shrubland
# if tropical and forest then tropical forest
# if temperate and forest then temperate and forest
# if boreal/tundra, then boreal
# if Mediterranean and forest, then temperate forest
# if Mediterranean and not forest then temperate woodland, savanna, shrubland
df = df %>% 
  mutate(biome2 = case_when((!grepl("Forest|Temperate", biome)) & (abs(df$y) >= tropical[1]) & (abs(df$y) < tropical[2]) ~ "Tropical woodland, savanna, and shrubland",
                            (!grepl("Forest|Tropical|Tundra", biome)) & (abs(df$y) >= temperate[1]) & (abs(df$y) < temperate[2]) ~ "Temperate woodland, savanna, and shrubland",
                            grepl("(?=.*Forest)(?=.*Tropical)", biome, perl=T) ~ "Tropical forest",
                            grepl("(?=.*Forest)(?=.*Temperate)", biome, perl=T) ~ "Temperate forest",
                            grepl("Boreal|Tundra", biome) ~ "Boreal",
                            grepl("Med", biome) & grepl("forest", ecoregion) & abs(df$y) >= temperate[1] & abs(df$y) < temperate[2] ~ "Temperate forest",
                            grepl("Med", biome) & !grepl("forest", ecoregion, ignore.case = T) & abs(df$y) >= temperate[1] & abs(df$y) < temperate[2] ~ "Temperate woodland, savanna, and shrubland",
                            .default = NA))

# test biome2 classifications
test.r = df %>%
  filter(biome2 == "Boreal" & taxa == "Birds") %>% 
  dplyr::select(x,y,canopy_height2) %>% 
  rast(crs = "+proj=cea +datum=WGS84")
unique(df$biome[df$lat=="Tropical"])

# only filter out rich < 5 for verticality analysis not here

# including bare land (canopy height2 is no bare land (vegheight >= 3m) or human modified land,
# canopy height is aggregated land 2023 without masking landcovers)
canopy_height3 = rast("data/derivative_data/resampled_env_rasters_50km/canopy_height_lang2023/canopy_height_agg.tif")
mapa = rast(extent = c(-20592508, 20588492, -5743602, 6573398), crs = "+proj=cea +datum=WGS84")
res(mapa) = 50000
canopy_height3 = project(canopy_height3, "+proj=cea +datum=WGS84")
canopy_height3 = resample(canopy_height3, mapa)
df$canopy_height3 = terra::extract(canopy_height3, df[,c("x", "y")])[,2] # including bare land (2 is not including bare land)

df = df %>% drop_na(canopy_height2)
df = df %>% filter(rich > 0)

# Explore some preliminary plots
ggplot(df, aes(canopy_height2, rich)) +
  geom_point(pch = ".", alpha = 0.2) +
  geom_smooth(method = "lm") +
  facet_wrap(~taxa, scales = "free_y") +
  theme_classic()

# log-log relationship
ggplot(df, aes(log(canopy_height2 +1), log(rich))) +
  geom_point(pch = ".", alpha = 0.2) +
  geom_smooth(method = "lm") +
  facet_wrap(~taxa, scales = "free_y") +
  theme_classic()

ggplot(df, aes(log(canopy_height2 +1), log(rich), color = biome2)) +
  geom_point(pch = ".", alpha = 0.2) +
  geom_smooth(method = "lm") +
  facet_wrap(~taxa, scales = "free_y") +
  theme_classic()

ggplot(df, aes(canopy_height2, rich, color = biome2)) +
  geom_point(pch = ".", alpha = 0.2) +
  geom_smooth(method = "lm") +
  facet_wrap(~taxa, scales = "free_y") +
  theme_classic()

ggplot(df, aes(canopy_height3, rich, color = biome2)) +
  geom_point(pch = ".", alpha = 0.2) +
  geom_smooth(method = "lm") +
  facet_wrap(~taxa, scales = "free_y") +
  theme_classic()

ggplot(df, aes(canopy_height2+1, rich, color = biome2)) +
  geom_point(pch = ".", alpha = 0.2) +
  geom_smooth(method = "glm", formula = y~log(x), method.args = list(family = "poisson")) +
  facet_wrap(~taxa, scales = "free_y") +
  theme_classic()

ggplot(df, aes(canopy_height2, rich, color = biome2)) +
  geom_point(pch = ".", alpha = 0.2) +
  geom_smooth(method = "glm", formula = y~x, method.args = list(family = "poisson")) +
  facet_wrap(~taxa, scales = "free_y") +
  theme_classic()


# Global richness by canopy height ----------------------------------------


# *- model for each taxonomic group ---------------------------------------

global.taxa = list()
global.taxa.pred = list()
global.taxa.mod = list()
for(i in unique(df$taxa)) {
  d = df %>% filter(taxa == i) %>% 
    drop_na(rich, canopy_height2) %>% 
    filter(rich > 0 & canopy_height2 > 0) %>% 
    mutate(log_ch = log(canopy_height2))
  
  # gaussian (shouldn't be a good fit)
  # m1 = glm(rich ~ canopy_height,
  #          data = d)
  # m1.resid = simulateResiduals(m1, plot = T)
  
  # model richness by canopy height
  # Poisson distribution does not fit data well
  # m2 = glmmTMB(rich ~ log(canopy_height2),
  #          data = d,
  #          family = poisson(link = "log"))
  # m2.resid = simulateResiduals(m2, plot = T) # overdispersion
  
  # negative binomial distribution
  # reptiles: somewhat underdispersed
  m3 = glmmTMB(rich ~ log(canopy_height2), data = d, family = nbinom2(link = "log"))
  m3.resid = simulateResiduals(m3, plot = T)
  testDispersion(m3)
  
  # Deviance residuals
  resid.deviance = residuals(m3, type = "deviance")
  d$resid.deviance = resid.deviance
  # d$resid.dharma = residuals(m3.resid)
  d$taxa = i
  
  global.taxa[[i]] = d
  
  newdat = data.frame(canopy_height2 = seq(min(d$canopy_height2), max(d$canopy_height2), 1))
  pred = ggeffects::predict_response(m3, newdat)
  pred$taxa = i
  
  global.taxa.pred[[i]] = pred
  global.taxa.mod[[i]] = m3
}

# global model for each taxonomic group
global.taxa = bind_rows(global.taxa) 
global.taxa.pred = bind_rows(global.taxa.pred)


# Global model all taxa combined ------------------------------------------
# model richness by canopy height for all taxa combined

df_combined = df %>% 
  dplyr::select(x,y,rich, canopy_height2, taxa, biome2) %>% 
  pivot_wider(names_from = taxa, values_from = rich, names_prefix = "richness_") %>% 
  filter(canopy_height2 > 0)
df_combined$total_richness = apply(df_combined[,5:8], 1, sum, na.rm = T)

fit.global.all = glmmTMB(total_richness ~ log(canopy_height2), data = df_combined, 
                         family = nbinom2(link = "log"))
fit.resid = simulateResiduals(fit.global.all, plot = T)
testDispersion(fit.global.all)
resid.deviance = residuals(fit.global.all, type = "deviance")
df_combined$resid.deviance = resid.deviance


# Models by taxon and biome2 ----------------------------------------------
# model richness by canopy height for each taxon in each of the coarser biome/lat classifications

biome.taxa = list()
biome.taxa.pred = list()
biome.taxa.mod = list()
for(i in unique(df$taxa)) {
  for(b in unique(df$biome2)) {
    
    nm = paste(i,b,sep = "_")
    
    # filter to given taxon and biome2
    d = df %>% filter(taxa == i & biome2 == b) %>% 
      drop_na(rich, canopy_height2) %>% 
      filter(rich > 0 & canopy_height2 > 0) %>% 
      mutate(log_ch = log(canopy_height2))
    
    # negative binomial distribution
    # reptiles: somewhat underdispersed
    m3 = glmmTMB(rich ~ log(canopy_height2), data = d, family = nbinom2(link = "log"))
    #m3.resid = simulateResiduals(m3, plot = T)
    #testDispersion(m3)
    
    # checked for Mammals, Amphibians, and Reptiles in boreal biome2
    # does not change parameter estimates, but model convergence no longer an issue
    # so no need to worry about returning model convergence issues for these three
    # m4 = update(m3, control = glmmTMBControl(optimizer = optim, optArgs=list(method = "BFGS")))
    # 
    # m5 = update(m3, family = poisson())
    
    # Deviance residuals
    resid.deviance = residuals(m3, type = "deviance")
    d$resid.deviance = resid.deviance
    # d$resid.dharma = residuals(m3.resid)
    d$taxa = i
    d$biome2 = b
    
    biome.taxa[[nm]] = d
    
    # predict model to range of canopy height
    newdat = data.frame(canopy_height2 = seq(min(d$canopy_height2), max(d$canopy_height2), 1))
    pred = ggeffects::predict_response(m3, newdat)
    # pred = as.data.frame(pred) %>% 
    #   rename(canopy_height2 = x, predicted_rich = predicted)
    pred$taxa = i
    pred$biome2 = b
    
    biome.taxa.pred[[nm]] = pred
    biome.taxa.mod[[nm]] = m3
  }
}

# global model for each taxonomic group
biome.taxa = bind_rows(biome.taxa) 
biome.taxa.pred = bind_rows(biome.taxa.pred)

# checking convergence issues - no need to worry - see notes above
# for(i in 1:length(biome.taxa.mod)) {
#   print(names(biome.taxa.mod)[i])
#   print(performance::check_convergence(biome.taxa.mod[[i]]))
# }



# PLOT RESULTS ------------------------------------------------------------


# Fig. 1 ------------------------------------------------------------------

# Legend: a) canopy height, b) vertebrate richness (all vertebrates combined)
# c) Deviance residuals for a model of total vertebrate richness by canopy height
# d) Relationships between species richness and canopy height for each taxon in each
# biome category. Dots indicate observed values and lines indicate negative binomial
# model fits of richness by canopy height for each combination of taxon and biome.

# *- Total canopy, richness, and resid maps -------------------------------

wd = vect("data/original/rnaturalearth_world.shp") %>% 
  crop(ext(-180,180,-60,90)) %>% 
  project("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m")

# Map of canopy richness, canopy height, and residuals for total_richness model
df_combined.r = df_combined %>% 
  dplyr::select(x,y,canopy_height2, total_richness, resid.deviance) %>% 
  rast(crs = "+proj=cea +datum=WGS84") %>% 
  project("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m")

total_rich.plt = ggplot() + 
  geom_spatvector(data = wd, fill = "black", color = NA) +
  geom_spatraster(data = df_combined.r, aes(fill = total_richness)) +
  scale_fill_continuous_sequential(palette = "Sunset", na.value = NA,
                                   guide = guide_colorbar(title = "")) +
  theme_void() +
  ggtitle("Vertebrate Richness") +
  theme(legend.position = "inside",
        legend.justification = c(0.1,0),
        legend.key.size = unit(2, "mm"),
        legend.text = element_text(size = 4),
        legend.title = element_text(size = 4),
        plot.title = element_text(hjust = 0.5, size = 6))
  
canopy_height.plt = ggplot() + 
  geom_spatvector(data = wd, fill = "black", color = NA) +
  geom_spatraster(data = df_combined.r, aes(fill = canopy_height2)) +
  scale_fill_continuous_sequential(palette = "YlGn", na.value = NA,
                                   guide = guide_colorbar(title = "m")) +
  theme_void() +
  ggtitle("Canopy Height") +
  theme(legend.position = "inside",
        legend.justification = c(0.1,0),
        legend.key.size = unit(2, "mm"),
        legend.text = element_text(size = 4),
        legend.title = element_text(size = 4),
        plot.title = element_text(hjust = 0.5, size = 6))

total_resid.plt = ggplot() + 
  geom_spatvector(data = wd, fill = "black", color = NA) +
  geom_spatraster(data = df_combined.r, aes(fill = resid.deviance)) +
  scale_fill_continuous_divergingx(palette = "Spectral", na.value = NA,
                                   guide = guide_colorbar(title = ""),
                                   rev = T) +
  theme_void() +
  ggtitle("Residuals") +
  theme(legend.position = "inside",
        legend.justification = c(0.1,0),
        legend.key.size = unit(2, "mm"),
        legend.text = element_text(size = 4),
        legend.title = element_text(size = 4),
        plot.title = element_text(hjust = 0.5, size = 6))


total.global.map = canopy_height.plt + total_rich.plt + total_resid.plt +
  plot_layout(nrow = 1)



# *- Biome per taxa response ----------------------------------------------

biome.taxa.pred.df = biome.taxa.pred %>% 
  as.data.frame() %>% 
  mutate(taxa = factor(taxa, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")),
         biome2 = factor(biome2, levels = c("Tropical forest", "Tropical woodland, savanna, and shrubland",
                                            "Temperate forest", "Temperate woodland, savanna, and shrubland",
                                            "Boreal")),
         col = case_when(biome2 == "Tropical forest" ~ "#008234",
                         biome2 == "Tropical woodland, savanna, and shrubland" ~ "#50DE00",
                         biome2 == "Temperate forest" ~ "#E08A00",
                         biome2 == "Temperate woodland, savanna, and shrubland" ~ "#FFCB76",
                         biome2 == "Boreal" ~ "#56B4E9"),
         col = factor(col, levels = c("#008234", "#50DE00","#E08A00", "#FFCB76", "#56B4E9")))

biome.taxa = biome.taxa %>% 
  mutate(taxa = factor(taxa, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")),
         biome2 = factor(biome2, levels = c("Tropical forest", "Tropical woodland, savanna, and shrubland",
                                            "Temperate forest", "Temperate woodland, savanna, and shrubland",
                                            "Boreal")),
         col = case_when(biome2 == "Tropical forest" ~ "#008234",
                         biome2 == "Tropical woodland, savanna, and shrubland" ~ "#50DE00",
                         biome2 == "Temperate forest" ~ "#E08A00",
                         biome2 == "Temperate woodland, savanna, and shrubland" ~ "#FFCB76",
                         biome2 == "Boreal" ~ "#56B4E9"),
         col = factor(col, levels = c("#008234", "#50DE00","#E08A00", "#FFCB76", "#56B4E9")))

img = list.files("figures/icons/", full.names = T, pattern = ".svg")
nms = gsub("svg", "png", img)


for(i in 1:4) {
  rsvg_png(img[i], file = nms[i], width = 300)
}
img = list.files("figures/icons", full.names = T, pattern = ".png")
img = list(img[1], img[4], img[3], img[2])
names(img) = c("Birds", "Mammals", "Reptiles", "Amphibians")
img = lapply(img, png::readPNG)
img = lapply(img, grid::rasterGrob, interpolate = T)

taxa = c("Birds", "Mammals", "Reptiles", "Amphibians")
plts = list()
for(i in taxa) {
  split_df = biome.taxa.pred.df %>% filter(taxa == i)
  split_df2 = biome.taxa %>% filter(taxa == i)
  
  base_plot = ggplot() +
    geom_point(data = split_df2, aes(canopy_height2, rich, color = col), pch = ".", alpha = 0.05) +
    geom_line(data = split_df, aes(x, predicted, color = col), linewidth = 0.5) +
    geom_ribbon(data = split_df, aes(x, ymin = conf.low, ymax = conf.high, fill = col),
                color = NA, alpha = 0.2) +
    scale_x_continuous("Canopy Height (m)") +
    scale_y_continuous("Richness") +
    scale_color_identity(guide = guide_legend(nrow = 2, ncol = 3, byrow = F), 
                         labels = levels(biome.taxa.pred.df$biome2),
                         breaks = levels(biome.taxa.pred.df$col)) +
    scale_fill_identity(guide = guide_legend(nrow = 2, ncol = 3, byrow = F), 
                        labels = levels(biome.taxa.pred.df$biome2),
                        breaks = levels(biome.taxa.pred.df$col)) +
    theme_classic() +
    theme(panel.grid = element_blank(),
          legend.title = element_blank(),
          axis.title.x = element_text(size = 6),
          axis.text = element_text(size = 4),
          legend.text = element_text(size = 5),
          legend.key.size = unit(3, "mm"),
          axis.line = element_line(linewidth = 0.25),
          axis.ticks = element_line(linewidth = 0.1))
  
  legend = ggpubr::get_legend(base_plot)
  
  base_plot = base_plot + theme(legend.position = "none")
  
  plts[[i]] = base_plot
  
  # Create an overlay plot with the image in the top-left
  image_layer <- ggplot() +
    annotation_custom(img[[i]], xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    theme_void()

  # Combine base plot with image overlay
  plts[[i]] = base_plot + inset_element(image_layer, left = 0, bottom = 0.8, right = 0.2, top = 1)
  
}

shared_y_title <- grid::textGrob("Richness", rot = 90, gp = grid::gpar(fontsize = 6))

main_plot <- (wrap_elements(shared_y_title) | plts[[1]] | plts[[2]] | plts[[3]] | plts[[4]]) +
  plot_layout(widths = c(0.05,1,1,1,1)) &
  theme(axis.title.y = element_blank())

main_plot = main_plot/wrap_elements(legend) +
  plot_layout(heights = c(1,0.3))

tags = list(c("a", "b", "c", "d", "", "", "", "", "", "", "", "", ""))
main_plot = (canopy_height.plt + total_rich.plt + total_resid.plt) / main_plot + 
  plot_layout(heights = c(0.8,1)) + 
  plot_annotation(tag_levels = tags) &
  theme(plot.tag.position = c(0.02, 0.98),
        plot.tag = element_text(size = 6))

ggsave("figures/ms_figures/fig1_canopy_richness.png", width = 180, height = 100, units = "mm", dpi = 300)









