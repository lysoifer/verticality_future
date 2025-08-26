# Model relationships between canopy height and species richness

library(DHARMa)
library(dplyr)
#library(foreach)
library(terra)
library(tidyterra)
library(broom.mixed)
#library(glmmTMB)
library(MASS)
library(ggplot2)
library(patchwork)
library(performance)
library(ggeffects)
library(colorspace)

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
# test.r = df %>%
#   filter(biome2 == "Boreal" & taxa == "Birds") %>% 
#   dplyr::select(x,y,canopy_height2) %>% 
#   rast(crs = "+proj=cea +datum=WGS84")
# unique(df$biome[df$lat=="Tropical"])
# 
# # only filter out rich < 5 for verticality analysis not here
# 
# # including bare land (canopy height2 is no bare land (vegheight >= 3m) or human modified land,
# # canopy height is aggregated land 2023 without masking landcovers)
# canopy_height3 = rast("data/derivative_data/resampled_env_rasters_50km/canopy_height_lang2023/canopy_height_agg.tif")
# mapa = rast(extent = c(-20592508, 20588492, -5743602, 6573398), crs = "+proj=cea +datum=WGS84")
# res(mapa) = 50000
# canopy_height3 = project(canopy_height3, "+proj=cea +datum=WGS84")
# canopy_height3 = resample(canopy_height3, mapa)
# df$canopy_height3 = terra::extract(canopy_height3, df[,c("x", "y")])[,2] # including bare land (2 is not including bare land)

df = df %>% drop_na(canopy_height2)
df = df %>% filter(rich > 0)

# Explore some preliminary plots
# ggplot(df, aes(canopy_height2, rich)) +
#   geom_point(pch = ".", alpha = 0.2) +
#   geom_smooth(method = "lm") +
#   facet_wrap(~taxa, scales = "free_y") +
#   theme_classic()
# 
# # log-log relationship
# ggplot(df, aes(log(canopy_height2 +1), log(rich))) +
#   geom_point(pch = ".", alpha = 0.2) +
#   geom_smooth(method = "lm") +
#   facet_wrap(~taxa, scales = "free_y") +
#   theme_classic()
# 
# ggplot(df, aes(log(canopy_height2 +1), log(rich), color = biome2)) +
#   geom_point(pch = ".", alpha = 0.2) +
#   geom_smooth(method = "lm") +
#   facet_wrap(~taxa, scales = "free_y") +
#   theme_classic()
# 
# ggplot(df, aes(canopy_height2, rich, color = biome2)) +
#   geom_point(pch = ".", alpha = 0.2) +
#   geom_smooth(method = "lm") +
#   facet_wrap(~taxa, scales = "free_y") +
#   theme_classic()
# 
# ggplot(df, aes(canopy_height3, rich, color = biome2)) +
#   geom_point(pch = ".", alpha = 0.2) +
#   geom_smooth(method = "lm") +
#   facet_wrap(~taxa, scales = "free_y") +
#   theme_classic()
# 
# ggplot(df, aes(canopy_height2+1, rich, color = biome2)) +
#   geom_point(pch = ".", alpha = 0.2) +
#   geom_smooth(method = "glm", formula = y~log(x), method.args = list(family = "poisson")) +
#   facet_wrap(~taxa, scales = "free_y") +
#   theme_classic()
# 
# ggplot(df, aes(canopy_height2, rich, color = biome2)) +
#   geom_point(pch = ".", alpha = 0.2) +
#   geom_smooth(method = "glm", formula = y~x, method.args = list(family = "poisson")) +
#   facet_wrap(~taxa, scales = "free_y") +
#   theme_classic()


# Global richness by canopy height ----------------------------------------


# *- model for each taxonomic group ---------------------------------------

global.taxa = list()
global.taxa.pred = list()
global.taxa.mod = list()
global.taxa.pd = list()
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
  m3 = MASS::glm.nb(rich~canopy_height2, data=d)
  m4 =  MASS::glm.nb(rich~log(canopy_height2), data=d)
  
  # calculate percent deviance explained
  mnull = MASS::glm.nb(rich~1, data = d) # intercept-only model
  m3.pd = (1-(m3$deviance/mnull$deviance))*100
  m4.pd = (1-(m4$deviance/mnull$deviance))*100
  

  # m3 is consistently performing better
  AIC(m3, m4)
  m3.pd
  m3.pd
  
  # m3.resid = simulateResiduals(m3, plot = T)
  # testDispersion(m3)
  
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
  global.taxa.pd[[i]] = m3.pd
}

# global model for each taxonomic group
global.taxa = bind_rows(global.taxa) 
global.taxa.pred = bind_rows(global.taxa.pred)

global.taxa.tidy = data.frame()
for(i in 1:length(global.taxa.mod)) {
  m = global.taxa.mod[[i]]
  m.tidy = tidy(m, conf.int = T, conf.level = 0.95) %>% 
    mutate(taxa = names(global.taxa.mod)[i])
  global.taxa.tidy = rbind(global.taxa.tidy, m.tidy)
}

# add percent deviance to tidy df
global.taxa.pd = data.frame(taxa = names(global.taxa.pd), estimate = unlist(global.taxa.pd), term = "CPDE (%)")
global.taxa.tidy = bind_rows(global.taxa.tidy, global.taxa.pd)



# Global model all taxa combined ------------------------------------------
# model richness by canopy height for all taxa combined

df_combined = df %>% 
  dplyr::select(x,y,rich, canopy_height2, taxa, biome2) %>% 
  pivot_wider(names_from = taxa, values_from = rich, names_prefix = "richness_") %>% 
  filter(canopy_height2 > 0)
df_combined$total_richness = apply(df_combined[,5:8], 1, sum, na.rm = T)

fit.global.all1 = MASS::glm.nb(total_richness ~ canopy_height2, data = df_combined)
fit.global.all2 = MASS::glm.nb(total_richness ~ log(canopy_height2), data = df_combined)

AIC(fit.global.all1, fit.global.all2) # rich~canopyheight better (fit.global.all1)
# percent deviance
mnull = MASS::glm.nb(total_richness~1, data = df_combined)
fit.global.all.pd = (1-(fit.global.all1$deviance/mnull$deviance))*100

nd = data.frame(canopy_height2 = c(min(df_combined$canopy_height2), 
                                   seq(1,floor(max(df_combined$canopy_height2)),1)))
fit.global.all.pred = predict_response(fit.global.all1, nd)

# fit.resid = simulateResiduals(fit.global.all1, plot = T)
# testDispersion(fit.global.all1)
resid.deviance = residuals(fit.global.all1, type = "deviance")
df_combined$resid.deviance = resid.deviance

fit.global.tidy = tidy(fit.global.all1, conf.int = T, conf.level = 0.95)
global.pd = data.frame(term = "CPDE (%)", estimate = fit.global.all.pd)
fit.global.tidy = bind_rows(fit.global.tidy, global.pd)


# BIOME INTERACTION MODEL - ALL TAXA COMBINED -----------------------------


df_combined$biome2 = as.factor(df_combined$biome2)

# fit poisson model first
fit1 = glm(total_richness ~ canopy_height2*biome2 - 1 - canopy_height2,
           data = df_combined, family = poisson())
# simulateResiduals(fit1, plot = T) # overdispersed

fit2 = MASS::glm.nb(total_richness ~ canopy_height2*biome2 - 1 - canopy_height2,
                    data = df_combined)
# simulateResiduals(fit2, plot = T) # looks good

fit3 = MASS::glm.nb(total_richness ~ log(canopy_height2)*biome2 - 1 - log(canopy_height2),
                    data = df_combined)
# simulateResiduals(fit3, plot = T) # looks good

AIC(fit2, fit3)

mnull = MASS::glm.nb(total_richness~1, data = df_combined)
fit.biome.all.pd = (1-(fit2$deviance/mnull$deviance))*100

chmax = df_combined %>% 
  group_by(biome2) %>% 
  summarize(chmax = round(max(canopy_height2)))

fit.biome.all.pred = predict_response(fit2, terms = c("canopy_height2 [1:43]", "biome2")) %>% 
  as.data.frame() %>% 
  rename(canopy_height = x, richness = predicted, biome2 = group) %>% 
  left_join(chmax, by = "biome2") %>% 
  filter(canopy_height <= chmax) # filter out predictions above max canopy height per biome


# BIOME INTERACTION MODELS BY TAXON -------------------------------------------------------------

biome.df = list()
biome.mod = list()
biome.pred = list()
biome.pd = list()
df$biome2 = factor(df$biome2)
for(i in unique(df$taxa)) {
  
  d = df %>% filter(taxa == i) %>% 
    drop_na(rich, canopy_height2) %>% 
    filter(rich > 0 & canopy_height2 > 0) %>% 
    mutate(log_ch = log(canopy_height2))
  
  fit1 = glm(rich ~ canopy_height2*biome2 - 1 - canopy_height2, data = d, family = poisson())
  #simulateResiduals(fit1, plot = T)
  
  fit2 = MASS::glm.nb(rich ~ canopy_height2*biome2 - 1 - canopy_height2, data = d)
  #simulateResiduals(fit2, plot = T)
  
  fit3 = MASS::glm.nb(rich ~ log(canopy_height2)*biome2 - 1 - log(canopy_height2), data = d)
  #simulateResiduals(fit3, plot = T)
  
  AIC(fit2, fit3)

  # percent deviance
  mnull = MASS::glm.nb(rich ~ 1, data = d)
  fit2.pd = (1-(fit2$deviance/mnull$deviance))*100
  fit3.pd = (1-(fit3$deviance/mnull$deviance))*100
  fit2.pd
  fit3.pd
  
  biome.pd[[i]] = fit2.pd
  
  # manual prediction (same as ggeffects::predict_response)
  # nd = expand.grid(c(1,10,20), levels(d$biome2), stringsAsFactors = F)
  # names(nd) = c("canopy_height2", "biome2")
  # test = predict(fit2, nd, se.fit = T)
  # nd$pred = test$fit
  # nd$pred.se = test$se.fit
  # nd = nd %>% 
  #   mutate(pred.resp = exp(pred),
  #          conf.high = exp(pred+pred.se*qnorm(0.975)),
  #          conf.low = exp(pred-pred.se*qnorm(0.975)))
  
  chmax = d %>% 
    group_by(biome2) %>% 
    summarize(chmax = round(max(canopy_height2)))
  
  pred = predict_response(fit2, terms = c("canopy_height2 [1:43]", "biome2")) %>% 
    as.data.frame() %>% 
    rename(canopy_height = x, richness = predicted, biome2 = group) %>% 
    mutate(taxa = i) %>% 
    left_join(chmax, by = "biome2") %>% 
    filter(canopy_height <= chmax)
  
  biome.mod[[i]] = fit2
  biome.pred[[i]] = pred
  biome.df[[i]] = d
  
}

biome.pred = bind_rows(biome.pred) %>% 
  mutate(taxa = factor(taxa, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")))


# Result tables -----------------------------------------------------------

# summary table for models of richness by canopy height for each biome
biome.mod.tidy = data.frame()
for(i in 1:length(biome.mod)) {
  nm = names(biome.mod)[i]
  tidymod = broom.mixed::tidy(biome.mod[[i]], conf.int = T, conf.level = 0.95)
  tidymod$taxa = nm
  biome.mod.tidy = rbind(biome.mod.tidy, tidymod)
}

biome.pd = data.frame(taxa = names(biome.pd), estimate = unlist(biome.pd), term = "CPDE (%)")
biome.mod.tidy = bind_rows(biome.mod.tidy, biome.pd)


# biome.taxa.mod.summ %>% 
#   filter(term == "canopy_height2") %>% 
#   ggplot() +
#   geom_pointrange(aes(x = estimate, xmin = conf.low, xmax = conf.high, y = taxa, color = biome2))

# PLOT RESULTS ------------------------------------------------------------

# Fig 1. REVISED

# *- richness by canopy height for all taxa by biome ----------------------

# combine global model df with biome model df for plotting
global.all.pred = fit.global.all.pred %>% 
  as.data.frame() %>% 
  rename(canopy_height = x, richness = predicted, biome2 = group) %>% 
  mutate(biome2 = "Global")

biome.all.pred = fit.biome.all.pred %>% 
  as.data.frame() %>% 
  bind_rows(global.all.pred) %>% 
  mutate(biome2 = factor(biome2, levels = c("Tropical forest", "Tropical woodland, savanna, and shrubland",
                                            "Temperate forest", "Temperate woodland, savanna, and shrubland",
                                            "Boreal", "Global")),
         col = case_when(biome2 == "Tropical forest" ~ "#008234",
                         biome2 == "Tropical woodland, savanna, and shrubland" ~ "#50DE00",
                         biome2 == "Temperate forest" ~ "#E08A00",
                         biome2 == "Temperate woodland, savanna, and shrubland" ~ "#FFCB76",
                         biome2 == "Boreal" ~ "#56B4E9",
                         biome2 == "Global" ~ "black"),
         col = factor(col, levels = c("#008234", "#50DE00","#E08A00", "#FFCB76", "#56B4E9", "black")))

df_combined = df_combined %>% 
  mutate(biome2 = factor(biome2, levels = c("Tropical forest", "Tropical woodland, savanna, and shrubland",
                                            "Temperate forest", "Temperate woodland, savanna, and shrubland",
                                            "Boreal")),
         col = case_when(biome2 == "Tropical forest" ~ "#008234",
                         biome2 == "Tropical woodland, savanna, and shrubland" ~ "#50DE00",
                         biome2 == "Temperate forest" ~ "#E08A00",
                         biome2 == "Temperate woodland, savanna, and shrubland" ~ "#FFCB76",
                         biome2 == "Boreal" ~ "#56B4E9"),
         col = factor(col, levels = c("#008234", "#50DE00","#E08A00", "#FFCB76", "#56B4E9")))


ch_rich.plt = ggplot(biome.all.pred) +
  geom_point(data = df_combined, aes(canopy_height2, total_richness, color = col),
             pch = ".", alpha = 0.05) + 
  geom_line(aes(canopy_height, richness, color = col)) +
  geom_ribbon(aes(canopy_height, ymin = conf.low, ymax = conf.high, fill =col),
              alpha = 0.2) +
  scale_color_identity(guide = guide_legend(nrow = 3, ncol = 2, byrow = T), 
                       labels = levels(biome.all.pred$biome2),
                       breaks = levels(biome.all.pred$col)) +
  scale_fill_identity(guide = guide_legend(nrow = 3, ncol = 2, byrow = T), 
                      labels = levels(biome.all.pred$biome2),
                      breaks = levels(biome.all.pred$col)) +
  scale_x_continuous("Canopy Height (m)") +
  scale_y_continuous("Richness") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.title = element_text(size = 6),
        axis.text = element_text(size = 4),
        legend.text = element_text(size = 5),
        legend.key.size = unit(3, "mm"),
        axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.1),
        legend.position = "bottom",
        plot.margin = margin(0,0,0,0))

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
        plot.title = element_text(hjust = 0.5, size = 6),
        plot.margin = margin(1,0,1,0))

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
        plot.title = element_text(hjust = 0.5, size = 6),
        plot.margin = margin(1,0,1,0))

total_resid.plt = ggplot() + 
  geom_spatvector(data = wd, fill = "black", color = NA) +
  geom_spatraster(data = df_combined.r, aes(fill = resid.deviance)) +
  scale_fill_continuous_divergingx(palette = "Spectral", na.value = NA,
                                   guide = guide_colorbar(title = ""),
                                   rev = T) +
  # annotate(geom = "text", x = Inf, y = -Inf, hjust = 1.7, vjust = -0.5,
  #          label = paste0("R\u00b2 = ", round(global.r2[[1]], 2)), size = 2) +
  theme_void() +
  ggtitle("Residuals") +
  theme(legend.position = "inside",
        legend.justification = c(0.1,0),
        legend.key.size = unit(2, "mm"),
        legend.text = element_text(size = 4),
        legend.title = element_text(size = 4),
        plot.title = element_text(hjust = 0.5, size = 6),
        plot.gin = margin(1,0,1,0))


total.global.map = canopy_height.plt / total_rich.plt / total_resid.plt +
  plot_layout(nrow = 3) & theme(plot.tag.position = c(0.1, 0.9),
                                plot.tag = element_text(size = 6))

ch_rich.plt  = ch_rich.plt + theme(plot.tag.position = c(0.03, 0.95),
                                   plot.tag = element_text(size = 6))

(ch_rich.plt | total.global.map + 
    plot_layout(widths = c(1,0.1))) +
  plot_annotation(tag_levels = "a")

ggsave('figures/ms_figures/fig1-2_canopy_richness.png', width = 140, height = 100,
       units = "mm", dpi = 300)


# Fig. 1 OLD ------------------------------------------------------------------

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
  # annotate(geom = "text", x = Inf, y = -Inf, hjust = 1.7, vjust = -0.5,
  #          label = paste0("R\u00b2 = ", round(global.r2[[1]], 2)), size = 2) +
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

biome.pred.df = biome.pred %>% 
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

biome.df = biome.df %>% 
  bind_rows() %>% 
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


# for(i in 1:4) {
#   rsvg_png(img[i], file = nms[i], width = 300)
# }
img = list.files("figures/icons", full.names = T, pattern = ".png")
img = list(img[1], img[4], img[3], img[2])
names(img) = c("Birds", "Mammals", "Reptiles", "Amphibians")
img = lapply(img, png::readPNG)
img = lapply(img, grid::rasterGrob, interpolate = T)

taxa = c("Birds", "Mammals", "Reptiles", "Amphibians")
plts = list()
for(i in taxa) {
  split_df = biome.pred.df %>% filter(taxa == i)
  split_df2 = biome.df %>% filter(taxa == i)
  
  base_plot = ggplot() +
    geom_point(data = split_df2, aes(canopy_height2, rich, color = col), pch = ".", alpha = 0.05) +
    geom_line(data = split_df, aes(canopy_height, richness, color = col), linewidth = 0.5) +
    geom_ribbon(data = split_df, aes(canopy_height, ymin = conf.low, ymax = conf.high, fill = col),
                color = NA, alpha = 0.2) +
    # annotate(geom = "text", label = paste0("R\u00b2 = ", round(biome.r2[[i]][[1]],2)),
    #          x = -Inf, y = Inf, hjust = -0.1, vjust = 4.9, size = 2) +
    scale_x_continuous("Canopy Height (m)") +
    scale_y_continuous("Richness") +
    scale_color_identity(guide = guide_legend(nrow = 2, ncol = 3, byrow = F), 
                         labels = levels(biome.pred.df$biome2),
                         breaks = levels(biome.pred.df$col)) +
    scale_fill_identity(guide = guide_legend(nrow = 2, ncol = 3, byrow = F), 
                        labels = levels(biome.pred.df$biome2),
                        breaks = levels(biome.pred.df$col)) +
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


# GLOBAL TOTAL VERTEBRATE PREDICTIONS -------------------------------------

ggplot() +
  geom_point(data = df_combined, aes(canopy_height2, total_richness), pch = ".", alpha = 0.2) +
  geom_line(data = fit.global.all.pred, aes(x, predicted), color = "skyblue2") +
  geom_ribbon(data = fit.global.all.pred, aes(x, ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "skyblue2", color = NA) +
  scale_x_continuous("Canopy Height (m)") +
  scale_y_continuous("Vertebrate Richness") +
  theme_classic()
ggsave("figures/ms_figures/supp_figs/total_richness~canopy.png", width = 90, height = 90, dpi = 300, units = "mm")



# RESIDUAL MAPS PER SPECIES -----------------------------------------------

global.taxa.map = global.taxa %>% 
  dplyr::select(x,y,resid.deviance, taxa) %>% 
  pivot_wider(names_from = "taxa", values_from = "resid.deviance") %>% 
  rast(crs = "+proj=cea +datum=WGS84") %>% 
  project("+proj=robin +datum=WGS84")

df$taxa = factor(df$taxa, levels = c("Birds", "Mammals", "Reptiles", "Amphibians"))

limits = range(global.taxa$resid.deviance)

taxa.resid.plts = list()
for(i in levels(df$taxa)) {
  p1 = ggplot() +
    geom_spatvector(data = wd, fill = "black", color = NA) +
    geom_spatraster(data = global.taxa.map[[i]]) +
    scale_fill_continuous_divergingx(palette = "Spectral", na.value = NA,
                                     guide = guide_colorbar(title = "Deviance Residuals"),
                                     rev = T, limits = limits) +
    theme_void() +
    theme(legend.position = "bottom",
          legend.title.position = "top",
          legend.title = element_text(hjust = 0.5, size = 8),
          legend.key.width = unit(10, "mm"),
          legend.key.height = unit(2, "mm"),
          legend.text = element_text(size = 6))
  
  leg = ggpubr::get_legend(p1)
  
  p1 = p1 + theme(legend.position="none")
  
  # Create an overlay plot with the image in the top-left
  image_layer <- ggplot() +
    annotation_custom(img[[i]], xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    theme_void()
  
  # Combine base plot with image overlay
  taxa.resid.plts[[i]] = p1 + inset_element(image_layer, left = 0.1,
                                            bottom = 0.1, right = 0.35, top = 0.35)
  
}

tags = list(c("a", "", "b", "", "c", "", "d", "", ""))
wrap_plots(taxa.resid.plts) / leg + 
  plot_layout(heights = c(1,0.1)) +
  plot_annotation(tag_levels = tags) &
  theme(plot.tag.position = c(0.1, 0.95),
        plot.tag = element_text(size = 8))
ggsave("figures/ms_figures/supp_figs/deviance_resids_taxa.png", width = 180, height = 100, 
       dpi = 300, unit= "mm")

# RATES OF INCREASE IN RICHNESS -------------------------------------------

ratedifs = biome.pred %>% 
  as.data.frame() %>% 
  dplyr::select(canopy_height, richness, biome2, taxa) %>% 
  filter(canopy_height == 15 | canopy_height == 25) %>% 
  pivot_wider(names_from = canopy_height, values_from = richness, names_prefix = "h") %>% 
  mutate(dif = h25-h15)

# maximum richness in a cell per biome/taxa category
# may be better to have max richness per biome/taxa category, but don't have that available right now
# would have to recalculate from range maps
maxrichcell = df %>% 
  group_by(taxa, biome2) %>% 
  summarise(maxrich = max(rich))

ratedifs = left_join(ratedifs, maxrichcell, by = c("taxa", "biome2")) %>% 
  mutate(relative_increase = dif/maxrich * 100,
         percent_increase = dif/h15*100)





# SCRATCH SPACE -----------------------------------------------------------

# df2 = df %>% filter(rich >= 5)
# 
# amph2 = df2 %>% filter(taxa == "Amphibians")
# 
# m1 = lm(rich ~ vert.mean.ses, data = amph2)
# m1.resid = resid(m1)
# plot(amph2$canopy_height2, m1.resid)
# plot(amph2$precip_dry, m1.resid)
# ggplot(amph2, aes(log(precip_dry +1), m1.resid)) + geom_point() + geom_smooth(method = "lm")
# ggplot(amph2, aes(canopy_height2, m1.resid)) + geom_point() + geom_smooth(method = "lm")
# ggplot(amph2, aes(tmax_warm, m1.resid)) + geom_point() + geom_smooth(method = "lm")
# ggplot(amph2, aes(tmin_cold, m1.resid)) + geom_point() + geom_smooth(method = "lm")
# ggplot(amph2, aes(vert.mean.ses, rich)) + 
#   geom_point(alpha = 0.05) + 
#   geom_smooth(method = "glm",  method.args = list(family = "poisson"))
# 
# # so it looks like precip and canopy height drive verticality, which influences richness
# # therefore, verticality mediates relationships between canopy height and richness
# # and factors limiting verticality will limit richness
# # but temperature influences richness independently of verticality
# 
# m2 = lm(rich ~ vert.mean.ses + canopy_height2, data = amph2)
# summary(m2)
# 
# m3 = lm(rich ~ canopy_height2, data = amph2)
# summary(m3)
# 
# m1 = MASS::glm.nb(rich ~ vert.mean.ses, data = amph2)
# m1.resid = residuals(m1, type = "deviance")
# amph2$m1.resid = m1.resid
# ggplot(amph2, aes(canopy_height2, m1.resid)) +
#   geom_point(alpha = 0.05) +
#   geom_smooth(method = "lm")
# 
# mod.resid = lm(m1.resid ~ canopy_height2, data = amph2)
# summary(mod.resid)
# 
# cor.test(m1.resid, amph2$canopy_height2)
# 
# 
# 
# 
