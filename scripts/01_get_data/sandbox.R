library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(terra)
library(tidyterra)
library(colorspace)

env = read.csv("data/derivative_data/env_data_50km.csv")
head(env)
wd = vect("data/original/rnaturalearth_world.shp") %>% 
  project("epsg:4326") %>% 
  crop(ext(-180,180,-60,83.64))

env = env %>% drop_na(-fhd.mean)
hist(env$veg_den)
ggplot(env, aes(veg_den, fhd.mean)) +
  geom_point(pch = ".") +
  geom_smooth(method = "lm")
cor.test(env$fhd.mean, env$veg_den)

# forested biomes only
env.for = env %>% 
  filter(grepl("Forest", biome)) %>% 
  filter(!grepl("Med", biome))
env.for.r = rast(env.for, crs = "+proj=cea +datum=WGS84")

ggplot() +
  geom_spatvector(data = wd, fill = "gray", color = NA) +
  geom_spatraster(data = env.for.r, aes(fill = veg_den)) +
  scale_fill_continuous_sequential(palette = "viridis", na.value = NA) +
  theme_void()

# calculate quantiles of vegetation density in forests per realm
forden = env.for %>% 
  group_by(realm) %>% 
  summarise(vegden_min = min(veg_den),
            vegden05 = quantile(veg_den, prob = 0.05),
            vegden10 = quantile(veg_den, prob = 0.1))



# Min vegden --------------------------------------------------------------

env.denmin = env %>% 
  left_join(forden, by = "realm") %>% 
  filter(veg_den >= vegden_min)

env.denmin.r = rast(env.denmin, crs = "+proj=cea +datum=WGS84")

ggplot() +
  geom_spatvector(data = wd, fill = "gray", color = NA) +
  geom_spatraster(data = env.denmin.r, aes(fill = veg_den)) +
  scale_fill_continuous_sequential(palette = "viridis", na.value = NA) +
  theme_void()


# Veg den05 ---------------------------------------------------------------

# filter out veg den below 5th quantile of forest per realm

env.den05 = env %>% 
  left_join(forden, by = "realm") %>% 
  filter(veg_den >= vegden05)

env.den05.r = rast(env.den05, crs = "+proj=cea +datum=WGS84")

ggplot() +
  geom_spatvector(data = wd, fill = "gray", color = NA) +
  geom_spatraster(data = env.den05.r, aes(fill = veg_den)) +
  scale_fill_continuous_sequential(palette = "viridis", na.value = NA) +
  theme_void()


# filter out veg den below 10th quantile of forest per realm

env.den10 = env %>% 
  left_join(forden, by = "realm") %>% 
  filter(veg_den >= vegden10)

env.den10.r = rast(env.den10, crs = "+proj=cea +datum=WGS84")

ggplot() +
  geom_spatvector(data = wd, fill = "gray", color = NA) +
  geom_spatraster(data = env.den10.r, aes(fill = veg_den)) +
  scale_fill_continuous_sequential(palette = "viridis", na.value = NA) +
  theme_void()



# Prelim analysis ---------------------------------------------------------

amph = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/amph_comdat.csv")
birds = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/birds_comdat.csv")
rept = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/reptiles_comdat.csv")
mammals = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/mammals_comdat.csv")

amph$taxa = "Amphibians"
birds$taxa = "Birds"
rept$taxa = "Reptiles"
mammals$taxa = "Mammals"

cl = c("Amphibians", "Birds", "Reptiles", "Mammals")
bclass = c("Tropical", "Temperate", "Boreal")

df = bind_rows(amph, birds, rept, mammals)

dflist = list(env.for, env.denmin, env.den05, env.den10)

out = list()
for(i in 1:length(dflist)) {
  subst = dflist[[i]][,c("x", "y")]
  d = inner_join(df, subst)
  d = d %>% filter(rich > 5)
  
  map = d %>% 
    dplyr::select(x,y,veg_den) %>% 
    rast(crs = "+proj=cea +datum=WGS84")
  mapplot = ggplot() +
    geom_spatvector(data = wd, fill = "grey", color = NA) +
    geom_spatraster(data = map, aes(fill = veg_den)) +
    scale_fill_continuous_sequential(palette = "viridis", na.value = NA) +
    theme_void()
  
  res = list()
  res2 = list()
  vertmods = list()
  vertmods.tidy = list()
  # run global model
  for(c in 1:4) {
    dsub = d %>% filter(taxa == cl[c])
    cfit = glmmTMB::glmmTMB(rich ~ log(canopy_height2 + 1), data = dsub, family = nbinom2())
    dsub$pred = predict(cfit, type = "response")
    res[[cl[c]]] = dsub
    
    # run models in lat bands
    for(b in bclass) {
      nm = paste0(cl[c], "_", b)
      dsub2 = dsub %>% filter(grepl(b, biome))
      cfit2 = glmmTMB::glmmTMB(rich ~ log(canopy_height2 + 1), data = dsub2, family = nbinom2())
      dsub2$pred = predict(cfit2, type = "response")
      dsub2$bclass = b
      res2[[nm]] = dsub2
    }
    
    # model verticality
    print(cl[c])
    dsub.scale = dsub %>% 
      mutate(log_precip_dry = log(precip_dry +1)) %>% 
      mutate_at(.vars = vars(canopy_height:clim_velocity, elev, veg_den, veg_complexity, log_precip_dry, precip_warm, canopy_height2, fhd.mean), .funs = scale)
    
    car::vif(lm(vert.mean.ses ~ canopy_height2 + fhd.mean + tmax_warm + tmin_cold + log_precip_dry + precip_warm, data = dsub.scale))
    
    f1 = formula(vert.mean.ses ~ poly(tmax_warm,2) + tmin_cold +
                   precip_warm + log_precip_dry +
                   canopy_height2 + fhd.mean + fhd.mean:canopy_height2 +
                   precip_warm:canopy_height2 + tmin_cold:canopy_height2 + 
                   poly(tmax_warm,2):canopy_height2 + log_precip_dry:canopy_height2)
    
    m = glmmTMB(f1, data = dsub.scale)
    
    # results
    mtidy = tidy(m, conf.int = T) %>% 
      mutate(sig = ifelse(conf.low * conf.high > 0, "sig", "notsig"),
             taxa = cl[c])
    
    vertmods[[cl[c]]] = m
    vertmods.tidy[[cl[c]]] = mtidy
    
  }
  
  mod.global = bind_rows(res)
  mod.biome = bind_rows(res2)
  vertmods.tidy = bind_rows(vertmods.tidy)
  
  # canopy height by richness
  p1 = ggplot() +
    geom_line(data = mod.global, aes(canopy_height2, pred)) +
    geom_line(data = mod.biome, aes(canopy_height2, pred, color = bclass)) +
    facet_wrap(~taxa, scales = "free_y") +
    theme_bw()
  
  p2 = ggplot(vertmod) +
    geom_pointrange(aes(x = estimate, xmin = conf.low, xmax = conf.high, y = term, color = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    facet_wrap(~taxa) +
    theme_bw()
  
  out[[i]] = list(p1, p2, vertmods, vertmods.tidy, mapplot)
  names(out[[i]]) = c("chplot", "vertplot", "vertmods", "vertmods.tidy", "map")
}

forest = out[[1]]
forest$chplot
forest$vertplot
forest$vertmods.tidy
forest$map
summary(forest$vertmods[[1]])

denmin = out[[2]]
denmin$map
summary(denmin$vertmods[[1]])

vegden05 = out[[3]]
vegden05$chplot
vegden05$vertplot
vegden05$vertmods.tidy
tidy(vegden05$vertmods[[1]])
vegden05$map

vegden10 = out[[4]]
vegden10$chplot
vegden10$vertplot
vegden10$map





chplots = list(forest$chplot, denmin$chplot, vegden05$chplot, vegden10$chplot)
nms = c("forest", "vegden.min", "vegden05", "vegden10")
plts = Map(function(p, title) {
  p + ggtitle(title)
}, chplots, nms)
wrap_plots(plts) + plot_layout(guides = "collect")


vertplots = list(forest$vertplot, denmin$vertplot, vegden05$vertplot, vegden10$vertplot)
plts = Map(function(p, title) {
  p + ggtitle(title)
}, vertplots, nms)
wrap_plots(plts) + plot_layout(guides = "collect", nrow = 4)







