source("scripts/00_functions/00_functions.R")
library(nlme)
library(car)
library(ncf)
library(sdmTMB)
library(lme4)
library(tidyverse)
library(terra)


options(na.action = "na.omit")


# VIF variable selection
#env = read.csv("data/derivative_data/env_data_50km_forest.csv")

# amphibian forest only 50km resolution
dat = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/amph_comdat.csv")

# make key value dataframe for categorical rasters
key.ecoregion = data.frame(ecoregion = unique(dat$ecoregion))
key.ecoregion$key.ecoregion = 1:nrow(key.ecoregion)

key.biome = data.frame(biome = unique(dat$biome))
key.biome$key.biome = 1:nrow(key.biome)

key.realm = data.frame(realm = unique(dat$realm))
key.realm$key.realm = 1:nrow(key.realm)

datr = dat  %>% relocate(x,y, .before = rich) %>% 
  left_join(key.biome, by = "biome") %>% 
  left_join(key.realm, by = "realm") %>% 
  left_join(key.ecoregion, by = "ecoregion") %>% 
  rast(type = "xyz", crs = "+proj=cea +datum=WGS84")
plot(datr$vert.mean)
plot(datr$vert.mean.ses)

# variable selection 
#env_vars = vif_func(in_frame = env[c(4, 7:12,17:18,22)], thresh = 5, trace = T)

# randomly select 5000
set.seed(12345)
d = dat %>% 
  # subset to richness >= 5
  filter(rich >= 5) %>% 
  filter(!is.na(vert.mean.ses)) %>% 
  drop_na() %>% 
  slice_sample(n=5000) %>% 
  mutate_at(.vars = vars(canopy_height:clim_velocity, elev, veg_den, veg_complexity), .funs = scale)

# look at the data
# plot vert.mean.ses
d %>% 
  dplyr::select(vert.mean.ses, biome:clim_velocity, elev, veg_den, veg_complexity) %>%
  pivot_longer(cols = 3:16, names_to = "var", values_to = "val") %>% 
  ggplot(aes(x = val, y = vert.mean.ses, color = biome)) +
    geom_point(pch = ".") +
    facet_wrap(~var, scales = "free") +
    theme_classic()

# plot vert.mean.ses
d %>% 
  dplyr::select(vert.mean.ses, biome, canopy_height:clim_velocity, elev, veg_den, veg_complexity) %>%
  pivot_longer(cols = 3:16, names_to = "var", values_to = "val") %>% 
  ggplot(aes(x = val, y = vert.mean.ses, color = biome, fill = biome)) +
  geom_point(pch = ".") +
  geom_smooth(method = "lm") +
  facet_wrap(~var, scales = "free") +
  theme_classic()

d %>% 
  dplyr::select(vert.mean, biome:clim_velocity, elev, veg_den, veg_complexity) %>%
  pivot_longer(cols = 3:16, names_to = "var", values_to = "val") %>% 
  ggplot(aes(x = val, y = vert.mean, color = biome)) +
  geom_point(pch = ".") +
  facet_wrap(~var, scales = "free") +
  theme_classic()

d %>% 
  dplyr::select(vert.mean, biome:clim_velocity, elev, veg_den, veg_complexity) %>%
  pivot_longer(cols = 3:16, names_to = "var", values_to = "val") %>% 
  ggplot(aes(x = val, y = vert.mean, color = biome, fill = biome)) +
  geom_point(pch = ".") +
  geom_smooth(pch = ".") +
  facet_wrap(~var, scales = "free") +
  theme_classic()

# model selection
# H1: structural hypothesis - vegetation structure drives patterns of verticality
d %>% dplyr::select(canopy_height, veg_den, veg_complexity) %>% pairs()

# test VIF for all variables
vif(lm(vert.mean.ses ~ canopy_height + veg_den + veg_complexity, data = d))
# veg complexity has the highest vif and = canopy_height * veg_den, so remove veg_complexity
vif(lm(vert.mean.ses ~ canopy_height + veg_den, data = d))
# vif < 3 for all remaining variables

mod.lm = gls(vert.mean.ses ~ canopy_height + veg_den, data = d, method = "REML")
resid.lm = residuals(mod.lm, type = "normalized")
cor.test = data.frame(x = d$x, y = d$y, resid = resid.lm) %>% 
  slice_sample(n=1000)
cor.lm = ncf::correlog(x = cor.test$x, y = cor.test$y, z = cor.test$resid, increment = 100000, resamp = 99, latlon = F)
ncf:::plot.correlog(cor.lm, xlim = c(0,4e6))

# add biome as a random effect
mod.biome = lme(vert.mean.ses ~ canopy_height + veg_den, random = ~1 | biome, data = d, method = "REML")
summary(mod.biome)
# look at residuals
resid.biome = residuals(mod.biome, type = "normalized")
cor.test = data.frame(x = d$x, y = d$y, resid = resid.biome) %>% 
  slice_sample(n=1000)
cor.biome = ncf::correlog(x = cor.test$x, y = cor.test$y, z = cor.test$resid, increment = 100000, resamp = 99, latlon = F)
ncf:::plot.correlog(cor.biome, xlim = c(0,4e6))

# add biome as a random effect nested within realm
mod.biorealm = lme(vert.mean.ses ~ canopy_height + veg_den, random = ~1 | realm/biome, data = d, method = "REML")
summary(mod.biorealm)
# look at residuals
resid.biorealm = residuals(mod.biorealm, type = "normalized")
cor.test = data.frame(x = d$x, y = d$y, resid = resid.biorealm) %>% 
  slice_sample(n=1000)
cor.biorealm = ncf::correlog(x = cor.test$x, y = cor.test$y, z = cor.test$resid, increment = 100000, resamp = 99, latlon = F)
ncf:::plot.correlog(cor.biorealm, xlim = c(0,4e6))

# still spatial autocorrelation - try blocking by ecoregion
mod.ecoregion = lme(vert.mean.ses ~ canopy_height + veg_den, random = ~1 | ecoregion, data = d, method = "REML")
summary(mod.ecoregion)
# look at residuals
resid.ecoregion = residuals(mod.ecoregion, type = "normalized")
cor.test = data.frame(x = d$x, y = d$y, resid = resid.ecoregion) %>% 
  slice_sample(n=1000)
cor.ecoregion = ncf::correlog(x = cor.test$x, y = cor.test$y, z = cor.test$resid, increment = 100000, resamp = 99, latlon = F)
ncf:::plot.correlog(cor.ecoregion, xlim = c(0,4e6))




# H2: climate hypothesis
d %>% dplyr::select(tmax_warm, tmin_cold, temp_sea, precip_wet, precip_dry, precip_sea) %>% pairs()

# tmin_cold and temp_sea are very strongly correlated, so remove temp_sea
f2 = vert.mean.ses ~ tmax_warm + tmin_cold + precip_sea + precip_wet + precip_dry
mod2.ecoregion = lme(f2 , random = ~1 | ecoregion, data = d, method = "REML")
summary(mod2.ecoregion)
# look at residuals
resid2.ecoregion = residuals(mod2.ecoregion, type = "normalized")
cor.test = data.frame(x = d$x, y = d$y, resid = resid2.ecoregion) %>% 
  slice_sample(n=1000)
cor2.ecoregion = ncf::correlog(x = cor.test$x, y = cor.test$y, z = cor.test$resid, increment = 100000, resamp = 99, latlon = F)
ncf:::plot.correlog(cor2.ecoregion, xlim = c(0,4e6))


# climate + veg
f3 = vert.mean.ses ~ tmax_warm + tmin_cold + precip_sea + precip_wet + precip_dry + canopy_height + veg_den
mod3.ecoregion = lme(f3 , random = ~1 | ecoregion, data = d, method = "REML")
summary(mod3.ecoregion)
# look at residuals
resid3.ecoregion = residuals(mod3.ecoregion, type = "normalized")
cor.test = data.frame(x = d$x, y = d$y, resid = resid3.ecoregion) %>% 
  slice_sample(n=1000)
cor3.ecoregion = ncf::correlog(x = cor.test$x, y = cor.test$y, z = cor.test$resid, increment = 100000, resamp = 99, latlon = F)
ncf:::plot.correlog(cor3.ecoregion, xlim = c(0,4e6), ylim = c(-0.5,0.5))


# H3: historical climate hypothesis


# THINNING POINTS  --------------------------------------------------------
#THIN POINTS 150000 KM APAPRT

# Try thinning points to 600km distance between points
dthin = dat %>% 
  filter(rich >= 5) %>% 
  filter(!is.na(vert.mean.ses)) %>% 
  drop_na() %>% 
  dplyr::select(x,y, rich) %>% 
  rast(crs = "+proj=cea +datum=WGS84")
dthin = spatSample(dthin, size = 10000, as.df = T, xy = T, values = T, na.rm = T, method = "random")
dthin.dist = dist(dthin)
min(dthin.dist)
dthin = terra::extract(datr, dthin[,c("x", "y")], xy = T)
dthin = dthin %>% 
  dplyr::select(!c(biome, ecoregion, realm)) %>% 
  left_join(key.biome, by = "key.biome") %>% 
  left_join(key.ecoregion, by = "key.ecoregion") %>% 
  left_join(key.realm, by = "key.realm")

d = dthin

# add biome as a random effect
mod.biome = lme(vert.mean.ses ~ canopy_height + veg_den, random = ~1 | biome, data = d, method = "REML")
summary(mod.biome)
# look at residuals
resid.biome = residuals(mod.biome, type = "normalized")
cor.test = data.frame(x = d$x, y = d$y, resid = resid.biome) %>% 
  slice_sample(n=1000)
cor.biome = ncf::correlog(x = cor.test$x, y = cor.test$y, z = cor.test$resid, increment = 100000, resamp = 99, latlon = F)
ncf:::plot.correlog(cor.biome, xlim = c(0,4e6))

# add biome as a random effect nested within realm
mod.biorealm = lme(vert.mean.ses ~ canopy_height + veg_den, random = ~1 | realm/biome, data = d, method = "REML")
summary(mod.biorealm)
# look at residuals
resid.biorealm = residuals(mod.biorealm, type = "normalized")
cor.test = data.frame(x = d$x, y = d$y, resid = resid.biorealm) %>% 
  slice_sample(n=1000)
cor.biorealm = ncf::correlog(x = cor.test$x, y = cor.test$y, z = cor.test$resid, increment = 100000, resamp = 99, latlon = F)
ncf:::plot.correlog(cor.biorealm, xlim = c(0,4e6))

# still spatial autocorrelation - try blocking by ecoregion
mod.ecoregion = lme(vert.mean.ses ~ canopy_height + veg_den, random = ~1 | ecoregion, data = d, method = "REML")
summary(mod.ecoregion)
# look at residuals
resid.ecoregion = residuals(mod.ecoregion, type = "normalized")
cor.test = data.frame(x = d$x, y = d$y, resid = resid.ecoregion) %>% 
  slice_sample(n=1000)
cor.ecoregion = ncf::correlog(x = cor.test$x, y = cor.test$y, z = cor.test$resid, increment = 100000, resamp = 99, latlon = F)
ncf:::plot.correlog(cor.ecoregion, xlim = c(0,4e6), ylim = c(-0.2,0.2))




# H2: climate hypothesis
d %>% dplyr::select(tmax_warm, tmin_cold, temp_sea, precip_wet, precip_dry, precip_sea) %>% pairs()

# tmin_cold and temp_sea are very strongly correlated, so remove temp_sea
f2 = vert.mean.ses ~ tmax_warm + tmin_cold + precip_sea + precip_wet + precip_dry
mod2.ecoregion = lme(f2 , random = ~1 | ecoregion, data = d, method = "REML")
summary(mod2.ecoregion)
# look at residuals
resid2.ecoregion = residuals(mod2.ecoregion, type = "normalized")
cor.test = data.frame(x = d$x, y = d$y, resid = resid2.ecoregion) %>% 
  slice_sample(n=1000)
cor2.ecoregion = ncf::correlog(x = cor.test$x, y = cor.test$y, z = cor.test$resid, increment = 100000, resamp = 99, latlon = F)
ncf:::plot.correlog(cor2.ecoregion, xlim = c(0,4e6), ylim = c(-0.2,0.2))


# climate + veg
f3 = vert.mean.ses ~ tmax_warm + tmin_cold + precip_sea + precip_wet + precip_dry + canopy_height + veg_den
mod3.ecoregion = lme(f3 , random = ~1 | ecoregion, data = d, method = "REML")
summary(mod3.ecoregion)
# look at residuals
resid3.ecoregion = residuals(mod3.ecoregion, type = "normalized")
cor.test = data.frame(x = d$x, y = d$y, resid = resid3.ecoregion) %>% 
  slice_sample(n=1000)
cor3.ecoregion = ncf::correlog(x = cor.test$x, y = cor.test$y, z = cor.test$resid, increment = 100000, resamp = 99, latlon = F)
ncf:::plot.correlog(cor3.ecoregion, xlim = c(0,4e6), ylim = c(-0.5,0.5))

test = d %>% group_by(ecoregion) %>% count()


# SDMTMB ------------------------------------------------------------------

# fit mixed effects model with spatial correlation structure
head(d)
d = d %>% 
  mutate_at(vars(canopy_height:veg_complexity), .funs = scale)
# convert xy coordinates into km
d$x = d$x/100000
d$y = d$y/100000
d$biome = factor(d$biome)
d$realm = factor(d$realm)
d$vert.mean.ses.scale = scale(d$vert.mean.ses)

mesh = make_mesh(d, c("x", "y"), cutoff = 3)
plot(mesh)

f1 = formula(vert.mean.ses ~ canopy_height + veg_den + (1|biome))
m.glm = sdmTMB(
  data = d,
  formula = f1,
  family = gaussian(),
  spatial = "off"
)
m.glm
sanity(m.glm)
resids.glm = data.frame(
  x = d$x, y = d$y, 
  resids = residuals(m.glm, type = "mle-mvn")) %>% 
  slice_sample(n = 1000)

ggplot(resids.glm, aes(x, y, col = resids)) +
  scale_color_gradient2() +
  geom_point() +
  coord_sf(crs = "+proj=cea +datum=WGS84")
cor.glm = ncf::correlog(x = resids.glm$x, y = resids.glm$y,
              z = resids.glm$resids, increment = 1,
              resamp = 99)
plot(cor.glm, xlim = c(0,50))

m.spatial = sdmTMB(
  data = d,
  formula = f1,
  mesh = mesh,
  spatial = "on",
  spatiotemporal = "off"
)
m.spatial
tidy(m.spatial, conf.int = T)
tidy(m.spatial, "ran_pars", conf.int = T)

d.resids = d
d.resids$resids = residuals(m.spatial, type = "mle-mvn")
qqnorm(d.resids$resids, ylim = c(-10,10))
qqline(d.resids$resids)
ggplot(d.resids, aes(x, y, col = resids)) +
  geom_point() +
  scale_color_gradient2()

preds = predict(m.spatial, newdata = d)

ggplot(preds, aes(x*10000, y*10000, color = est)) +
  geom_tile() +
  scale_color_continuous_divergingx("spectral") +
  coord_sf(crs = "+proj=cea +datum=WGS84") +
  ggtitle("fixed + random effects")

ggplot(preds, aes(x*10000, y*10000, color = est_non_rf)) +
  geom_tile() +
  scale_color_continuous_divergingx("spectral") +
  coord_sf(crs = "+proj=cea +datum=WGS84") +
  ggtitle("fixed effects only")

ggplot(preds, aes(x*10000, y*10000, color = omega_s)) +
  geom_tile() +
  scale_color_continuous_divergingx("spectral") +
  coord_sf(crs = "+proj=cea +datum=WGS84") +
  ggtitle("spatial random effects only")

ggplot(d.resids, aes(x*10000, y*10000, color = resids)) +
  geom_tile() +
  scale_color_continuous_divergingx("spectral") +
  coord_sf(crs = "+proj=cea +datum=WGS84") +
  ggtitle("fixed + random effects")

d.resids.test = d.resids %>% slice_sample(n=1000)
cor.spatial = ncf::correlog(x = d.resids.test$x,
                            y = d.resids.test$y,
                            z = d.resids.test$resids,
                            increment = 10,
                            resamp = 99)
plot(cor.spatial)

# look at conditional effects
nd = data.frame(canopy_height = seq(min(d.resids$canopy_height),
                            max(d.resids$canopy_height),
                            length.out = 100),
                veg_den = mean(d.resids$veg_den),
                biome = "Tropical & Subtropical Moist Broadleaf Forests")
nd$biome = factor(nd$biome)
p = predict(m.spatial, newdata = nd, se_fit = T, re_form = NA)

ggplot(p, aes(canopy_height, est,
              ymin = est - 1.96*est_se,
              ymax = est + 1.93*est_se)) +
  geom_line() +
  geom_ribbon(alpha = 0.4) +
  scale_x_continuous() +
  coord_cartesian(expand = F) +
  labs(x = "Canopy Height", y = "SES Mean Verticality")
