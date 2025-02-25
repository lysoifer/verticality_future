library(sdmTMB)
library(ggeffects)
library(foreach)
library(doParallel)

amph = readRDS("results/sdmTMB_models2/predictions/amphibians_sesvert.rds")
amphmod = amph$mod

birds = readRDS("results/sdmTMB_models2/predictions/birds_sesvert.RData")
birdsmod = birds$mod

mammals = readRDS("results/sdmTMB_models2/predictions/mammals_sesvert.rds")
mammalsmod = mammals$mod

repts = readRDS("results/sdmTMB_models2/predictions/reptiles_sesvert.RData")
reptsmod = repts$mod

range(birdsmod$data$log_precip_dry)
range(birdsmod$data$tmax_warm)
range(birdsmod$data$tmin_cold)
range(birdsmod$data$canopy_height)
range(birdsmod$data$veg_den)
range(birdsmod$data$precip_wet)
range(birdsmod$data$precip_warm)
range(birdsmod$data$log_clim_velocity)

var = list(
  "log_precip_dry [-4.6:2.3 by=0.1]",
  "tmax_warm [-3.6:2.7 by=0.1]",
  "tmin_cold [-2.54:1.37 by=0.1]",
  "canopy_height [-1.3:3.5 by=0.1]",
  "veg_den [-1.3:9.0 by=0.1]",
  "precip_wet [-1.4:11.5 by=0.1]",
  "precip_warm [-1.3:15.4 by=0.1]",
  "log_clim_velocity [-4.3:2.5 by=0.1]"
  )

amphplots = foreach(v = 1:length(var)) %do% {
  g = predict_response(amphmod, var[[v]])
  g = cbind(g)
  g$term = var[[v]]
  g
}

saveRDS(amphplots, "results/sdmTMB_models2/marginal_effects/amph.rds")

# mammals

mammalsplots = foreach(v = 1:length(var)) %do% {
  g = predict_response(mammalsmod, var[[v]])
  g = cbind(g)
  g$term = var[[v]]
  g$taxa = "mammals"
  g
}

saveRDS(mammalsplots, "results/sdmTMB_models2/marginal_effects/mammals.rds")

# reptiles

reptsplots = foreach(v = 1:length(var)) %do% {
  g = predict_response(reptsmod, var[[v]])
  g = cbind(g)
  g$term = var[[v]]
  g$taxa = "reptiles"
  g
}


saveRDS(reptsplots, "results/sdmTMB_models2/marginal_effects/reptiles.rds")


# birds
birdsplots = foreach(v = 1:length(var)) %do% {
  g = predict_response(birdsmod, var[[v]])
  g = cbind(g)
  g$term = var[[v]]
  g$taxa = "birds"
  g
}

saveRDS(birdsplots, "results/sdmTMB_models2/marginal_effects/birds.rds")


# plot conditional effects
ampheffects = readRDS("results/sdmTMB_models2/marginal_effects/amph.rds")
ampheffects = rbindlist(ampheffects) %>% 
  mutate(taxa = "Amphibians")

birdeffects = readRDS("results/sdmTMB_models2/marginal_effects/birds.rds")
birdeffects = rbindlist(birdeffects) %>% 
  mutate(taxa = "Birds")

repteffects = readRDS("results/sdmTMB_models2/marginal_effects/reptiles.rds")
repteffects = rbindlist(repteffects) %>% 
  mutate(taxa = "Reptiles")

mammaleffects = readRDS("results/sdmTMB_models2/marginal_effects/mammals.rds")
mammaleffects = rbindlist(mammaleffects) %>% 
  mutate(taxa = "Mammals")

df = bind_rows(ampheffects, birdeffects, repteffects, mammaleffects) %>% 
  separate(term, into = c("term", "by1", "by2"), sep = " ") %>% 
  mutate(term = case_when(term == "canopy_height" ~ "Canopy\nheight",
                          term == "log_clim_velocity" ~ "log(Climate\nvelocity)",
                          term == "log_precip_dry" ~ "log(Precip\ndry)",
                          term == "precip_warm" ~ "Precip\n warm",
                          term == "precip_wet" ~ "Precip\n wet",
                          term == "tmax_warm" ~ "Tmax warm",
                          term == "tmin_cold" ~ "Tmin cold",
                          term == "veg_den" ~ "Vegetation\ndensity"),
         taxa = factor(taxa, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")))

ggplot(df, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  facet_grid(taxa ~ term, switch = "x") +
  scale_y_continuous("Predicted") +
  theme_bw() +
  theme(strip.placement = "outside",
        strip.background.x = element_blank(),
        axis.title.x = element_blank())

ggsave(filename = "figures/supp_figs/conditional_effects.png", height = 120, width = 180, units = "mm", dpi = 300)


# canopy interaction term -------------------------------------------------

var = list(
  c("canopy_height [-4:4 by=0.1]", "log_precip_dry [-4:4 by=2]"),
  c("canopy_height [-4:4 by=0.1]", "tmax_warm [-4:4 by=2]"),
  c("canopy_height [-4:4 by=0.1]", "tmin_cold[-4:4 by=2]"),
  c("canopy_height [-4:4 by=0.1]", "precip_warm[-4:4 by=2]")
)

term = c(
  "canopy-log_precip_dry",
  "canopy-tmax_warm",
  "canopy-tmin_cold",
  "canopy-precip_warm"
)

amphplots = list()
for(v in 1:3) {
  g = predict_response(amphmod, var[[v]])
  g = cbind(g)
  g$term = term[v]
  g$taxa = "Amphibians"
  amphplots[[v]] = g
  saveRDS(amphplots, "results/sdmTMB_models2/marginal_effects/canopy_interactions/amph_interactions.rds")
}

mammalsplots = list()
for(v in 1:length(var)) {
  g = predict_response(mammalsmod, var[[v]])
  g = cbind(g)
  g$term = term[v]
  g$taxa = "Mammals"
  mammalsplots[[v]] = g
  saveRDS(mammalsplots, "results/sdmTMB_models2/marginal_effects/canopy_interactions/mammals_interactions.rds")
}

birdsplots = list()
for(v in 1:length(var)) {
  g = predict_response(birdsmod, var[[v]])
  g = cbind(g)
  g$term = term[v]
  g$taxa = "Birds"
  birdsplots[[v]] = g
  saveRDS(birdsplots, "results/sdmTMB_models2/marginal_effects/canopy_interactions/birds_interactions.rds")
}

reptsplots = list()
for(v in 1:length(var)) {
  g = predict_response(reptsmod, var[[v]])
  g = cbind(g)
  g$term = term[v]
  g$taxa = "Reptiles"
  reptsplots[[v]] = g
  saveRDS(reptsplots, "results/sdmTMB_models2/marginal_effects/canopy_interactions/repts_interactions.rds")
}

amph.df = readRDS("results/sdmTMB_models2/marginal_effects/canopy_interactions/amph_interactions.rds")
amph.df = rbindlist(amph.df) %>% mutate(taxa = "Amphibians")

rept.df = readRDS("results/sdmTMB_models2/marginal_effects/canopy_interactions/repts_interactions.rds")
rept.df = rbindlist(rept.df) %>% mutate(taxa = "Reptiles")

mammals.df = readRDS("results/sdmTMB_models2/marginal_effects/canopy_interactions/mammals_interactions.rds")
mammals.df = rbindlist(amph.df) %>% mutate(taxa = "Mammals")

birds.df = readRDS("results/sdmTMB_models2/marginal_effects/canopy_interactions/birds_interactions.rds")
birds.df = rbindlist(amph.df) %>% mutate(taxa = "Birds")

df = bind_rows(amph.df, rept.df, mammals.df, birds.df) %>% 
  filter(x >= -1.3 & x <= 3.5) %>% 
  filter(group == -2 | group == 0 | group == 2) %>% 
  mutate(term = case_when(term == "canopy-log_precip_dry" ~ "log(Precip dry)",
                          term == "canopy-tmax_warm" ~ "Tmax warm",
                          term == "canopy-tmin_cold" ~ "Tmin cold",
                          term ==  "canopy-precip_warm" ~ "Precip warm"))

ggplot(df, aes(x, predicted, color = group, fill = group)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  facet_grid(term ~ taxa) +
  scale_y_continuous("Predicted") +
  scale_x_continuous("Scaled canopy height") +
  scale_color_discrete(guide = guide_legend("Scaled\ninteraction\nterm")) +
  scale_fill_discrete(guide = guide_legend("Scaled\ninteraction\nterm")) +
  theme_bw()

ggsave(filename = "figures/supp_figs/conditional_effects_interactions.png",
       height = 120, width = 150, units = "mm", dpi = 300)

min = amphmod$data %>% 
  dplyr::select(canopy_height, precip_wet, log_precip_dry, precip_warm, tmax_warm, tmin_cold, veg_den, log_clim_velocity) %>% 
  summarise(across(canopy_height:log_clim_velocity, ~ min(.x, na.rm = TRUE)))

max = amphmod$data %>% 
  dplyr::select(canopy_height, precip_wet, log_precip_dry, precip_warm, tmax_warm, tmin_cold, veg_den, log_clim_velocity) %>% 
  summarise(across(canopy_height:log_clim_velocity, ~ max(.x, na.rm = TRUE)))

df = as.data.frame(matrix(c(min, max), ncol = 2))
df$var = names(min)

amph.df[[4]] %>% 
  ggplot(aes(x = x, y = predicted)) +
  geom_line() +
  theme_bw()


amph.df[[1]] %>% 
  ggplot(aes(x = x, y = predicted, color = group)) +
  geom_line() +
  theme_bw()

