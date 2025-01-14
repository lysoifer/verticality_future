# Plot coefficients for sdmtmb models using 50 km resolution data

library(tidyverse)
library(lemon)
library(ggh4x)
library(spatialreg)
library(sdmTMB)

# VERTMEAN ----------------------------------------------------------------

# ADD variable importance weights to the graphs ##

load("results/sdmTMB_models/amphibians_meanvert.RData")
amph.mod = mod

load("results/sdmTMB_models/mammals_meanvert.RData")
mammals.mod = mod

load("results/sdmTMB_models/reptiles_meanvert.RData")
rept.mod = mod

load("results/sdmTMB_models/birds_meanvert.RData")
birds.mod = mod

# * - get coefficients with 95% CI for each class -----------------------------------------
amph.coefs = tidy(amph.mod, conf.int = T, conf.level = 0.95, effects = "fixed")
amph.coefs$class = "Amphibians"

mammals.coefs = tidy(mammals.mod, conf.int = T, conf.level = 0.95, effects = "fixed")
mammals.coefs$class = "Mammals"

reptiles.coefs = tidy(rept.mod, conf.int = T, conf.level = 0.95, effects = "fixed")
reptiles.coefs$class = "Reptiles"

birds.coefs = tidy(birds.mod, conf.int = T, conf.level = 0.95, effects = "fixed")
birds.coefs$class = "Birds"


meanvert.coefs = bind_rows(amph.coefs, mammals.coefs, reptiles.coefs, birds.coefs)%>% 
  mutate(class = factor(class, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")),
         type = "Mean Verticality",
         sig = ifelse(conf.low*conf.high > 0, "sig", "notsig"))


# * - plot --------------------------------------------------------------------

meanvert.coefs %>% 
  filter(term != "(Intercept)") %>% 
  ggplot() +
  geom_pointrange(aes(x = estimate, y = reorder(term, estimate), xmin = conf.low, xmax = conf.high)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  #geom_text(data = meanvert.r2, aes(label = paste0("R\u00b2 = ", r2)), x = 0.042, y = 12, inherit.aes = F, hjust = 1, vjust = 0) +
  #annotate(geom = "text", x = -0.005, y = 7, label = paste0("R\u00b2 = ", r2)) +
  #scale_color_manual(values = c("black", "red3")) +
  #scale_y_discrete("") +
  #coord_cartesian(clip = "off", ylim = c(0,13)) +
  facet_wrap(facets = ~class) +
  ggtitle("Mean Verticality") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(color = "black", fill = "white"),
        panel.grid = element_line(color = "gray90"),
        strip.background = element_rect(color = "black"))

# SES VERTMEAN ----------------------------------------------------------------

load("results/sdmTMB_models/amphibians_sesvert.RData")
amph.mod = mod

load("results/sdmTMB_models/mammals_sesvert.RData")
mammals.mod = mod

load("results/sdmTMB_models/reptiles_sesvert.RData")
rept.mod = mod

load("results/sdmTMB_models/birds_sesvert.RData")
birds.mod = mod



# * - get coefficients for each class -----------------------------------------
amph.coefs = tidy(amph.mod, conf.int = T, conf.level = 0.95, effects = "fixed")
amph.coefs$class = "Amphibians"

mammals.coefs = tidy(mammals.mod, conf.int = T, conf.level = 0.95, effects = "fixed")
mammals.coefs$class = "Mammals"

reptiles.coefs = tidy(rept.mod, conf.int = T, conf.level = 0.95, effects = "fixed")
reptiles.coefs$class = "Reptiles"

birds.coefs = tidy(birds.mod, conf.int = T, conf.level = 0.95, effects = "fixed")
birds.coefs$class = "Birds"


sesvert.coefs = bind_rows(amph.coefs, mammals.coefs, reptiles.coefs, birds.coefs)%>% 
  mutate(class = factor(class, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")),
         type = "SES Mean Verticality",
         sig = ifelse(conf.low*conf.high > 0, "sig", "notsig"))



# * - plot --------------------------------------------------------------------

sesvert.coefs %>% 
  filter(term != "(Intercept)") %>%
  ggplot(aes(x = estimate, y = reorder(term, estimate), xmin = conf.low, xmax = conf.high)) +
  geom_pointrange() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_discrete("") +
  facet_wrap(facets = ~class) +
  ggtitle("SES Mean Verticality") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(color = "black", fill = "white"),
        panel.grid = element_line(color = "gray90"),
        strip.background = element_rect(color = "black"))


# plot meanvert and sesmeanvert -------------------------------------------

coefs = bind_rows(sesvert.coefs, meanvert.coefs) %>% 
  filter(term != "(Intercept)") %>%
  mutate(term = factor(term, levels = c(
    "log_clim_velocity", "canopy_height", "veg_den",
    "precip_wet", "log_precip_dry",
    "I(tmax_warm^2)", "tmax_warm", "tmin_cold")),
    class = factor(class, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")))

coefs = bind_rows(sesvert.coefs, meanvert.coefs) %>% 
  filter(term != "(Intercept)") %>%
  mutate(term = case_when(term == "log_clim_velocity" ~ "Climate velocity",
                          term == "canopy_height" ~ "Canopy height",
                          term == "veg_den" ~ "Vegetation density",
                          term == "precip_wet" ~ "Wet season precip",
                          term == "log_precip_dry" ~ "log(dry season precip)",
                          term == "I(tmax_warm^2)" ~ "(Max temp wamest month)\u00b2",
                          term == "tmax_warm" ~ "Max temp warmest month",
                          term == "tmin_cold" ~ "Min temp coldest month"),
         term = factor(term, levels = c(
    "Climate velocity", "Canopy height", "Vegetation density",
    "Wet season precip", "log(dry season precip)",
    "(Max temp wamest month)\u00b2", "Max temp warmest month", "Min temp coldest month")),
    class = factor(class, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")))


library(ggh4x)
coef.plt = coefs %>%
  ggplot(aes(x = estimate, y = term, xmin = conf.low, xmax = conf.high, color = sig)) +
  geom_pointrange(size = 0.15) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  #geom_text(data = r2, aes(label = paste0("R\u00b2 = ", r2), x = Estimate), y = 2, inherit.aes = F, hjust = 1, vjust = 0) +
  #geom_text(data = r2, aes(label = paste0("R\u00b2 trend = ", r2.trend), x = Estimate), y = 0.9, inherit.aes = F, hjust = 1, vjust = 0) +
  #scale_color_manual(values = c("black", "red3")) +
  #scale_size_continuous(range = c(0.02,0.5)) +
  scale_y_discrete("") +
  scale_x_continuous("Estimate") +
  scale_color_manual(values = c("black", "red3")) +
  facet_grid2(rows = vars(type), cols = vars(class), scales = "free_x", independent = "x") +
  scale_x_facet(ROW == 1, limits = c(-0.08, 0.25)) +
  scale_x_facet(ROW == 2, limits = c(-0.2, 0.8)) +
  #facet_grid(class ~ type, scales = "free") +
  #ggtitle("SES Mean Verticality") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(color = "black", fill = "white"),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(color = "black"),
        legend.position = "none",
        legend.box.margin = margin(0,0,0,0),
        legend.box.spacing = unit(0.1, "mm"),
        legend.title.position = "top",
        legend.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 6))

png("figures/main_figs/fig2_sdmTMB_coefs.png", width = 180, height = 98, res = 600, units = "mm")
coef.plt
dev.off()
