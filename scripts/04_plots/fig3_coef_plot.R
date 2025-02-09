# Plot coefficients for sdmtmb models using 50 km resolution data

library(tidyverse)
library(lemon)
library(ggh4x)
library(spatialreg)
library(sdmTMB)


# SES Verticality Model Coefficients --------------------------------------

amph = readRDS("results/sdmTMB_models2/amphibians_sesvert.rds")
birds = readRDS("results/sdmTMB_models2/birds_sesvert.rds")
mammals = readRDS("results/sdmTMB_models2/mammals_sesvert.rds")
repts = readRDS("results/sdmTMB_models2/reptiles_sesvert.rds")

mods = list(amph$bestmod, birds$bestmod, mammals$bestmod, repts$bestmod)
mods = lapply(mods, tidy, effects = "fixed", conf.int = T)
mods[[1]]$taxa = "Amphibians"
mods[[2]]$taxa = "Birds"
mods[[3]]$taxa = "Mammals"
mods[[4]]$taxa = "Reptiles"

mods = bind_rows(mods)
mods = mods %>% 
  filter(term != "(Intercept)") %>%
  mutate(taxa = factor(taxa, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")),
         term = case_when(term == "canopy_height:tmin_cold" ~ "canopy height:Tmin[cold]",
                          term == "canopy_height:tmax_warm" ~ "canopy height:Tmax[warm]",
                          term == "canopy_height:log_precip_dry" ~ "canopy height:log(Precip[dry]",
                          term == "precip_warm:canopy_height" ~ "Precip[warm]:canopy height",
                          term == "log_clim_velocity" ~ "log(climate velocity)", 
                          term == "veg_den" ~ "Vegetation density",
                          term == "canopy_height" ~ "Canopy height",
                          term == "I(canopy_height^2)" ~ "(Canopy height)\u00b2",
                          term == "I(precip_warm^2)" ~ "(Precip[warm])\u00b2",
                          term == "precip_warm" ~ "Precip[warm]",
                          term == "precip_wet" ~ "Precip[wet]",
                          term == "I(log_precip_dry^2)" ~ "(log(Precip[dry]))\u00b2",
                          term == "log_precip_dry" ~ "log(Precip[dry])",
                          term == "I(tmax_warm^2)" ~ "(Tmax[warm])\u00b2",
                          term == "tmax_warm" ~ "Tmax[warm]",
                          term == "I(tmin_cold^2)" ~ "(Tmin[cold])\u00b2",
                          term == "tmin_cold" ~ "Tmin[cold]"),
        term = factor(term, levels = c(
          "canopy height:Tmin[cold]",
          "canopy height:Tmax[warm]",
          "canopy height:log(Precip[dry]",
          "Precip[warm]:canopy height",
          "log(climate velocity)",
          "Vegetation density",
          "(Canopy height)\u00b2",
          "Canopy height",
          "(Precip[warm])\u00b2",
          "Precip[warm]",
          "Precip[wet]",
          "(log(Precip[dry]))\u00b2",
          "log(Precip[dry])",
          "(Tmax[warm])\u00b2",
          "Tmax[warm]",
          "(Tmin[cold])\u00b2",
          "Tmin[cold]")),
        sig = ifelse(conf.low*conf.high > 0, "sig", "notsig"))




library(ggh4x)
coef.plt = mods %>%
  ggplot(aes(x = estimate, y = term, xmin = conf.low, xmax = conf.high, color = sig)) +
  geom_pointrange(size = 0.15) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_discrete("",
                   labels = c("canopy height:Tmin[cold]" = expression("Canopy height:Tmin"[cold]),
                              "canopy height:Tmax[warm]" = expression("Canopy height:Tmax"[warm]),
                              "canopy height:log(Precip[dry]" = expression("Canopy height:log(Precip"["dry"]*")"),
                              "Precip[warm]:canopy height" = expression("Precip"[warm]*":Canopy height"),
                              "(Precip[warm])\u00b2" = expression("(Precip"["warm"]*")\u00b2"),
                              "Precip[warm]" = expression("Precip"[warm]),
                              "Precip[wet]" = expression("Precip"[wet]),
                              "(log(Precip[dry]))\u00b2" = expression("(log(Precip"["dry"]*"))\u00b2"),
                              "log(Precip[dry])" = expression("log(Precip"["dry"]*")"),
                              "(Tmax[warm])\u00b2" = expression("(Tmax"["warm"]*")\u00b2"),
                              "Tmax[warm]" = expression("Tmax"[warm]),
                              "(Tmin[cold])\u00b2" = expression("(Tmin"["cold"]*")\u00b2"),
                              "Tmin[cold]" = expression("Tmin"[cold]))) +
  #scale_x_continuous("Effect size of standardized predictors") +
  scale_color_manual(values = c("black", "red3")) +
  facet_wrap(vars(taxa)) +
  # facet_grid2(rows = vars(type), cols = vars(class), scales = "free_x", independent = "x") +
  # scale_x_facet(ROW == 1, limits = c(-0.08, 0.18)) +
  # scale_x_facet(ROW == 2, limits = c(-0.2, 0.7)) +
  labs(x = "Effect size of standardized predictors") +
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

png("figures/main_figs/fig3/fig3_sdmTMB_coefs.png", width = 180, height = 180, res = 600, units = "mm")
coef.plt
dev.off()