library(tidyverse)
library(lemon)
library(ggh4x)
library(spatialreg)

# VERTMEAN 50km resolution sdmTMB models ----------------------------------------------------------------

# ADD variable importance weights to the graphs ##

load("results/sdmTMB_models/amphibians_meanvert.RData")
amph.mod = m.final

load("results/sdmTMB_models/mammals_meanvert.RData")
mammals.mod = m.final

load("results/sdmTMB_models/reptiles_meanvert.RData")
rept.mod = m.final

load("results/sdmTMB_models/birds_meanvert.RData")
birds.mod = m.final

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
         type = "Mean Verticality")


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
amph.mod = m.final

load("results/sdmTMB_models/mammals_sesvert.RData")
mammals.mod = m.final

load("results/sdmTMB_models/reptiles_sesvert.RData")
rept.mod = m.final

load("results/sdmTMB_models/birds_sesvert.RData")
birds.mod = m.final



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
         type = "SES Mean Verticality")



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

library(ggh4x)
coef.plt = coefs %>%
  ggplot(aes(x = estimate, y = term, xmin = conf.low, xmax = conf.high)) +
  geom_pointrange(size = 0.15) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  #geom_text(data = r2, aes(label = paste0("R\u00b2 = ", r2), x = Estimate), y = 2, inherit.aes = F, hjust = 1, vjust = 0) +
  #geom_text(data = r2, aes(label = paste0("R\u00b2 trend = ", r2.trend), x = Estimate), y = 0.9, inherit.aes = F, hjust = 1, vjust = 0) +
  #scale_color_manual(values = c("black", "red3")) +
  #scale_size_continuous(range = c(0.02,0.5)) +
  scale_y_discrete("") +
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
        legend.position = "bottom",
        legend.box.margin = margin(0,0,0,0),
        legend.box.spacing = unit(0.1, "mm"),
        legend.title.position = "top",
        legend.title = element_text(hjust = 0.5))

png("figures/main_figs/sdmTMB_coefs.png", width = 220, height = 120, res = 300, units = "mm")
coef.plt
dev.off()


# FOREST ONLY COEF PLOTS (first iteration without forest based SES and missing some wooded habitats) --------------------------------------------------

# * - VERTMEAN ----------------------------------------------------------------

load("results/sar_mods_forest_only/amphibians/amphibians_sar_vertmean2.RData")

amph.mods = best_mods_pred
amph.avgmod = vert.avg

load("results/sar_mods_forest_only/mammals/mammals_sar_vertmean2.RData")
mammals.mods = best_mods_pred
mammals.avgmod = vert.avg

load("results/sar_mods_forest_only/reptiles/reptiles_sar_vertmean2.RData")
rept.mods = best_mods_pred
rept.avgmod = vert.avg

load("results/sar_mods_forest_only/birds/birdselton_sar_vertmean.RData")
birds.mods = best_mods_pred
birds.avgmod = vert.full


# * - get coefficients for each class -----------------------------------------
amph.sum = summary(amph.avgmod)
amph.coefs = amph.sum$coefmat.full %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "var") %>% 
  rename(SE = "Std. Error", z_value = "z value", p_value = "Pr(>|z|)") %>% 
  mutate(class = "Amphibians")

mammals.sum = summary(mammals.avgmod)
mammals.coefs = mammals.sum$coefmat.full %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "var") %>% 
  rename(SE = "Std. Error", z_value = "z value", p_value = "Pr(>|z|)") %>% 
  mutate(class = "Mammals")

reptiles.sum = summary(rept.avgmod)
repts.coefs = reptiles.sum$coefmat.full %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "var") %>% 
  rename(SE = "Std. Error", z_value = "z value", p_value = "Pr(>|z|)") %>% 
  mutate(class = "Reptiles")

birds.sum = summary(birds.avgmod)
birds.coefs = birds.sum$Coef %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "var") %>% 
  rename(SE = "Std. Error", z_value = "z value", p_value = "Pr(>|z|)") %>% 
  mutate(class = "Birds")

meanvert.coefs = bind_rows(amph.coefs, mammals.coefs, repts.coefs, birds.coefs)%>% 
  mutate(class = factor(class, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")),
         type = "Mean Verticality")


# * - make r2 table -----------------------------------------------------------

amph.r2 = mean(sapply(amph.mods, "[[", 4))
mammals.r2 = mean(sapply(mammals.mods, "[[", 5))
repts.r2 = mean(sapply(rept.mods, "[[", 4))
birds.r2 = mean(sapply(birds.mods, "[[", 4))

meanvert.r2 = data.frame(class = c("Amphibians", "Mammals", "Reptiles", "Birds"),
                         r2 = c(amph.r2, mammals.r2, repts.r2, round(birds.r2,2)))%>% 
  mutate(class = factor(class, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")),
         type = "Mean Verticality")



# * - plot --------------------------------------------------------------------

meanvert.coefs %>% 
  filter(var != "(Intercept)") %>% 
  ggplot(aes(x = Estimate, y = reorder(var, Estimate), xmin = Estimate-SE*1.96, xmax = Estimate+SE*1.96, color = p_value<0.05)) +
  geom_pointrange() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text(data = meanvert.r2, aes(label = paste0("R\u00b2 = ", r2)), x = 0.042, y = 1.2, inherit.aes = F, hjust = 1, vjust = 0) +
  #annotate(geom = "text", x = -0.005, y = 7, label = paste0("R\u00b2 = ", r2)) +
  scale_color_manual(values = c("black", "red3")) +
  scale_y_discrete("") +
  facet_wrap(facets = ~class) +
  ggtitle("Mean Verticality") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(color = "black", fill = "white"),
        panel.grid = element_line(color = "gray90"),
        strip.background = element_rect(color = "black"))


# SES VERTMEAN ----------------------------------------------------------------



# * - get coefficients for each class -----------------------------------------

amph.coefs = amph.sum$coefmat.full %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "var") %>% 
  rename(SE = "Std. Error", z_value = "z value", p_value = "Pr(>|z|)") %>% 
  mutate(class = "Amphibians")

mammals.sum = summary(mammals.avgmod)
mammals.coefs = mammals.sum$coefmat.full %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "var") %>% 
  rename(SE = "Std. Error", z_value = "z value", p_value = "Pr(>|z|)") %>% 
  mutate(class = "Mammals")

reptiles.sum = summary(rept.avgmod)
repts.coefs = reptiles.sum$coefmat.full %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "var") %>% 
  rename(SE = "Std. Error", z_value = "z value", p_value = "Pr(>|z|)") %>% 
  mutate(class = "Reptiles")

birds.sum = summary(birds.avgmod)
birds.coefs = birds.sum$Coef %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "var") %>% 
  rename(SE = "Std. Error", z_value = "z value", p_value = "Pr(>|z|)") %>% 
  mutate(class = "Birds")

sesvert.coefs = bind_rows(amph.coefs, mammals.coefs, repts.coefs, birds.coefs) %>% 
  mutate(class = factor(class, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")),
         type = "SES Mean Verticality")


# * - make r2 table -----------------------------------------------------------

amph.r2 = mean(sapply(amph.mods, "[[", 4))
mammals.r2 = mean(sapply(mammals.mods, "[[", 4))
repts.r2 = mean(sapply(rept.mods, "[[", 4))
birds.r2 = mean(sapply(birds.mods, "[[", 4))

sesvert.r2 = data.frame(class = c("Amphibians", "Mammals", "Reptiles", "Birds"),
                        r2 = c(amph.r2, mammals.r2, repts.r2, round(birds.r2,2))) %>% 
  mutate(class = factor(class, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")),
         type = "SES Mean Verticality")




# * - plot --------------------------------------------------------------------

sesvert.coefs %>% 
  filter(var != "(Intercept)") %>%
  ggplot(aes(x = Estimate, y = reorder(var, Estimate), xmin = Estimate-SE*1.96, xmax = Estimate+SE*1.96, color = p_value<0.05)) +
  geom_pointrange() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text(data = sesvert.r2, aes(label = paste0("R\u00b2 = ", r2)), x = 0.6, y = 1.2, inherit.aes = F, hjust = 1, vjust = 0) +
  scale_color_manual(values = c("black", "red3")) +
  scale_y_discrete("") +
  facet_wrap(facets = ~class) +
  ggtitle("SES Mean Verticality") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(color = "black", fill = "white"),
        panel.grid = element_line(color = "gray90"),
        strip.background = element_rect(color = "black"))


# plot meanvert and sesmeanvert -------------------------------------------

coefs = bind_rows(sesvert.coefs, meanvert.coefs) %>% 
  filter(var != "(Intercept)") %>%
  mutate(var = factor(var, levels = c(
    "log_elev", "log_clim_velocity", "veg_complexity",
    "precip_wet", "log_precip_wet", "precip_dry", "log_precip_dry",
    "precip_sea", "log_precip_sea",
    "tmax_warm", "tmin_cold")))
r2 = bind_rows(sesvert.r2, meanvert.r2) %>% 
  mutate(Estimate = ifelse(type == "SES Mean Verticality", 0.75, 0.06))

coef.plt = coefs  %>% 
  ggplot(aes(x = Estimate, y = var, xmin = Estimate-SE*1.96, xmax = Estimate+SE*1.96, color = p_value<0.05)) +
  geom_pointrange() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text(data = r2, aes(label = paste0("R\u00b2 = ", r2), x = Estimate), y = 0.9, inherit.aes = F, hjust = 1, vjust = 0) +
  scale_color_manual(values = c("black", "red3")) +
  scale_y_discrete("") +
  facet_grid2(rows = vars(type), cols = vars(class), scales = "free_x", independent = "x") +
  scale_x_facet(ROW == 1, limits = c(-0.06, 0.06)) +
  scale_x_facet(ROW == 2, limits = c(-0.75, 0.75)) +
  #facet_grid(class ~ type, scales = "free") +
  #ggtitle("SES Mean Verticality") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(color = "black", fill = "white"),
        panel.grid = element_line(color = "gray90"),
        strip.background = element_rect(color = "black"),
        legend.position = "bottom",
        legend.box.margin = margin(0,0,0,0),
        legend.box.spacing = unit(0.1, "mm"))

png("figures/sar_coef_plots/sar_mods_forest_only.png", width = 220, height = 120, res = 300, units = "mm")
coef.plt
dev.off()

# VERTMEAN ----------------------------------------------------------------

# ADD variable importance weights to the graphs ##

load("results/sar_mods_forestOnly_forestSES/amphibians/amphibians_sar_meanvert.RData")

#amph.mods = best_mods_pred
amph.avgmod = vert.avg
amph.r2 = r2
amph.r2.trend = r2.trend

load("results/sar_mods_forestOnly_forestSES/mammals/mammals_sar_meanvert.RData")
#mammals.mods = best_mods_pred
mammals.avgmod = vert.avg
mammals.r2 = r2
mammals.r2.trend = r2.trend

load("results/sar_mods_forestOnly_forestSES/reptiles/reptiles_sar_meanvert.RData")
#rept.mods = best_mods_pred
rept.avgmod = vert.avg
rept.r2 = r2
rept.r2.trend = r2.trend

load("results/sar_mods_forestOnly_forestSES/birds/birds_sar_meanvert.RData")
#birds.mods = best_mods_pred
birds.avgmod = vert.avg
birds.r2 = r2
birds.r2.trend = r2.trend

# * - get coefficients for each class -----------------------------------------
amph.sum = summary(amph.avgmod)
amph.coefs = amph.sum$coefmat.full %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "var") %>% 
  rename(SE = "Std. Error", z_value = "z value", p_value = "Pr(>|z|)") %>% 
  mutate(class = "Amphibians")

mammals.sum = summary(mammals.avgmod)
mammals.coefs = mammals.sum$coefmat.full %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "var") %>% 
  rename(SE = "Std. Error", z_value = "z value", p_value = "Pr(>|z|)") %>% 
  mutate(class = "Mammals")

reptiles.sum = summary(rept.avgmod)
repts.coefs = reptiles.sum$coefmat.full %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "var") %>% 
  rename(SE = "Std. Error", z_value = "z value", p_value = "Pr(>|z|)") %>% 
  mutate(class = "Reptiles")

birds.sum = summary(birds.avgmod)
birds.coefs = birds.sum$coefmat.full %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "var") %>% 
  rename(SE = "Std. Error", z_value = "z value", p_value = "Pr(>|z|)") %>% 
  mutate(class = "Birds")

meanvert.coefs = bind_rows(amph.coefs, mammals.coefs, repts.coefs, birds.coefs)%>% 
  mutate(class = factor(class, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")),
         type = "Mean Verticality")

# add variable importance weights
varimp.amph = data.frame(varimp = amph.avgmod$sw, class = "Amphibians") %>% 
  rownames_to_column("var")
varimp.mammals = data.frame(varimp = mammals.avgmod$sw, class = "Mammals") %>% 
  rownames_to_column("var")
varimp.rept = data.frame(varimp = rept.avgmod$sw, class = "Reptiles") %>% 
  rownames_to_column("var")
varimp.birds = data.frame(varimp = birds.avgmod$sw, class = "Birds") %>% 
  rownames_to_column("var")
varimp = rbind(varimp.amph, varimp.mammals, varimp.rept, varimp.birds)

meanvert.coefs = left_join(meanvert.coefs, varimp, by = c("class", "var"))
meanvert.coefs$varimp = as.numeric(meanvert.coefs$varimp)
# * - make r2 table -----------------------------------------------------------

# amph.r2 = mean(sapply(amph.mods, "[[", 4))
# mammals.r2 = mean(sapply(mammals.mods, "[[", 5))
# repts.r2 = mean(sapply(rept.mods, "[[", 4))
# birds.r2 = mean(sapply(birds.mods, "[[", 4))

meanvert.r2 = data.frame(class = c("Amphibians", "Mammals", "Reptiles", "Birds"),
                r2 = c(amph.r2, mammals.r2, rept.r2, birds.r2),
                r2.trend = c(amph.r2.trend, mammals.r2.trend, rept.r2.trend, birds.r2.trend))%>% 
  mutate(class = factor(class, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")),
         type = "Mean Verticality",
         r2 = round(r2, 2),
         r2.trend = round(r2.trend, 2))



# * - plot --------------------------------------------------------------------

meanvert.coefs %>% 
  filter(var != "(Intercept)") %>% 
  ggplot() +
  geom_pointrange(aes(x = Estimate, y = reorder(var, Estimate), xmin = Estimate-SE*1.96, xmax = Estimate+SE*1.96, 
                      color = p_value<0.05, size = varimp)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text(data = meanvert.r2, aes(label = paste0("R\u00b2 = ", r2)), x = 0.042, y = 12, inherit.aes = F, hjust = 1, vjust = 0) +
  #annotate(geom = "text", x = -0.005, y = 7, label = paste0("R\u00b2 = ", r2)) +
  scale_color_manual(values = c("black", "red3")) +
  scale_y_discrete("") +
  coord_cartesian(clip = "off", ylim = c(0,13)) +
  facet_wrap(facets = ~class) +
  ggtitle("Mean Verticality") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(color = "black", fill = "white"),
        panel.grid = element_line(color = "gray90"),
        strip.background = element_rect(color = "black"))



# SES VERTMEAN ----------------------------------------------------------------

load("results/sar_mods_forestOnly_forestSES/amphibians/amphibians_sar_sesvert.RData")

#amph.mods = best_mods_pred
amph.avgmod = sesvert.avg
amph.r2 = r2
amph.r2.trend = r2.trend

load("results/sar_mods_forestOnly_forestSES/mammals/mammals_sar_sesvert.RData")
#mammals.mods = best_mods_pred
mammals.avgmod = sesvert.avg
mammals.r2 = r2
mammals.r2.trend = r2.trend


load("results/sar_mods_forestOnly_forestSES/reptiles/reptiles_sar_sesvert.RData")
#rept.mods = best_mods_pred
rept.avgmod = sesvert.avg
rept.r2 = r2
rept.r2.trend = r2.trend

load("results/sar_mods_forestOnly_forestSES/birds/birds_sar_sesvert.RData")
#birds.mods = best_mods_pred
birds.avgmod = sesvert.avg
birds.r2 = r2
birds.r2.trend = r2.trend


# * - get coefficients for each class -----------------------------------------
amph.sum = summary(amph.avgmod)
amph.coefs = amph.sum$coefmat.full %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "var") %>% 
  rename(SE = "Std. Error", z_value = "z value", p_value = "Pr(>|z|)") %>% 
  mutate(class = "Amphibians")

mammals.sum = summary(mammals.avgmod)
mammals.coefs = mammals.sum$coefmat.full %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "var") %>% 
  rename(SE = "Std. Error", z_value = "z value", p_value = "Pr(>|z|)") %>% 
  mutate(class = "Mammals")

reptiles.sum = summary(rept.avgmod)
repts.coefs = reptiles.sum$coefmat.full %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "var") %>% 
  rename(SE = "Std. Error", z_value = "z value", p_value = "Pr(>|z|)") %>% 
  mutate(class = "Reptiles")

birds.sum = summary(birds.avgmod)
birds.coefs = birds.sum$coefmat.full %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "var") %>% 
  rename(SE = "Std. Error", z_value = "z value", p_value = "Pr(>|z|)") %>% 
  mutate(class = "Birds")
 

sesvert.coefs = bind_rows(amph.coefs, mammals.coefs, repts.coefs, birds.coefs) %>% 
  mutate(class = factor(class, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")),
         type = "SES Mean Verticality")

varimp.amph = data.frame(varimp = amph.avgmod$sw, class = "Amphibians") %>% 
  rownames_to_column("var")
varimp.mammals = data.frame(varimp = mammals.avgmod$sw, class = "Mammals") %>% 
  rownames_to_column("var")
varimp.rept = data.frame(varimp = rept.avgmod$sw, class = "Reptiles") %>% 
  rownames_to_column("var")
varimp.birds = data.frame(varimp = birds.avgmod$sw, class = "Birds") %>% 
  rownames_to_column("var")
varimp = rbind(varimp.amph, varimp.mammals, varimp.rept, varimp.birds)
varimp$varimp = as.numeric(varimp$varimp)

sesvert.coefs = left_join(sesvert.coefs, varimp, by = c("class", "var"))
# * - make r2 table -----------------------------------------------------------

# amph.r2 = mean(sapply(amph.mods, "[[", 4))
# mammals.r2 = mean(sapply(mammals.mods, "[[", 4))
# repts.r2 = mean(sapply(rept.mods, "[[", 4))
# birds.r2 = mean(sapply(birds.mods, "[[", 4))

sesvert.r2 = data.frame(class = c("Amphibians", "Mammals", "Reptiles", "Birds"),
                r2 = c(amph.r2, mammals.r2, rept.r2, birds.r2),
                r2.trend = c(amph.r2.trend, mammals.r2.trend, rept.r2.trend, birds.r2.trend)) %>% 
  mutate(class = factor(class, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")),
         type = "SES Mean Verticality",
         r2 = round(r2, 2),
         r2.trend = round(r2.trend, 2))




# * - plot --------------------------------------------------------------------

sesvert.coefs %>% 
  filter(var != "(Intercept)") %>%
  ggplot(aes(x = Estimate, y = reorder(var, Estimate), xmin = Estimate-SE*1.96, xmax = Estimate+SE*1.96, color = p_value<0.05)) +
  geom_pointrange() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text(data = sesvert.r2, aes(label = paste0("R\u00b2 = ", r2)), x = 0.6, y = 1.2, inherit.aes = F, hjust = 1, vjust = 0) +
  scale_color_manual(values = c("black", "red3")) +
  scale_y_discrete("") +
  facet_wrap(facets = ~class) +
  ggtitle("SES Mean Verticality") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(color = "black", fill = "white"),
        panel.grid = element_line(color = "gray90"),
        strip.background = element_rect(color = "black"))


# plot meanvert and sesmeanvert -------------------------------------------

coefs = bind_rows(sesvert.coefs, meanvert.coefs) %>% 
  filter(var != "(Intercept)") %>%
  mutate(var = factor(var, levels = c(
    "log_elev", "log_clim_velocity", "veg_complexity",
    "precip_wet", "log_precip_wet", "precip_dry", "log_precip_dry",
    "precip_sea", "log_precip_sea",
    "I(tmax_warm^2)", "tmax_warm", "tmin_cold")),
    class = factor(class, levels = c("Birds", "Mammals", "Reptiles", "Amphibians"))) %>% 
  rename("Variable Importance" = varimp)
r2 = bind_rows(sesvert.r2, meanvert.r2) %>% 
  mutate(Estimate = ifelse(type == "SES Mean Verticality", 1.0, 0.08),
         class = factor(class, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")))

library(ggh4x)
coef.plt = coefs %>%
  ggplot(aes(x = Estimate, y = var, xmin = Estimate-SE*1.96, xmax = Estimate+SE*1.96, color = p_value<0.05)) +
  geom_pointrange(aes(size = `Variable Importance`)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text(data = r2, aes(label = paste0("R\u00b2 = ", r2), x = Estimate), y = 2, inherit.aes = F, hjust = 1, vjust = 0) +
  geom_text(data = r2, aes(label = paste0("R\u00b2 trend = ", r2.trend), x = Estimate), y = 0.9, inherit.aes = F, hjust = 1, vjust = 0) +
  scale_color_manual(values = c("black", "red3")) +
  scale_size_continuous(range = c(0.02,0.5)) +
  scale_y_discrete("") +
  facet_grid2(rows = vars(type), cols = vars(class), scales = "free_x", independent = "x") +
  scale_x_facet(ROW == 1, limits = c(-0.03, 0.08)) +
  scale_x_facet(ROW == 2, limits = c(-0.4, 1.0)) +
  #facet_grid(class ~ type, scales = "free") +
  #ggtitle("SES Mean Verticality") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(color = "black", fill = "white"),
        panel.grid = element_line(color = "gray90"),
        strip.background = element_rect(color = "black"),
        legend.position = "bottom",
        legend.box.margin = margin(0,0,0,0),
        legend.box.spacing = unit(0.1, "mm"),
        legend.title.position = "top",
        legend.title = element_text(hjust = 0.5))

# png("figures/sar_coef_plots/sar_mods.png", width = 220, height = 120, res = 300, units = "mm")
# coef.plt
# dev.off()

png("figures/sar_coef_plots/sar_mods_forestOnly_forestSES.png", width = 250, height = 120, res = 300, units = "mm")
coef.plt
dev.off()

png("figures/main_figs/sar_mods_forestOnly_forestSES.png", width = 250, height = 120, res = 300, units = "mm")
coef.plt
dev.off()

# FOREST ONLY COEF PLOTS (first iteration without forest based SES and missing some wooded habitats) --------------------------------------------------

# * - VERTMEAN ----------------------------------------------------------------

load("results/sar_mods_forest_only/amphibians/amphibians_sar_vertmean2.RData")

amph.mods = best_mods_pred
amph.avgmod = vert.avg

load("results/sar_mods_forest_only/mammals/mammals_sar_vertmean2.RData")
mammals.mods = best_mods_pred
mammals.avgmod = vert.avg

load("results/sar_mods_forest_only/reptiles/reptiles_sar_vertmean2.RData")
rept.mods = best_mods_pred
rept.avgmod = vert.avg

load("results/sar_mods_forest_only/birds/birdselton_sar_vertmean.RData")
birds.mods = best_mods_pred
birds.avgmod = vert.full


# * - get coefficients for each class -----------------------------------------
amph.sum = summary(amph.avgmod)
amph.coefs = amph.sum$coefmat.full %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "var") %>% 
  rename(SE = "Std. Error", z_value = "z value", p_value = "Pr(>|z|)") %>% 
  mutate(class = "Amphibians")

mammals.sum = summary(mammals.avgmod)
mammals.coefs = mammals.sum$coefmat.full %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "var") %>% 
  rename(SE = "Std. Error", z_value = "z value", p_value = "Pr(>|z|)") %>% 
  mutate(class = "Mammals")

reptiles.sum = summary(rept.avgmod)
repts.coefs = reptiles.sum$coefmat.full %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "var") %>% 
  rename(SE = "Std. Error", z_value = "z value", p_value = "Pr(>|z|)") %>% 
  mutate(class = "Reptiles")

birds.sum = summary(birds.avgmod)
birds.coefs = birds.sum$Coef %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "var") %>% 
  rename(SE = "Std. Error", z_value = "z value", p_value = "Pr(>|z|)") %>% 
  mutate(class = "Birds")

meanvert.coefs = bind_rows(amph.coefs, mammals.coefs, repts.coefs, birds.coefs)%>% 
  mutate(class = factor(class, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")),
         type = "Mean Verticality")


# * - make r2 table -----------------------------------------------------------

amph.r2 = mean(sapply(amph.mods, "[[", 4))
mammals.r2 = mean(sapply(mammals.mods, "[[", 5))
repts.r2 = mean(sapply(rept.mods, "[[", 4))
birds.r2 = mean(sapply(birds.mods, "[[", 4))

meanvert.r2 = data.frame(class = c("Amphibians", "Mammals", "Reptiles", "Birds"),
                         r2 = c(amph.r2, mammals.r2, repts.r2, round(birds.r2,2)))%>% 
  mutate(class = factor(class, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")),
         type = "Mean Verticality")



# * - plot --------------------------------------------------------------------

meanvert.coefs %>% 
  filter(var != "(Intercept)") %>% 
  ggplot(aes(x = Estimate, y = reorder(var, Estimate), xmin = Estimate-SE*1.96, xmax = Estimate+SE*1.96, color = p_value<0.05)) +
  geom_pointrange() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text(data = meanvert.r2, aes(label = paste0("R\u00b2 = ", r2)), x = 0.042, y = 1.2, inherit.aes = F, hjust = 1, vjust = 0) +
  #annotate(geom = "text", x = -0.005, y = 7, label = paste0("R\u00b2 = ", r2)) +
  scale_color_manual(values = c("black", "red3")) +
  scale_y_discrete("") +
  facet_wrap(facets = ~class) +
  ggtitle("Mean Verticality") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(color = "black", fill = "white"),
        panel.grid = element_line(color = "gray90"),
        strip.background = element_rect(color = "black"))


# SES VERTMEAN ----------------------------------------------------------------

load("results/sar_mods_forest_only/amphibians/amphibians_sar_sesvert2.RData")

amph.mods = best_mods_pred
amph.avgmod = sesvert.avg

load("results/sar_mods_forest_only/mammals/mammals_sar_sesvert2.RData")
mammals.mods = best_mods_pred
mammals.avgmod = sesvert.avg

load("results/sar_mods_forest_only/reptiles/reptiles_sar_sesvert2.RData")
rept.mods = best_mods_pred
rept.avgmod = sesvert.avg

load("results/sar_mods_forest_only/birds/birdselton_sar_sesvert2.RData")
birds.mods = best_mods_pred
birds.avgmod = sesvert.full


# * - get coefficients for each class -----------------------------------------
amph.sum = summary(amph.avgmod)
amph.coefs = amph.sum$coefmat.full %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "var") %>% 
  rename(SE = "Std. Error", z_value = "z value", p_value = "Pr(>|z|)") %>% 
  mutate(class = "Amphibians")

mammals.sum = summary(mammals.avgmod)
mammals.coefs = mammals.sum$coefmat.full %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "var") %>% 
  rename(SE = "Std. Error", z_value = "z value", p_value = "Pr(>|z|)") %>% 
  mutate(class = "Mammals")

reptiles.sum = summary(rept.avgmod)
repts.coefs = reptiles.sum$coefmat.full %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "var") %>% 
  rename(SE = "Std. Error", z_value = "z value", p_value = "Pr(>|z|)") %>% 
  mutate(class = "Reptiles")

birds.sum = summary(birds.avgmod)
birds.coefs = birds.sum$Coef %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "var") %>% 
  rename(SE = "Std. Error", z_value = "z value", p_value = "Pr(>|z|)") %>% 
  mutate(class = "Birds")

sesvert.coefs = bind_rows(amph.coefs, mammals.coefs, repts.coefs, birds.coefs) %>% 
  mutate(class = factor(class, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")),
         type = "SES Mean Verticality")


# * - make r2 table -----------------------------------------------------------

amph.r2 = mean(sapply(amph.mods, "[[", 4))
mammals.r2 = mean(sapply(mammals.mods, "[[", 4))
repts.r2 = mean(sapply(rept.mods, "[[", 4))
birds.r2 = mean(sapply(birds.mods, "[[", 4))

sesvert.r2 = data.frame(class = c("Amphibians", "Mammals", "Reptiles", "Birds"),
                        r2 = c(amph.r2, mammals.r2, repts.r2, round(birds.r2,2))) %>% 
  mutate(class = factor(class, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")),
         type = "SES Mean Verticality")




# * - plot --------------------------------------------------------------------

sesvert.coefs %>% 
  filter(var != "(Intercept)") %>%
  ggplot(aes(x = Estimate, y = reorder(var, Estimate), xmin = Estimate-SE*1.96, xmax = Estimate+SE*1.96, color = p_value<0.05)) +
  geom_pointrange() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text(data = sesvert.r2, aes(label = paste0("R\u00b2 = ", r2)), x = 0.6, y = 1.2, inherit.aes = F, hjust = 1, vjust = 0) +
  scale_color_manual(values = c("black", "red3")) +
  scale_y_discrete("") +
  facet_wrap(facets = ~class) +
  ggtitle("SES Mean Verticality") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(color = "black", fill = "white"),
        panel.grid = element_line(color = "gray90"),
        strip.background = element_rect(color = "black"))


# plot meanvert and sesmeanvert -------------------------------------------

coefs = bind_rows(sesvert.coefs, meanvert.coefs) %>% 
  filter(var != "(Intercept)") %>%
  mutate(var = factor(var, levels = c(
    "log_elev", "log_clim_velocity", "veg_complexity",
    "precip_wet", "log_precip_wet", "precip_dry", "log_precip_dry",
    "precip_sea", "log_precip_sea",
    "tmax_warm", "tmin_cold")))
r2 = bind_rows(sesvert.r2, meanvert.r2) %>% 
  mutate(Estimate = ifelse(type == "SES Mean Verticality", 0.75, 0.06))

coef.plt = coefs  %>% 
  ggplot(aes(x = Estimate, y = var, xmin = Estimate-SE*1.96, xmax = Estimate+SE*1.96, color = p_value<0.05)) +
  geom_pointrange() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text(data = r2, aes(label = paste0("R\u00b2 = ", r2), x = Estimate), y = 0.9, inherit.aes = F, hjust = 1, vjust = 0) +
  scale_color_manual(values = c("black", "red3")) +
  scale_y_discrete("") +
  facet_grid2(rows = vars(type), cols = vars(class), scales = "free_x", independent = "x") +
  scale_x_facet(ROW == 1, limits = c(-0.06, 0.06)) +
  scale_x_facet(ROW == 2, limits = c(-0.75, 0.75)) +
  #facet_grid(class ~ type, scales = "free") +
  #ggtitle("SES Mean Verticality") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(color = "black", fill = "white"),
        panel.grid = element_line(color = "gray90"),
        strip.background = element_rect(color = "black"),
        legend.position = "bottom",
        legend.box.margin = margin(0,0,0,0),
        legend.box.spacing = unit(0.1, "mm"))

png("figures/sar_coef_plots/sar_mods_forest_only.png", width = 220, height = 120, res = 300, units = "mm")
coef.plt
dev.off()

