library(tidyverse)
library(lemon)
library(ggh4x)

# VERTMEAN ----------------------------------------------------------------

load("results/amphibians_sar_vertmean2.RData")

amph.mods = best_mods_pred
amph.avgmod = vert.avg

load("results/mammals_sar_vertmean2.RData")
mammals.mods = best_mods_pred
mammals.avgmod = vert.avg

load("results/reptiles_sar_vertmean2.RData")
rept.mods = best_mods_pred
rept.avgmod = vert.avg

load("results/birdselton_sar_vertmean.RData")
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

load("results/amphibians_sar_sesvert2.RData")

amph.mods = best_mods_pred
amph.avgmod = sesvert.avg

load("results/mammals_sar_sesvert2.RData")
mammals.mods = best_mods_pred
mammals.avgmod = sesvert.avg

load("results/reptiles_sar_sesvert2.RData")
rept.mods = best_mods_pred
rept.avgmod = sesvert.avg

load("results/birdselton_sar_sesvert2.RData")
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
    "I(tmax_warm^2)", "tmax_warm", "tmin_cold")))
r2 = bind_rows(sesvert.r2, meanvert.r2) %>% 
  mutate(Estimate = ifelse(type == "SES Mean Verticality", 0.65, 0.022))

coef.plt = coefs %>%
  ggplot(aes(x = Estimate, y = var, xmin = Estimate-SE*1.96, xmax = Estimate+SE*1.96, color = p_value<0.05)) +
  geom_pointrange() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text(data = r2, aes(label = paste0("R\u00b2 = ", r2), x = Estimate), y = 0.9, inherit.aes = F, hjust = 1, vjust = 0) +
  scale_color_manual(values = c("black", "red3")) +
  scale_y_discrete("") +
  facet_grid2(rows = vars(type), cols = vars(class), scales = "free_x", independent = "x") +
  scale_x_facet(ROW == 1, limits = c(-0.022, 0.022)) +
  scale_x_facet(ROW == 2, limits = c(-0.62, 0.62)) +
  #facet_grid(class ~ type, scales = "free") +
  #ggtitle("SES Mean Verticality") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(color = "black", fill = "white"),
        panel.grid = element_line(color = "gray90"),
        strip.background = element_rect(color = "black"),
        legend.position = "bottom",
        legend.box.margin = margin(0,0,0,0),
        legend.box.spacing = unit(0.1, "mm"))

png("figures/sar_coef_plots/sar_mods.png", width = 220, height = 120, res = 300, units = "mm")
coef.plt
dev.off()



# FOREST ONLY COEF PLOTS --------------------------------------------------

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

