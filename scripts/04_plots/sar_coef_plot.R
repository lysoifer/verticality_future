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

coefs = bind_rows(sesvert.coefs, meanvert.coefs)
r2 = bind_rows(sesvert.r2, meanvert.r2) %>% 
  mutate(Estimate = ifelse(type == "SES Mean Verticality", 0.65, 0.022))

coefs %>% 
  filter(var != "(Intercept)") %>%
  ggplot(aes(x = Estimate, y = reorder(var, Estimate), xmin = Estimate-SE*1.96, xmax = Estimate+SE*1.96, color = p_value<0.05)) +
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
        strip.background = element_rect(color = "black"))


