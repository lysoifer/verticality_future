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
         term = case_when(term == "tmin_cold:canopy_height" ~ "Tmin[cold]:canopy height",
                          term == "tmax_warm:canopy_height" ~ "Tmax[warm]:canopy height",
                          term == "I(tmax_warm^2):canopy_height" ~ "Tmax[warm]\u00b2:canopy height",
                          term == "log_precip_dry:canopy_height" ~ "log(Precip[dry]):canopy height",
                          term == "precip_warm:canopy_height" ~ "Precip[warm]:canopy height",
                          term == "log_clim_velocity" ~ "log(climate velocity)", 
                          term == "veg_den" ~ "Vegetation density",
                          term == "canopy_height" ~ "Canopy height",
                          term == "precip_warm" ~ "Precip[warm]",
                          term == "precip_wet" ~ "Precip[wet]",
                          term == "log_precip_dry" ~ "log(Precip[dry])",
                          term == "I(tmax_warm^2)" ~ "(Tmax[warm])\u00b2",
                          term == "tmax_warm" ~ "Tmax[warm]",
                          term == "tmin_cold" ~ "Tmin[cold]"),
        term = factor(term, levels = c(
          "Tmin[cold]:canopy height",
          "Tmax[warm]:canopy height",
          "Tmax[warm]\u00b2:canopy height",
          "log(Precip[dry]):canopy height",
          "Precip[warm]:canopy height",
          "log(climate velocity)",
          "Vegetation density",
          "Canopy height",
          "Precip[warm]",
          "Precip[wet]",
          "log(Precip[dry])",
          "(Tmax[warm])\u00b2",
          "Tmax[warm]",
          "Tmin[cold]")),
        sig = ifelse(conf.low*conf.high > 0, "sig", "notsig"))




# plot model coefficients
mods = mods %>%
  mutate(lty = ifelse(sig == "sig", "solid", "dashed"),
         taxacol = case_when(taxa == "Birds" ~ "royalblue3",
                            taxa == "Mammals" ~ "plum3",
                            taxa == "Reptiles" ~ "palegreen3",
                            taxa == "Amphibians" ~ "gold2"),
         fillcol = ifelse(sig == "sig", taxacol, "white"),
         taxa = factor(taxa, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")),
         taxacol = factor(taxacol, levels = c("royalblue3", "plum3", "palegreen3", "gold2")))

mods %>% 
  ggplot(aes(x = estimate, y = term, xmin = conf.low, xmax = conf.high, fill = fillcol, color = taxacol, linetype = lty)) +
  geom_pointrange(size = 0.5, pch = 21, position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_discrete("",
                   labels = c("Tmin[cold]:canopy height" = expression("Tmin"[cold]*":Canopy height"),
                              "Tmax[warm]:canopy height" = expression("Tmax"[warm]*":Canopy height"),
                              "Tmax[warm]\u00b2:canopy height" = expression("(Tmax"["warm"]*")\u00b2:Canopy height"),
                              "log(Precip[dry]):canopy height" = expression("log(Precip"["dry"]*"):Canopy height"),
                              "Precip[warm]:canopy height" = expression("Precip"[warm]*":Canopy height"),
                              "Precip[warm]" = expression("Precip"[warm]),
                              "Precip[wet]" = expression("Precip"[wet]),
                              "log(Precip[dry])" = expression("log(Precip"["dry"]*")"),
                              "(Tmax[warm])\u00b2" = expression("(Tmax"["warm"]*")\u00b2"),
                              "Tmax[warm]" = expression("Tmax"[warm]),
                              "Tmin[cold]" = expression("Tmin"[cold]))) +
  scale_x_continuous(limits = c(-0.2,0.4)) +
  scale_fill_identity(guide = "legend",
                      labels = levels(mods$taxa),
                      breaks = levels(mods$taxacol)) +
  scale_color_identity() +
  scale_linetype_identity() +
  # facet_grid2(rows = vars(type), cols = vars(class), scales = "free_x", independent = "x") +
  # scale_x_facet(ROW == 1, limits = c(-0.08, 0.18)) +
  # scale_x_facet(ROW == 2, limits = c(-0.2, 0.7)) +
  labs(x = "Effect size of standardized predictors") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(color = "black", fill = "white"),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(color = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.842,0.135),
        legend.box.margin = margin(0,0,0,0),
        legend.box.spacing = unit(0.1, "mm"),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.background = element_rect(color = "black", fill = "white"),
        axis.text.x = element_text(size = 6),
        panel.grid.major.x = element_blank())

ggsave("figures/main_figs/fig3_sdmTMB_coefs_only.png", width = 140, height = 120, units = "mm", dpi = 300)

# Add interaction plots to coef plot --------------------------------------

# amph.df = readRDS("results/sdmTMB_models2/marginal_effects/canopy_interactions/amph_interactions.rds")
# amph.df = rbindlist(amph.df) %>% mutate(taxa = "Amphibians")
# 
# rept.df = readRDS("results/sdmTMB_models2/marginal_effects/canopy_interactions/repts_interactions.rds")
# rept.df = rbindlist(rept.df) %>% mutate(taxa = "Reptiles")
# 
# mammals.df = readRDS("results/sdmTMB_models2/marginal_effects/canopy_interactions/mammals_interactions.rds")
# mammals.df = rbindlist(mammals.df) %>% mutate(taxa = "Mammals")
# 
# birds.df = readRDS("results/sdmTMB_models2/marginal_effects/canopy_interactions/birds_interactions.rds")
# birds.df = rbindlist(birds.df) %>% mutate(taxa = "Birds")
# 
# df = bind_rows(amph.df, rept.df, mammals.df, birds.df) %>% 
#   filter(x >= -1.3 & x <= 3.5) %>% 
#   filter(group == -2 | group == 0 | group == 2) %>% 
#   mutate(term = case_when(term == "canopy-log_precip_dry" ~ "log(Precip dry)",
#                           term == "canopy-tmax_warm" ~ "Tmax warm",
#                           term == "canopy-tmin_cold" ~ "Tmin cold",
#                           term ==  "canopy-precip_warm" ~ "Precip warm"),
#          taxa = factor(taxa, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")))
# 
# # select only significant interactions
# df = df %>% 
#   filter((taxa == "Birds" & (term == "Precip warm" | term == "Tmax warm" | term == "Tmin cold")) |
#            (taxa == "Mammals" & term != "Tmin cold") |
#            (taxa == "Reptiles" & term == "Tmin cold") |
#            (taxa == "Amphibians" & (term == "Tmax warm" | term == "Tmin cold")))
# 
# 
# groups = as.data.frame(df %>% distinct(taxa, term))
# taxa = c()
# term = c()
# plots = list()
# 
# for(i in 1:nrow(groups)) {
#   t = as.character(groups[i,1])
#   tr = as.character(groups[i,2])
#   d = df %>% filter(taxa == t & term == tr)
#   p = ggplot(d, aes(x, predicted, color = group, fill = group)) +
#     geom_line() +
#     geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
#     scale_y_continuous("Predicted SES Verticality") +
#     scale_x_continuous("Scaled canopy height") +
#     scale_color_discrete(guide = guide_legend("Scaled interaction term")) +
#     scale_fill_discrete(guide = guide_legend("Scaled interaction term")) +
#     annotate(geom = "text", label = tr, x = -Inf, y = Inf, hjust = -0.1, vjust = 2) +
#     theme_bw() +
#     theme(legend.position = "none")
#   taxa = c(taxa, t)
#   term = c(term, tr)
#   plots[[i]] = p
# }
# 
# plot_data = tibble(taxa = taxa, term = term, plots = plots)
# 
# coef.taxa = c()
# coef.plt = list()
# for(i in c("Amphibians", "Birds", "Mammals", "Reptiles")) {
#   d = mods %>% filter(taxa == i)
#   p = d %>% 
#     ggplot(aes(x = estimate, y = term, xmin = conf.low, xmax = conf.high, color = sig)) +
#     geom_pointrange(size = 0.15) +
#     geom_vline(xintercept = 0, linetype = "dashed") +
#     scale_x_continuous(limits = c(-0.2,0.3), breaks = seq(-0.2,0.3, 0.1)) +
#     scale_y_discrete("",
#                      labels = c("Tmin[cold]:canopy height" = expression("Tmin"[cold]*":Canopy height"),
#                                 "Tmax[warm]:canopy height" = expression("Tmax"[warm]*":Canopy height"),
#                                 "Tmax[warm]\u00b2:canopy height" = expression("(Tmax"["warm"]*")\u00b2:Canopy height"),
#                                 "log(Precip[dry]):canopy height" = expression("log(Precip"["dry"]*"):Canopy height"),
#                                 "Precip[warm]:canopy height" = expression("Precip"[warm]*":Canopy height"),
#                                 "Precip[warm]" = expression("Precip"[warm]),
#                                 "Precip[wet]" = expression("Precip"[wet]),
#                                 "log(Precip[dry])" = expression("log(Precip"["dry"]*")"),
#                                 "(Tmax[warm])\u00b2" = expression("(Tmax"["warm"]*")\u00b2"),
#                                 "Tmax[warm]" = expression("Tmax"[warm]),
#                                 "Tmin[cold]" = expression("Tmin"[cold]))) +
#     scale_color_manual(values = c("black", "red3")) +
#     labs(x = "Effect size of\nstandardized predictors") +
#     theme(plot.title = element_text(hjust = 0.5),
#           panel.background = element_rect(color = "black", fill = NA),
#           panel.grid.major = element_line(color = "gray90"),
#           panel.grid.minor = element_blank(),
#           strip.background = element_rect(color = "black"),
#           legend.position = "none",
#           legend.box.margin = margin(0,0,0,0),
#           legend.box.spacing = unit(0.1, "mm"),
#           legend.title.position = "top",
#           legend.title = element_text(hjust = 0.5),
#           axis.text = element_text(size = 10))
#   
#   coef.taxa = c(coef.taxa, i)
#   coef.plt[[i]] = p
# }
# 
# coef.plotdata = tibble(taxa = coef.taxa, plots = coef.plt)
# 
# pcoef = coef.plotdata$plots
# presp = plot_data$plots
# 
# 
# row1 = wrap_elements(panel = textGrob("Birds", rot = 90))
# row2 = wrap_elements(panel = textGrob("Mammals", rot = 90))
# row3 = wrap_elements(panel = textGrob("Reptiles", rot = 90))
# row4 = wrap_elements(panel = textGrob("Amphibians", rot = 90))
# 
# 
# coefplots = pcoef[[2]] + pcoef[[3]] + pcoef[[4]] + pcoef[[1]] + plot_layout(ncol = 1, axes = "collect_x")
# rowlabs = row1 + row2 + row3 + row4 + plot_layout(ncol = 1)
# 
# design = "
# abc
# def
# ghi
# jkl"
# 
# respplots = presp[[7]] + presp[[8]] + presp[[9]] +
#  presp[[5]] + presp[[6]] + presp[[4]] +
#   presp[[3]] + plot_spacer() + plot_spacer() +
#   presp[[1]] + presp[[2]] + plot_spacer() +
#   plot_layout(design = design, byrow = T, axes = "collect")
# 
# coefplots.grob = patchworkGrob(coefplots)
# rowlabs.grob = patchworkGrob(rowlabs)
# respplots.grob = patchworkGrob(respplots)
# 
# d = df %>% filter(taxa == t & term == tr)
# legend = ggplot(d, aes(x, predicted, color = group, fill = group)) +
#   geom_line() +
#   geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
#   scale_y_continuous("Predicted SES Verticality") +
#   scale_x_continuous("Scaled canopy height") +
#   scale_color_discrete(guide = guide_legend("Scaled interaction term")) +
#   scale_fill_discrete(guide = guide_legend("Scaled interaction term")) +
#   annotate(geom = "text", label = tr, x = -Inf, y = Inf, hjust = -0.1, vjust = 2) +
#   theme_bw() +
#   theme(legend.position = "bottom")
# legend = ggpubr::get_legend(legend)
# 
# 
# layout = rbind(c(1,2,3),
#                c(NA,NA, 4))
# gridExtra::grid.arrange(rowlabs.grob, coefplots.grob, respplots.grob, legend, ncol = 3, widths = c(0.1, 1.2, 2.3), heights = c(4, 0.1), layout_matrix = layout)
# 
# g = gridExtra::arrangeGrob(rowlabs.grob, coefplots.grob, respplots.grob, legend, ncol = 3, widths = c(0.2, 1.2, 2.3), heights = c(4, 0.1), layout_matrix = layout)
# ggsave("figures/main_figs/fig3/fig3_sdmTMB_coefs-interactions.png", width = 300, height = 250, dpi = 300, units = "mm", plot = g)
# 



