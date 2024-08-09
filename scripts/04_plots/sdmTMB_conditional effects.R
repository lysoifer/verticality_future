library(sdmTMB)
library(foreach)
library(doParallel)
library(visreg)
library(ggeffects)

mods = list.files("results/sdmTMB_models/", pattern = ".RData", full.names = T)

plot_cond2d = function(modpath) {
  load(modpath)
  p2d.ses = visreg2d(mod, xvar = "log_precip_dry", yvar = "tmax_warm", scale = "response", plot.type = "gg", data = mod$data)
  p2d.ses2 = visreg2d(mod, xvar = "log_precip_dry", yvar = "tmin_cold", scale = "response", plot.type = "gg", data = mod$data)
  return(list(p2d.ses, p2d.ses2))
}

plts = foreach(i = 1:length(mods)) %do% {
 p = plot_cond2d(mods[i])
 p
}
names(plts) = mods

row1 = wrap_elements(panel = textGrob(expression(bold("Birds")), rot = 90, gp = gpar(fonface = "bold")))
row2 = wrap_elements(panel = textGrob(expression(bold("Mammals")), rot = 90, gp = gpar(fonface = "bold")))
row3 = wrap_elements(panel = textGrob(expression(bold("Reptiles")), rot = 90, gp = gpar(fonface = "bold")))
row4 = wrap_elements(panel = textGrob(expression(bold("Amphibians")), rot = 90, gp = gpar(fonface = "bold")))

plts.meanvert = list(row1, b1 = plts[[3]][[1]], b2 = plts[[3]][[2]],
                     row2, m1 = plts[[5]][[1]], m2 = plts[[5]][[2]],
                     row3, r1 = plts[[7]][[1]], r2 = plts[[7]][[2]],
                     row4, a1 = plts[[1]][[1]], a2 = plts[[1]][[2]])

plot = wrap_plots(plts.meanvert) +
  plot_layout(ncol = 3, nrow = 4, byrow = T, 
              widths = c(0.1,1,1),
              axis_titles = "collect_x") +
  plot_annotation(title = "Mean Verticality", 
                  tag_levels = list(c("", "A", "B", 
                                 "", "C", "D",
                                 "", "E", "F",
                                 "", "G", "H"))) &
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.tag.position = c(0.02,0.98))

png("figures/supp_figs/cond_effects_meanvert.png", width = 225, height = 300, res = 300, units = "mm")
plot
dev.off()


plts.sesvert = list(row1, b1 = plts[[4]][[1]], b2 = plts[[4]][[2]],
                     row2, m1 = plts[[6]][[1]], m2 = plts[[6]][[2]],
                     row3, r1 = plts[[8]][[1]], r2 = plts[[8]][[2]],
                     row4, a1 = plts[[2]][[1]], a2 = plts[[2]][[2]])

plot = wrap_plots(plts.sesvert) +
  plot_layout(ncol = 3, nrow = 4, byrow = T, 
              widths = c(0.1,1,1),
              axis_titles = "collect_x") +
  plot_annotation(title = "SES Verticality", 
                  tag_levels = list(c("", "A", "B", 
                                      "", "C", "D",
                                      "", "E", "F",
                                      "", "G", "H"))) &
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.tag.position = c(0.02,0.98))

png("figures/supp_figs/cond_effects_sesvert.png", width = 225, height = 300, res = 300, units = "mm")
plot
dev.off()