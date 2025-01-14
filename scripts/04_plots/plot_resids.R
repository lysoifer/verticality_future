library(DHARMa)
library(grid)
library(patchwork)
source("scripts/00_functions/dharma_plot_resids.R")

plot_resids = function(mod, pred, response_var, taxa, xlab) {
  # mod: sdmTMB model to check residuals on
  # response variable used in the model
  # file path without file extension
  
  #pred = predict(mod, type = "response")
  s = simulate(mod, nsim = 500, seed = 12345, type = "mle-mvn")
  
  # simulate quantile residuals
  r <- DHARMa::createDHARMa(
    simulatedResponse = s,
    observedResponse = pred[,response_var],
    fittedPredictedResponse = pred[,"est"]
  )
  
  par(bg = "#00000000")
  plotResiduals2(r, rank = F, form = pred[,response_var], xlab = xlab, main = taxa)
  p = recordPlot()
  
  return(p)
  
}

load("results/sdmTMB_models/amphibians_sesvert.RData")
amph.sesvert = plot_resids(mod, pred, response_var = "vert.mean.ses", taxa = "Amphibians", xlab = "SES verticality")

load("results/sdmTMB_models/amphibians_meanvert.RData")
amph.meanvert = plot_resids(mod, pred, response_var = "vert.mean", taxa = "Amphibians", xlab = "Mean verticality")

load("results/sdmTMB_models/reptiles_sesvert.RData")
rept.sesvert =  plot_resids(mod, pred, response_var = "vert.mean.ses", taxa = "Reptiles", xlab = "SES verticality")
                        
load("results/sdmTMB_models/reptiles_meanvert.RData")
rept.meanvert = plot_resids(mod, pred, response_var = "vert.mean", taxa = "Reptiles", xlab = "Mean verticality")

load("results/sdmTMB_models/mammals_sesvert.RData")
mammals.sesvert = plot_resids(mod, pred, response_var = "vert.mean.ses", taxa = "Mammals", xlab = "SES verticality")

load("results/sdmTMB_models/mammals_meanvert.RData")
mammals.meanvert = plot_resids(mod, pred, response_var = "vert.mean", taxa = "Mammals", xlab = "Mean verticality")

load("results/sdmTMB_models/birds_sesvert.RData")
birds.sesvert =  plot_resids(mod, pred, response_var = "vert.mean.ses", taxa = "Birds", xlab = "SES verticality")

load("results/sdmTMB_models/birds_meanvert.RData")
birds.meanvert = plot_resids(mod, pred, response_var = "vert.mean", taxa = "Birds", xlab = "Mean verticality")



library(ggpubr)
library(gridGraphics)
png("figures/supp_figs/resid_plots.png", height = 400, width = 300, res = 300, units = "mm")
ggarrange(plotlist = list(birds.sesvert, birds.meanvert,
                          mammals.sesvert, mammals.meanvert,
                          rept.sesvert, rept.meanvert,
                          amph.sesvert, amph.meanvert), nrow = 4, ncol = 2)
dev.off()


