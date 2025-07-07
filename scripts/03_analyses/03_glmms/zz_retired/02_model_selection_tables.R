# Make tables of model selection objectsion
# Date: July 18, 2024
# Author: Lydia Soifer

library(sdmTMB)
library(foreach)


# sanity check models -----------------------------------------------------

load("results/sdmTMB_models/amphibians_meanvert.RData")
load("results/sdmTMB_models/model_selection/amphibians_sesvert.RData")
sanity(mod) # ok

load("results/sdmTMB_models/amphibians_sesvert.RData")
sanity(mod) # ok

load("results/sdmTMB_models/reptiles_meanvert.RData")
sanity(mod) #ok

load("results/sdmTMB_models/reptiles_sesvert.RData")
sanity(mod) #ok

load("results/sdmTMB_models/mammals_meanvert.RData")
sanity(mod) # ok (sigma G smaller than 0.01 indicates random intercept isn't adding much to the model)

load("results/sdmTMB_models/mammals_sesvert.RData")
sanity(mod) # ok

load("results/sdmTMB_models/birds_meanvert.RData")
sanity(mod) # ok

load("results/sdmTMB_models/birds_sesvert.RData")
load("results/sdmTMB_models/model_selection/birds_sesvert.RData")
sanity(mod) # ok

rm(mod, pred, pred.f)

# AIC based model selection -----------------------------------------------
# files should be loaded and removed one at a time becasue they are large
# check convergence for best model from AIC

files = list.files("results/sdmTMB_models/model_selection/", pattern = ".RData", full.names = T)

comps = foreach(i = 1:length(files), .combine = "rbind") %do% {
  load(files[i])
  comp.aic = compMods_aic[[2]]
  comp.cv = compMods_cv[[3]]
  comp = inner_join(comp.aic, comp.cv, by = c("model", "taxon", "response_var"))
  gc()
}

write.csv(comps, "tables/supplement/model_selection_compare.csv", row.names = F)













