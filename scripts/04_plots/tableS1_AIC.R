amph = readRDS("results/sdmTMB_models2/amphibians_sesvert.rds")
birds = readRDS("results/sdmTMB_models2/birds_sesvert.rds")
reptiles = readRDS("results/sdmTMB_models2/reptiles_sesvert.rds")
mammals = readRDS("results/sdmTMB_models2/mammals_sesvert.rds")

amph = amph$compMods_aic$modsel
birds = birds$compMods_aic$modsel
reptiles = reptiles$compMods_aic$modsel
mammals = mammals$compMods_aic$modsel

df = rbind(birds, mammals, reptiles, amph) %>% 
  mutate(taxon = factor(taxon, levels = c("Birds", "Mammals", "Reptiles", "Amphibians")))


mods = data.frame(mods = c(
  "f1",
       "f1, . ~ .-precip_warm:canopy_height",
       "f1, . ~ .-log_precip_dry:canopy_height",
       "f1, . ~ .-tmin_cold:canopy_height",
       "f1, . ~ .-tmax_warm:canopy_height - I(tmax_warm^2):canopy_height",
       "f1, . ~ .-precip_warm:canopy_height - log_precip_dry:canopy_height",
       "f1, . ~ .-precip_warm:canopy_height - tmin_cold:canopy_height",
       "f1, . ~ .-precip_warm:canopy_height - tmax_warm:canopy_height - I(tmax_warm^2):canopy_height",
       "f1, . ~ .-log_precip_dry:canopy_height - tmin_cold:canopy_height",
       "f1, . ~ .-log_precip_dry:canopy_height - tmax_warm:canopy_height - I(tmax_warm^2):canopy_height",
       "f1, . ~ .-tmin_cold:canopy_height - tmax_warm:canopy_height - I(tmax_warm^2):canopy_height",
       "f1, . ~ .-precip_warm:canopy_height - log_precip_dry:canopy_height - tmin_cold:canopy_height",
       "f1, . ~ .-precip_warm:canopy_height - log_precip_dry:canopy_height - tmax_warm:canopy_height - I(tmax_warm^2):canopy_height",
       "f1, . ~ .-log_precip_dry:canopy_height - tmin_cold:canopy_height - tmax_warm:canopy_height - I(tmax_warm^2):canopy_height",
       "f1, . ~ .-precip_warm:canopy_height - log_precip_dry:canopy_height - tmin_cold:canopy_height - tmax_warm:canopy_height - I(tmax_warm^2):canopy_height"
  ),
  modnum = as.character(1:15))

df = left_join(df, mods, by = c("mod" = "modnum"))

write.csv(df, "results/sdmTMB_models2/tableS1_AIC.csv", row.names = F)
