# Models of richness by verticality
library(tidyverse)
library(sdmTMB)
library(lme4)
library(terra)
library(car)
library(data.table)
library(DHARMa)
source("scripts/00_functions/manuscript_functions.R")


amph = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/amph_comdat.csv")
birds = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/birds_comdat.csv")
rept = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/reptiles_comdat.csv")
mammals = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/mammals_comdat.csv")

amph$taxa = "Amphibians"
birds$taxa = "Birds"
rept$taxa = "Reptiles"
mammals$taxa = "Mammals"

amph = amph %>% 
  # subset to richness >= 5
  filter(rich >= 5) %>% 
  filter(!is.na(vert.mean.ses)) %>% 
  drop_na() %>% 
  dplyr::select(x,y,rich, vert.mean, vert.mean.ses, ecoregion, biome, realm) %>% 
  mutate(x = x/1e5,
         y = y/1e5,
         ecoregion = factor(ecoregion),
         biome = factor(biome),
         biorealm = factor(paste(biome, realm, sep = "_")),
         realm = factor(realm))

birds = birds %>% 
  # subset to richness >= 5
  filter(rich >= 5) %>% 
  filter(!is.na(vert.mean.ses)) %>% 
  drop_na() %>% 
  dplyr::select(x,y,rich, vert.mean, vert.mean.ses, ecoregion, biome, realm) %>% 
  mutate(x = x/1e5,
         y = y/1e5,
         ecoregion = factor(ecoregion),
         biome = factor(biome),
         biorealm = factor(paste(biome, realm, sep = "_")),
         realm = factor(realm))

rept = rept %>% 
  # subset to richness >= 5
  filter(rich >= 5) %>% 
  filter(!is.na(vert.mean.ses)) %>% 
  drop_na() %>% 
  dplyr::select(x,y,rich, vert.mean, vert.mean.ses, ecoregion, biome, realm) %>% 
  mutate(x = x/1e5,
         y = y/1e5,
         ecoregion = factor(ecoregion),
         biome = factor(biome),
         biorealm = factor(paste(biome, realm, sep = "_")),
         realm = factor(realm))

mammals = mammals %>% 
  # subset to richness >= 5
  filter(rich >= 5) %>% 
  filter(!is.na(vert.mean.ses)) %>% 
  drop_na() %>% 
  dplyr::select(x,y,rich, vert.mean, vert.mean.ses, ecoregion, biome, realm) %>% 
  mutate(x = x/1e5,
         y = y/1e5,
         ecoregion = factor(ecoregion),
         biome = factor(biome),
         biorealm = factor(paste(biome, realm, sep = "_")),
         realm = factor(realm))


f1 = formula(rich ~ vert.mean.ses)


# SES VERT ----------------------------------------------------------------

vert_richness_mod = function(f, dat, family, wts = NULL, edge.ratio) {
  
  # fit mesh
  if(nrow(dat) > 1000) {samp = dat %>% sample_n(1000)} else {samp = dat}
  samp.cor = ncf::correlog(x = samp$x, y = samp$y, z = samp$rich, increment = 50000/1e5, resamp = 99)
  
  fitmesh = fit_mesh_richness(f1, dat, range = samp.cor$x.intercept, v = v, family = family, edge.ratio = edge.ratio)
  
  mesh = fitmesh$meshes[[length(fitmesh$meshes)]]
  
  all_ok = FALSE
  iter = 1
  while(!all_ok & iter <= 2) {
    print(iter)
    if(iter == 1) {
      mod = sdmTMB(f, 
                   data = dat,
                   weights = wts,
                   mesh = mesh,
                   spatial = "on",
                   reml = T,
                   family = family,
                   control = sdmTMBcontrol(eval.max = 6000, iter.max = 3000))
    } else {
      mod = sdmTMB(f, 
                   data = dat,
                   weights = wts,
                   mesh = mesh,
                   spatial = "on",
                   reml = T,
                   family = family,
                   control = sdmTMBcontrol(eval.max = 6000, iter.max = 3000,
                                           start = list(ln_phi = ln_phi,
                                                        ln_kappa = ln_kappa,
                                                        ln_tau_O = ln_tau_o)))
    }
    all_ok = sanity(mod)$all_ok
    pars = get_pars(mod)
    ln_phi = pars$ln_phi
    ln_kappa = pars$ln_kappa
    ln_tau_o = pars$ln_tau_O
    iter = iter+1
    gc()
  }
  return(mod)
}

# * richness ~ SES verticality ------------------------------------------------

# models account for spatial autocorrelation, but do not include realm as a random intercept or any spatially varying coefficient
# if models error, try starting to run a different model and then rerun the one erroring - not sure why this works, but it does

amph.mod = vert_richness_mod(f1, dat = amph, family = gaussian(), edge.ratio = 0.35) # had to increase edge ratio to get model to converge
sanity(amph.mod)
saveRDS(amph.mod, file = "results/sdmTMB_models/richness_vert_models/sesvert/amph_mod.rds")
mod = readRDS("results/sdmTMB_models/richness_vert_models/sesvert/amph_mod.rds")

birds.mod = vert_richness_mod(f1, dat = birds, family = gaussian(), edge.ratio = 0.4) # had to increase edge ratio to get model to converge
sanity(birds.mod)
saveRDS(birds.mod, file = "results/sdmTMB_models/richness_vert_models/sesvert/birds_mod.rds")
mod = readRDS("results/sdmTMB_models/richness_vert_models/sesvert/birds_mod.rds")

rept.mod = vert_richness_mod(f1, dat = rept, family = gaussian(), edge.ratio = 0.2)
sanity(rept.mod)
saveRDS(rept.mod, file = "results/sdmTMB_models/richness_vert_models/sesvert/rept_mod.rds")

mammals.mod = vert_richness_mod(f1, dat = mammals, family = gaussian(), edge.ratio = 0.3) # had to increase edge ratio to fix gradient issues
sanity(mammals.mod)
saveRDS(mammals.mod, file = "results/sdmTMB_models/richness_vert_models/sesvert/mammals_mod.rds")

summary(amph.mod)
summary(birds.mod)
summary(rept.mod)
summary(mammals.mod)

# now run the models for individual biomes
for(b in unique(amph$biome)) {
  d = amph %>% filter(biome == b)
  san = F
  edge.ratio = 0.2
  while(!san & edge.ratio < 1) {
    m = vert_richness_mod(f1, dat = d, family = gaussian(), edge.ratio = edge.ratio)
    san = sanity(m)$all_ok
    edge.ratio = edge.ratio + 0.05
    print(edge.ratio)
  }
  
  b = gsub(" ", "_", b)
  b = gsub("/", "_", b)
  saveRDS(m, paste0("results/sdmTMB_models/richness_vert_models/sesvert/biomes/amphibians/", b, ".rds"))
}

for(b in unique(amph$biome)) {
  print(b)
  d = amph %>% filter(biome == b)
  print(nrow(d))
}



for(b in unique(birds$biome)) {
  d = birds %>% filter(biome == b)
  san = F
  edge.ratio = 0.25
  
  fname = gsub(" ", "_", b)
  fname = gsub("/", "_", fname)
  fname = paste0("results/sdmTMB_models/richness_vert_models/sesvert/biomes/birds/", fname, ".rds")
  if(!file.exists(fname) & b != "Temperate Grasslands, Savannas & Shrublands" & b != "Montane Grasslands & Shrublands" & 
     b != "Mediterranean Forests, Woodlands & Scrub" & nrow(d) >= 30) {
    while(!san & edge.ratio < 1) {
      m = vert_richness_mod(f1, dat = d, family = gaussian(), edge.ratio = edge.ratio)
      san = sanity(m)$all_ok
      edge.ratio = edge.ratio + 0.05
      print(edge.ratio)
    }
    saveRDS(m, fname)
  }
  
}

for(b in unique(birds$biome)) {
  print(b)
  d = birds %>% filter(biome == b)
  print(nrow(d))
}

# mammals
for(b in unique(mammals$biome)) {
  d = mammals %>% filter(biome == b)
  san = F
  edge.ratio = 0.2
  
  fname = gsub(" ", "_", b)
  fname = gsub("/", "_", fname)
  fname = paste0("results/sdmTMB_models/richness_vert_models/sesvert/biomes/mammals/", fname, ".rds")
  if(!file.exists(fname) & nrow(d) >= 30) {
    while(!san & edge.ratio < 1) {
      m = vert_richness_mod(f1, dat = d, family = gaussian(), edge.ratio = edge.ratio)
      san = sanity(m)$all_ok
      edge.ratio = edge.ratio + 0.05
      print(edge.ratio)
    }
    saveRDS(m, fname)
  }
}

for(b in unique(mammals$biome)) {
  print(b)
  d = mammals %>% filter(biome == b)
  print(nrow(d))
}

# reptiles
for(b in unique(rept$biome)) {
  
  d = rept %>% filter(biome == b)
  san = F
  edge.ratio = 0.2
  
  fname = gsub(" ", "_", b)
  fname = gsub("/", "_", fname)
  fname = paste0("results/sdmTMB_models/richness_vert_models/sesvert/biomes/reptiles/", fname, ".rds")
  
  if(!file.exists(fname) & nrow(d) >= 30) {
    while(!san & (edge.ratio < 1)) {
      m = vert_richness_mod(f1, dat = d, family = gaussian(), edge.ratio = edge.ratio)
      san = sanity(m)$all_ok
      edge.ratio = edge.ratio + 0.05
      print(edge.ratio)
    }
    #saveRDS(m, fname)
  }
}


for(b in unique(reptiles$biome)) {
  print(b)
  d = reptiles %>% filter(biome == b)
  print(nrow(d))
}








