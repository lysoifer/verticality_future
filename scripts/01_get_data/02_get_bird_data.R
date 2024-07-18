# Lydia Soifer
# May 1, 2024
# Get arboreality and environmental data and add to the matrix from get_pres_abs_data

library(terra)
library(data.table)
library(tictoc)
library(chelsaDL)
library(tidyverse)
library(doParallel)
library(parallel)
library(foreach)
library(mFD)
source("scripts_lydia/00_functions.R")

# REFERENCE PROJECTION MAP
mapa <- rast(xmin = -20592508, xmax = 20588492, ymin = -5743602, ymax = 6573398,
             crs = "+proj=cea +datum=WGS84")
res(mapa) <- 111000

mundi = vect("data/Mundi_contour/Mundi_contour.shp")
mundi = project(mundi, "+proj=cea +datum=WGS84")
mundi.rast = rasterize(mundi, mapa)
mundi.df = as.data.frame(mundi.rast, xy = T) %>% 
  dplyr::select(x,y)

# load breeding bird occurrence data (cover0 includes all species, even if ranges cover less than 1/3 grid cell)
occr <- fread("data/derivative_data/pres_abs_cover0/pres_abs_birds_breeding_resident.csv")
colnames(occr)[1:2] = c("x", "y")
colnames(occr) = gsub(" ", "_", colnames(occr))

# remove ocean pixels
occr = right_join(occr, mundi.df, by = c("x", "y"))

# list of species from IUCN range maps
spp <- colnames(occr)[-1:-2]
spp <- data.frame(sciname = gsub(" ","_",spp))

# load env data
env = read.csv("data/derivative_data/env_data.csv")

# Trait and Verticality data --------------------------------------------------------------

# Using Elton traits  
# Hamish Wilman, Jonathan Belmaker, Jennifer Simpson, Carolina de la Rosa, 
# Marcelo M. Rivadeneira, and Walter Jetz. 2014. EltonTraits 1.0: Species-level 
# foraging attributes of the world's birds and mammals. 
# Ecology 95:2027. http://dx.doi.org/10.1890/13-1917.1

### OPEN LARGE DATASET
trait <- read.csv("data/original/trait_data/elton_traits/BirdFuncDat.txt", sep = "\t")

names(trait)

# Get traits of interest - ForStrat-Value, BodyMass-Value
trait = trait %>% 
  dplyr::select(PassNonPass:BLFamilyEnglish, Scientific:PelagicSpecialist, Nocturnal:BodyMass.Value) %>% 
  rename(body_size = BodyMass.Value,
         sciname = Scientific,
         family = BLFamilyLatin) %>% 
  mutate(sciname = gsub(" ","_", sciname)) %>% 
  drop_na()

names(trait)

# add verticality scores to trait dataset
vert = read.csv("data/original/Moura_arboreality_vertebrates.csv")

trait = vert %>% 
  rename(sciname = "Scientific.Name") %>% 
  mutate(sciname = gsub(" ", "_", sciname)) %>% 
  dplyr::select(sciname, Order, Fos:Aer, Verticality) %>% 
  inner_join(trait, by = "sciname")

# what would verticality by elton traits look like in comparison to Moura traits
trait = trait %>%
  mutate(verticality.elton = (ForStrat.watbelowsurf * 0.5 + ForStrat.wataroundsurf * 0.5 +
           ForStrat.ground * 0.5 + ForStrat.understory * 0.667 + ForStrat.midhigh * 0.8337 +
           ForStrat.canopy * 1 + ForStrat.aerial * 1)/100)

plot(trait$verticality.elton, trait$Verticality)

# remove aerial only species
trait = trait %>% 
  mutate(aer_only = case_when((Fos + Ter + Aqu + Arb + Aer)==1 & Aer == 1 ~ 1))

View(trait %>% filter(aer_only == 1))
unique(trait %>% filter(aer_only == 1) %>% dplyr::select(BLFamilyEnglish))

# aerial only species include species in swifts, woodswallows, petrels and shearwaters
# nightjars, tyrant-flycatchers, swallows and martins, rollers, treeswifts,
# osprey/kites/hawks/eagles, bee-eaters, butcherbirds, and coursers and pratincoles

# remove the 141 aerial only species
trait = trait %>% 
  filter(is.na(aer_only)) %>% 
  dplyr::select(!aer_only)

# get species and traits in common between iucn and trait dataset
trait = inner_join(spp, trait, by = "sciname")

# 7811 species in iucn database with body mass values from Elton traits

# subset occurrence data to species that occur in both iucn and trait datasets
occr = occr %>% 
  dplyr::select(x,y,trait$sciname)

# species level data ------------------------------------------------------
birds_breed_spdat = get_sp_dat(trait_dat = trait, occ = occr, env = env)
write.csv(birds_breed_spdat, "data/derivative_data/species_data/birds_breedresident_spdat.csv", row.names = F)


# community level data per grid cell --------------------------------------
birds_breed_spdat = read.csv("data/derivative_data/species_data/birds_breedresident_spdat.csv")

birds_breed_comdat = get_gridcell_dat(trait_dat = birds_breed_spdat, occ = occr, env = env, rich_min = 5)
write.csv(birds_breed_comdat, "data/derivative_data/gridcell_data/birds_comdat/birds_breedresident_comdat.csv")



# nonbreeding birds -------------------------------------------------------

# load breeding bird occurrence data (cover0 includes all species, even if ranges cover less than 1/3 grid cell)
occr <- fread("data/derivative_data/pres_abs_cover0/pres_abs_birds_nonbreeding_resident.csv")
colnames(occr)[1:2] = c("x", "y")
colnames(occr) = gsub(" ", "_", colnames(occr))

# remove ocean pixels
occr = right_join(occr, mundi.df, by = c("x", "y"))

# list of species from IUCN range maps
spp <- colnames(occr)[-1:-2]
spp <- data.frame(sciname = gsub(" ","_",spp))

# load env data
env = read.csv("data/derivative_data/env_data.csv")

# Trait and Verticality data --------------------------------------------------------------

# Using Elton traits  
# Hamish Wilman, Jonathan Belmaker, Jennifer Simpson, Carolina de la Rosa, 
# Marcelo M. Rivadeneira, and Walter Jetz. 2014. EltonTraits 1.0: Species-level 
# foraging attributes of the world's birds and mammals. 
# Ecology 95:2027. http://dx.doi.org/10.1890/13-1917.1

### OPEN LARGE DATASET
trait <- read.csv("data/original/trait_data/elton_traits/BirdFuncDat.txt", sep = "\t")

names(trait)

# Get traits of interest - ForStrat-Value, BodyMass-Value
trait = trait %>% 
  dplyr::select(PassNonPass:BLFamilyEnglish, Scientific:PelagicSpecialist, Nocturnal:BodyMass.Value) %>% 
  rename(body_size = BodyMass.Value,
         sciname = Scientific,
         family = BLFamilyLatin) %>% 
  mutate(sciname = gsub(" ","_", sciname)) %>% 
  drop_na()

names(trait)

# add verticality scores to trait dataset
vert = read.csv("data/original/Moura_arboreality_vertebrates.csv")

trait = vert %>% 
  rename(sciname = "Scientific.Name") %>% 
  mutate(sciname = gsub(" ", "_", sciname)) %>% 
  dplyr::select(sciname, Order, Fos:Aer, Verticality) %>% 
  inner_join(trait, by = "sciname")

# what would verticality by elton traits look like in comparison to Moura traits
trait = trait %>% 
  mutate(verticality.elton = (ForStrat.watbelowsurf * 0.5 + ForStrat.wataroundsurf * 0.5 +
                                ForStrat.ground * 0.5 + ForStrat.understory * 0.667 + ForStrat.midhigh * 0.8337 +
                                ForStrat.canopy * 1 + ForStrat.aerial * 1)/100)

plot(trait$verticality.elton, trait$Verticality)

# get species and traits in common between iucn and trait dataset
trait = inner_join(spp, trait, by = "sciname")

# 7914 species in iucn database with body mass values from Elton traits

# subset occurrence data to species that occur in both iucn and trait datasets
occr = occr %>% 
  dplyr::select(x,y,trait$sciname)

# species level data ------------------------------------------------------
birds_nonbreed_spdat = get_sp_dat(trait_dat = trait, occ = occr, env = env)
write.csv(birds_nonbreed_spdat, "data/derivative_data/species_data/birds_nonbreedresident_spdat.csv", row.names = F)


# community level data per grid cell --------------------------------------

birds_nonbreed_comdat = get_gridcell_dat(trait_dat = birds_nonbreed_spdat, occ = occr, env = env, rich_min = 5)
write.csv(birds_nonbreed_comdat, "data/derivative_data/gridcell_data/birds_comdat/birds_nonbreedresident_comdat.csv")


# check correlation between breeding and nonbreeding ----------------------

breed = read.csv("data/derivative_data/gridcell_data/birds_comdat/birds_breedresident_comdat.csv", row.names = "X")
nonbreed = read.csv("data/derivative_data/gridcell_data/birds_comdat/birds_nonbreedresident_comdat.csv", row.names = "X")

colnames(breed)[3:ncol(breed)] = paste0("breed_", colnames(breed)[3:ncol(breed)])
colnames(nonbreed)[3:ncol(nonbreed)] = paste0("nonbreed_", colnames(nonbreed)[3:ncol(nonbreed)])

bb = inner_join(breed, nonbreed, by = c("x", "y"))

plot(bb$breed_richness, bb$nonbreed_richness, 
     xlab = "richness (breeding + resident)", ylab = "richness (non-breeding + resident)")
cor.test(bb$breed_richness, bb$nonbreed_richness)

plot(bb$breed_body.size, bb$nonbreed_body.size, 
     xlab = "body size (breeding + resident)", ylab = "body size (non-breeding + resident)")
cor.test(bb$breed_body.size, bb$nonbreed_body.size)

plot(bb$breed_vert, bb$nonbreed_vert, 
     xlab = "mean verticality (breeding + resident)", ylab = "mean vert (non-breeding + resident)")
cor(bb$breed_vert, bb$nonbreed_vert)

plot(bb$breed_p.arb, bb$nonbreed_p.arb, 
     xlab = "% arboreal (breeding + resident)", ylab = "% arboreal (non-breeding + resident)")
cor(bb$breed_p.arb, bb$nonbreed_p.arb)


# breeding + elton verticality --------------------------------------------

# load breeding bird occurrence data (cover0 includes all species, even if ranges cover less than 1/3 grid cell)
occr <- fread("data/derivative_data/pres_abs_cover0/pres_abs_birds_breeding_resident.csv")
colnames(occr)[1:2] = c("x", "y")
colnames(occr) = gsub(" ", "_", colnames(occr))

# remove ocean pixels
#occr = right_join(occr, mundi.df, by = c("x", "y"))

# list of species from IUCN range maps
spp <- colnames(occr)[-1:-2]
spp <- data.frame(sciname = gsub(" ","_",spp))

# load env data
env = read.csv("data/derivative_data/env_data.csv")

# Trait and Verticality data --------------------------------------------------------------

# Using Elton traits  
# Hamish Wilman, Jonathan Belmaker, Jennifer Simpson, Carolina de la Rosa, 
# Marcelo M. Rivadeneira, and Walter Jetz. 2014. EltonTraits 1.0: Species-level 
# foraging attributes of the world's birds and mammals. 
# Ecology 95:2027. http://dx.doi.org/10.1890/13-1917.1

### OPEN LARGE DATASET
trait <- read.csv("data/original/trait_data/elton_traits/BirdFuncDat.txt", sep = "\t")
trait.moura = read.csv("data/original/trait_data/TetrapodTraits_1.0.0/TetrapodTraits_1.0.0.csv")

names(trait)

# Get traits of interest - ForStrat-Value, BodyMass-Value
trait = trait %>% 
  dplyr::select(PassNonPass:BLFamilyEnglish, Scientific:PelagicSpecialist, Nocturnal:BodyMass.Value) %>% 
  rename(body_size = BodyMass.Value,
         sciname = Scientific,
         family = BLFamilyLatin) %>% 
  mutate(sciname = gsub(" ","_", sciname)) %>% 
  drop_na()

names(trait)

# add verticality scores to trait dataset
vert = read.csv("data/original/Moura_arboreality_vertebrates.csv")

trait = vert %>% 
  rename(sciname = "Scientific.Name") %>% 
  mutate(sciname = gsub(" ", "_", sciname)) %>% 
  dplyr::select(sciname, Order, Fos:Aer, Verticality) %>% 
  inner_join(trait, by = "sciname")

# what would verticality by elton traits look like in comparison to Moura traits
trait = trait %>% 
  mutate(verticality.elton = (ForStrat.watbelowsurf * 0.5 + ForStrat.wataroundsurf * 0.5 +
                                ForStrat.ground * 0.5 + ForStrat.understory * 0.667 + ForStrat.midhigh * 0.8337 +
                                ForStrat.canopy * 1 + ForStrat.aerial * 1)/100)

plot(trait$verticality.elton, trait$Verticality)

trait = trait %>% 
  dplyr::select(!Verticality) %>% 
  rename(Verticality = verticality.elton)

# get species and traits in common between iucn and trait dataset
trait = inner_join(spp, trait, by = "sciname")

# 7914 species in iucn database with body mass values from Elton traits

# subset occurrence data to species that occur in both iucn and trait datasets
occr = occr %>% 
  dplyr::select(x,y,trait$sciname)

# species level data ------------------------------------------------------
birds_breed_spdat = get_sp_dat(trait_dat = trait, occ = occr, env = env)
write.csv(birds_breed_spdat, "data/derivative_data/species_data/birds_breedresident_spdat_eltonvert.csv", row.names = F)



# community level data per grid cell --------------------------------------

birds_breed_spdat = read.csv("data/derivative_data/species_data/birds_breedresident_spdat_eltonvert.csv")

# join nocturnal/diurnal from moura dataset
trait.moura = trait.moura %>% 
  rename(sciname = Scientific.Name) %>% 
  mutate(sciname = gsub(" ", "_", sciname)) %>% 
  dplyr::select(sciname, Arb, Fos, Ter, Aqu, Aer, Diu, Noc, Nocturnality, EcoTer, EcoMar, EcoFresh)

birds_breed_spdat = birds_breed_spdat %>% 
  dplyr::select(!c(Arb, Fos, Ter, Aqu, Aer)) %>% 
  left_join(trait.moura, by = "sciname") %>% 
  filter(!(EcoMar == 1 & EcoTer == 0 & EcoFresh == 0)) %>% 
  filter(!(EcoMar == 1 & EcoTer == 0 & EcoFresh == 1)) %>% 
  filter(!(Aer == 1 & Arb == 0 & Ter == 0 & Aqu == 0 & Fos == 0))

# test data
# subset occurrence data to species that occur in both iucn and trait datasets
occr = occr %>% 
  dplyr::select(x,y,birds_breed_spdat$sciname)

# filter out rows that don't have realm information
occr = env %>% 
  dplyr::select(x,y, realm) %>% 
  left_join(occr, by = c("x", "y")) %>% 
  drop_na(realm) %>% 
  dplyr::select(!realm)

# remove species that don't occur in any cells
spocc = apply(occr, 2, sum)
sppres = which(spocc > 0)
occr = occr[,sppres]

# subset traits to species that occur in the occr dataframe
birds_breed_spdat = birds_breed_spdat %>% 
  filter(sciname %in% colnames(occr))

# test to make sure species in trait and occr dataframes are in the same order
sum(birds_breed_spdat$sciname != colnames(occr[3:ncol(occr)]))

# test = occr[100:200,]
# s = apply(test, 2, sum)
# test = test[,-which(s==0)]
# 
# test_traits = birds_breed_spdat %>% filter(sciname %in% colnames(test[3:ncol(test)]))
# 
# testrun = get_gridcell_dat_parallel(trait_dat = test_traits, occ = test, env = env, rich_min = 5, ncore = 7, nsim = 3, eltonbirds = T)

birds_breed_comdat = get_gridcell_dat_parallel(trait_dat = birds_breed_spdat, occ = occr, env = env, rich_min = 5, ncore = 7, nsim = 100, eltonbirds = T)
write.csv(birds_breed_comdat, "data/derivative_data/gridcell_data/birds_comdat/birds_breedresident_comdat_eltonvert.csv")

# 6 hours
tic()
birds_trait_volume = get_trait_volume(trait_dat = birds_breed_spdat, occ = occr, env = env, rich_min = 5, ncore = 7, nsim = 100, eltonbirds = TRUE)
write.csv(birds_trait_volume, "data/derivative_data/gridcell_data/birds_comdat/birds_comdat_eltonvert_volume.csv", row.names = F)
toc()

birds_elton = fread("data/derivative_data/gridcell_data/birds_comdat/birds_breedresident_comdat_eltonvert.csv")
birds_elton = left_join(birds_elton, birds_trait_volume, by = c("x", "y"))

write.csv(birds_elton, "data/derivative_data/gridcell_data/birds_comdat/birds_breedresident_comdat_elton.csv", row.names = F)


# Moura vs Elton verticality score results --------------------------------

# for breeding birds
moura = read.csv("data/derivative_data/gridcell_data/birds_comdat/birds_breedresident_comdat.csv")
elton = read.csv("data/derivative_data/gridcell_data/birds_comdat/birds_breedresident_comdat_eltonvert.csv")

sum(moura$x != elton$x)
sum(moura$y != elton$y)

plot(moura$ses.vert, elton$ses.vert)
abline(a = 0, b=1, col = "red4", lwd = 4)
plot(moura$vert, elton$vert)
abline(a = 0, b=1, col = "red4", lwd = 4)



# diurnal vs noctural comdat ----------------------------------------------

# load breeding bird occurrence data (cover0 includes all species, even if ranges cover less than 1/3 grid cell)
occr <- fread("data/derivative_data/pres_abs_cover0/pres_abs_birds_breeding_resident.csv")
colnames(occr)[1:2] = c("x", "y")
colnames(occr) = gsub(" ", "_", colnames(occr))

# remove ocean pixels
occr = right_join(occr, mundi.df, by = c("x", "y"))

# list of species from IUCN range maps
spp <- colnames(occr)[-1:-2]
spp <- data.frame(sciname = gsub(" ","_",spp))

# load env data
env = read.csv("data/derivative_data/env_data.csv")

# load species data
birds_breed_spdat = read.csv("data/derivative_data/species_data/birds_breedresident_spdat.csv")

occr = occr %>% 
  dplyr::select(x,y,birds_breed_spdat$sciname)

#  * - community level data per grid cell --------------------------------------
birds_breed_noc_comdat = get_gridcell_dat(trait_dat = birds_breed_spdat, occ = occr, env = env, filter_col = "Nocturnal", filter_val = 1, rich_min = 5)
write.csv(birds_breed_noc_comdat, "data/derivative_data/gridcell_data/birds_comdat/birds_breedresident_noc_comdat.csv")

birds_breed_diu_comdat = get_gridcell_dat(trait_dat = birds_breed_spdat, occ = occr, env = env, filter_col = "Nocturnal", filter_val = 0, rich_min = 5)
write.csv(birds_breed_diu_comdat, "data/derivative_data/gridcell_data/birds_comdat/birds_breedresident_diu_comdat.csv")

