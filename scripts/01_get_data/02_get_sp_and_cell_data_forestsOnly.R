# Make species level and gridcell level dataframes for each class
library(tidyverse)
library(rnaturalearth)
source("scripts/00_functions/00_functions.R")
library(terra)
library(data.table)
library(tictoc)
library(foreach)
library(doParallel)

### HARD CODES ###
# phylogenetic imputation threshold for binary diel activity and verticality variables
# greater than threshold = 1
# less than or equal to threshold = 0
impute_thresh = 0.7

# extent of occurrence layer
mapa <- rast(xmin = -20592508, xmax = 20588492, ymin = -5743602, ymax = 6573398,
             crs = "+proj=cea +datum=WGS84")
res(mapa) = 111000

# terrestrial outline and gridcells
# world = rnaturalearth::ne_countries()
# world = vect(world)
# world = as.polygons(world)
# world = aggregate(world, dissolve = T)
# writeVector(world, "data/original/rnaturalearth_world.shp")
world = vect("data/original/rnaturalearth_world.shp")
values(world) <- 1
world.cea = project(world, "+proj=cea +datum=WGS84")
world.rast = rasterize(world, mapa)
values(world.rast) <- 1
world.rast = mask(world.rast, world.cea)
world.df = as.data.frame(world.rast, xy = T) %>% 
  dplyr::select(x,y)

# load env data
env = read.csv("data/derivative_data/env_data.csv")

# clip env data to wooded habitats only
env.forest = env %>% 
  filter(grepl("forest|woodland|taiga|Forest|várzea|Yungas|savanna|thicket|mallee", ecoregion))

# Moura traits
traits = read.csv("data/original/trait_data/TetrapodTraits_1.0.0/TetrapodTraits_1.0.0.csv")

traits %>% 
  filter(Palearctic == 1) %>% 
  ggplot(aes(x = factor(Verticality), y = RangeSize)) +
  geom_point()+
  geom_boxplot() +
  coord_cartesian(ylim = c(0,500))

plot(traits$Verticality, traits$RangeSize)

colnames(traits)

traits = traits %>%
  dplyr::select(Scientific.Name:Class, BodyMass_g, Diu:Noc, Nocturnality, Fos:Aer, Verticality,
                EcoTer:EcoMar, RangeSize:Elevation)

sum(is.na(traits$EcoMar))

# drop marine only species
# drop marine + freshwater species
traits = traits %>% 
  filter(!(EcoMar == 1 & EcoTer == 0 & EcoFresh == 0)) %>% 
  filter(!(EcoMar == 1 & EcoTer == 0 & EcoFresh == 1))

# drop aerial only species
traits = traits %>% 
  filter(!(Aer == 1 & Arb == 0 & Ter == 0 & Aqu == 0 & Fos == 0))

# convert imputed diel activity and verticality to binary variables using a threshold
traits = traits %>% 
  mutate(Noc = ifelse(Noc > impute_thresh, 1, 0),
         Diu = ifelse(Diu > impute_thresh, 1, 0),
         Ter = ifelse(Ter > impute_thresh, 1, 0),
         Fos = ifelse(Fos > impute_thresh, 1, 0),
         Aqu = ifelse(Fos > impute_thresh, 1, 0),
         Arb = ifelse(Arb > impute_thresh, 1, 0),
         Aer = ifelse(Aer > impute_thresh, 1, 0))

# remove species that are categorized as 0 for all categories in diel activity or verticality
traits = traits %>% 
  filter(!(Noc == 0 & Diu == 0)) %>% 
  filter(!(Ter == 0 & Fos == 0 & Aqu == 0 & Arb == 0 & Aer == 0))

# test %>% filter(Noc == 0 & Diu == 0) 
# nrow(test %>% filter(Fos == 0 & Ter == 0 & Arb == 0 & Aqu == 0))

traits = traits %>% 
  rename(sciname = Scientific.Name,
         body_size = BodyMass_g) %>% 
  mutate(sciname = gsub(" ","_", sciname))

# Amphibians --------------------------------------------------------------

# load amphibian occurrence data (cover0 includes all species, even if ranges cover less than 1/3 grid cell)
occr <- fread("data/derivative_data/pres_abs_cover0/pres_abs_amph.csv")
colnames(occr)[1:2] = c("x", "y")
colnames(occr) = gsub(" ", "_", colnames(occr))

# list of species from IUCN range maps
spp <- colnames(occr)[-1:-2]
spp <- data.frame(sciname = gsub(" ","_",spp))

# get species and traits in common between iucn and trait dataset
traits_amph = inner_join(spp, traits, by = "sciname")

# subset occurrence data to species that occur in both iucn and trait datasets
occr = occr %>% 
  dplyr::select(x,y,traits_amph$sciname)

# filter out rows that don't have realm information
occr = env.forest %>% 
  dplyr::select(x,y, realm) %>% 
  left_join(occr, by = c("x", "y")) %>% 
  drop_na(realm) %>% 
  dplyr::select(!realm)

# remove species that don't occur in any cells
spocc = apply(occr, 2, sum)
sppres = which(spocc > 0)
occr = occr[,sppres]

# subset traits to species that occur in the occr dataframe
traits_amph = traits_amph %>% 
  filter(sciname %in% colnames(occr))

# test to make sure species in trait and occr dataframes are in the same order
sum(traits_amph$sciname != colnames(occr[3:ncol(occr)]))

# * - species level data ------------------------------------------------------
amph_spdat = get_sp_dat(trait_dat = traits_amph, occ = occr, env = env)
write.csv(amph_spdat, "data/derivative_data/species_data/amph_spdat_moura_forestsOnly.csv", row.names = F)

amph_spdat = read.csv("data/derivative_data/species_data/amph_spdat_moura_forestsOnly.csv")

# takes 13 min

tic()
amph_comdat = get_gridcell_vert_parallel(trait_dat = amph_spdat, occ = occr, env = env, rich_min = 5, ncore = 7, nsim = 100)
write.csv(amph_comdat, "data/derivative_data/gridcell_data/amphibians_comdat/amph_comdat_parallel_forestsOnly.csv", row.names = F)
toc()


# Mammals -----------------------------------------------------------------

# load breeding bird occurrence data (cover0 includes all species, even if ranges cover less than 1/3 grid cell)
occr <- fread("data/derivative_data/pres_abs_cover0/pres_abs_mammals.csv")
colnames(occr)[1:2] = c("x", "y")
colnames(occr) = gsub(" ", "_", colnames(occr))

# list of species from IUCN range maps
spp <- colnames(occr)[-1:-2]
spp <- data.frame(sciname = gsub(" ","_",spp))

# get species and traits in common between iucn and trait dataset
traits_mammals = inner_join(spp, traits, by = "sciname")

# subset occurrence data to species that occur in both iucn and trait datasets
occr = occr %>% 
  dplyr::select(x,y,traits_mammals$sciname)

# filter out rows that don't have realm information
occr = env.forest %>% 
  dplyr::select(x,y, realm) %>% 
  left_join(occr, by = c("x", "y")) %>% 
  drop_na(realm) %>% 
  dplyr::select(!realm)

# remove species that don't occur in any cells
spocc = apply(occr, 2, sum)
sppres = which(spocc > 0)
occr = occr[,sppres]

# subset traits to species that occur in the occr dataframe
traits_mammals = traits_mammals %>% 
  filter(sciname %in% colnames(occr))

# test to make sure species in trait and occr dataframes are in the same order
sum(traits_mammals$sciname != colnames(occr[3:ncol(occr)]))

# * - species level data ------------------------------------------------------
mammals_spdat = get_sp_dat(trait_dat = traits_mammals, occ = occr, env = env, filter_col = NA)
write.csv(mammals_spdat, "data/derivative_data/species_data/mammals_spdat_moura_forestsOnly.csv", row.names = F)

mammals_spdat = read.csv("data/derivative_data/species_data/mammals_spdat_moura_forestsOnly.csv")

# takes 1210 sec
tic()
mammals_comdat = get_gridcell_vert_parallel(trait_dat = mammals_spdat, occ = occr, env = env, rich_min = 5, ncore = 7, nsim = 100)
write.csv(mammals_comdat, "data/derivative_data/gridcell_data/mammals_comdat/mammals_comdat_parallel_forestsOnly.csv", row.names = F)
toc()

# Reptiles --------------------------------------------------------------

# load reptile occurrence data (cover0 includes all species, even if ranges cover less than 1/3 grid cell)
occr <- fread("data/derivative_data/pres_abs_cover0/pres_abs_rept.csv")
colnames(occr)[1:2] = c("x", "y")
colnames(occr) = gsub(" ", "_", colnames(occr))

# list of species from IUCN range maps
spp <- colnames(occr)[-1:-2]
spp <- data.frame(sciname = gsub(" ","_",spp))

# get species and traits in common between iucn and trait dataset
traits_rept = inner_join(spp, traits, by = "sciname")

# subset occurrence data to species that occur in both iucn and trait datasets
occr = occr %>% 
  dplyr::select(x,y,traits_rept$sciname)

# filter out rows that don't have realm information
occr = env.forest %>% 
  dplyr::select(x,y, realm) %>% 
  left_join(occr, by = c("x", "y")) %>% 
  drop_na(realm) %>% 
  dplyr::select(!realm)

# remove species that don't occur in any cells
spocc = apply(occr, 2, sum)
sppres = which(spocc > 0)
occr = occr[,sppres]

# subset traits to species that occur in the occr dataframe
traits_rept = traits_rept %>% 
  filter(sciname %in% colnames(occr))

# test to make sure species in trait and occr dataframes are in the same order
sum(traits_rept$sciname != colnames(occr[3:ncol(occr)]))

# * - species level data ------------------------------------------------------
rept_spdat = get_sp_dat(trait_dat = traits_rept, occ = occr, env = env, filter_col = NA)
write.csv(rept_spdat, "data/derivative_data/species_data/rept_spdat_moura_forestsOnly.csv", row.names = F)

rept_spdat = read.csv("data/derivative_data/species_data/rept_spdat_moura_forestsOnly.csv")

# takes 2825 sec
tic()
rept_comdat = get_gridcell_vert_parallel(trait_dat = rept_spdat, occ = occr, env = env, rich_min = 5, ncore = 7, nsim = 100)
write.csv(rept_comdat, "data/derivative_data/gridcell_data/reptiles_comdat/rept_comdat_parallel_forestsOnly.csv", row.names = F)
toc()



# Birds -------------------------------------------------------------------

# load reptile occurrence data (cover0 includes all species, even if ranges cover less than 1/3 grid cell)
occr <- fread("data/derivative_data/pres_abs_cover0/pres_abs_birds_breeding_resident.csv")
colnames(occr)[1:2] = c("x", "y")
colnames(occr) = gsub(" ", "_", colnames(occr))

# list of species from IUCN range maps
spp <- colnames(occr)[-1:-2]
spp <- data.frame(sciname = gsub(" ","_",spp))

# Elton trait data

traits.birds <- read.csv("data/original/trait_data/elton_traits/BirdFuncDat.txt", sep = "\t")
names(traits.birds)

# Get traits of interest - ForStrat-Value, BodyMass-Value
traits.birds = trait.birds %>% 
  dplyr::select(PassNonPass:BLFamilyEnglish, Scientific:PelagicSpecialist, Nocturnal:BodyMass.Value) %>% 
  rename(body_size = BodyMass.Value,
         sciname = Scientific,
         family = BLFamilyLatin) %>% 
  mutate(sciname = gsub(" ","_", sciname)) %>% 
  drop_na()

names(trait)

# calculate verticality score based on foraging strata
traits.birds = traits.birds %>% 
  mutate(Verticality = (ForStrat.watbelowsurf * 0.5 + ForStrat.wataroundsurf * 0.5 +
                                ForStrat.ground * 0.5 + ForStrat.understory * 0.667 + ForStrat.midhigh * 0.8337 +
                                ForStrat.canopy * 1 + ForStrat.aerial * 1)/100)

# get species and traits in common between iucn and trait dataset
traits.birds = inner_join(spp, traits.birds, by = "sciname")

# 7914 species in iucn database with body mass values from Elton traits

traits.birds.moura = traits %>% 
  dplyr::select(sciname, EcoMar, EcoTer, EcoFresh)
# drop marine only species
# drop marine + freshwater species
trait = trait %>% 
  left_join(traits.birds.moura, by = "sciname") %>% 
  filter(!(EcoMar == 1 & EcoTer == 0 & EcoFresh == 0)) %>% 
  filter(!(EcoMar == 1 & EcoTer == 0 & EcoFresh == 1))

# drop aerial only species
trait = trait %>% 
  filter(!(ForStrat.aerial == 1 & ForStrat.canopy == 0 & ForStrat.midhigh == 0 & ForStrat.understory == 0 & ForStrat.ground == 0 & ForStrat.wataroundsurf == 0 & ForStrat.watbelowsurf == 0))

# leaves 6907 species

# subset occurrence data to species that occur in both iucn and trait datasets
occr = occr %>% 
  dplyr::select(x,y,trait$sciname)

# filter out rows that don't have realm information
occr = env.forest %>% 
  dplyr::select(x,y, realm) %>% 
  left_join(occr, by = c("x", "y")) %>% 
  drop_na(realm) %>% 
  dplyr::select(!realm)

# remove species that don't occur in any cells
spocc = apply(occr, 2, sum)
sppres = which(spocc > 0)
occr = occr[,sppres]

# subset traits to species that occur in the occr dataframe
trait = trait %>% 
  filter(sciname %in% colnames(occr))

# test to make sure species in trait and occr dataframes are in the same order
sum(trait$sciname != colnames(occr[3:ncol(occr)]))

# * - species level data ------------------------------------------------------
bird_spdat = get_sp_dat(trait_dat = trait, occ = occr, env = env, filter_col = NA, traitbase = "elton")
write.csv(bird_spdat, "data/derivative_data/species_data/birds_spdat_elton_forestsOnly.csv", row.names = F)

bird_spdat = read.csv("data/derivative_data/species_data/birds_spdat_elton_forestsOnly.csv")

# takes 2825 sec
tic()
bird_comdat = get_gridcell_vert_parallel(trait_dat = bird_spdat, occ = occr, env = env, rich_min = 5, ncore = 7, nsim = 100)
write.csv(bird_comdat, "data/derivative_data/gridcell_data/birds_comdat/birds_comdat_parallel_elton_forestsOnly.csv", row.names = F)
toc()



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



