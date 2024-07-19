# Make species level and gridcell level dataframes for each class
# Only included grid cells with wooded habitat
# 50 km resolution

library(tidyverse)
library(rnaturalearth)
source("scripts/00_functions/00_functions.R")
source("scripts/00_functions/00_functions2.R")
source("scripts/00_functions/manuscript_functions.R")
library(terra)
library(data.table)
library(tictoc)
library(foreach)
library(doParallel)
library(picante)

### HARD CODES ###
# phylogenetic imputation threshold for binary diel activity and verticality variables
# greater than threshold = 1
# less than or equal to threshold = 0
impute_thresh = 0.7

# extent of occurrence layer
mapa <- rast(xmin = -20592508, xmax = 20588492, ymin = -5743602, ymax = 6573398,
             crs = "+proj=cea +datum=WGS84")
#res(mapa) = 111000
res(mapa) = 50000

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
env = read.csv("data/derivative_data/env_data_50km.csv")

# clip env data to wooded habitats only
env.forest = env %>% 
  filter(grepl("forest|woodland|taiga|Forest|várzea|Yungas|savanna|thicket|mallee", ecoregion))
fwrite(env.forest, "data/derivative_data/env_data_50km_forest.csv", row.names = F)

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
#occr <- fread("data/derivative_data/pres_abs_cover0/pres_abs_amph.csv")
occr = fread("data/derivative_data/pres_abs_cover0/resolution_50km/pres_abs_amph.csv")
colnames(occr)[1:2] = c("x", "y")
colnames(occr) = gsub(" ", "_", colnames(occr))

# list of species from IUCN range maps
spp <- colnames(occr)[-1:-2]
spp <- data.frame(sciname = gsub(" ","_",spp))

# load phylogeny
#tree = read.tree("data/original/phylotrees/amph_shl_new_Consensus_7238.tre")
#tree.sp = data.frame(sciname = tree$tip.label)

# get species and traits in common between iucn and trait dataset
# and between iucn and phylogenetic tree
traits_amph = inner_join(spp, traits, by = "sciname")
#traits_amph = inner_join(traits_amph, tree.sp, by = "sciname")

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
amph_spdat = get_sp_dat(trait_dat = traits_amph, occ = occr, env = env, traitbase = "moura")
write.csv(amph_spdat, "data/derivative_data/species_data/amph_spdat_moura_forestsOnly.csv", row.names = F)

amph_spdat = read.csv("data/derivative_data/species_data/amph_spdat_moura_forestsOnly.csv")

# takes 7.8 minutes
tic()
amph_comdat = gridcell_dat(occ = occr, tr = amph_spdat, env = env.forest, reps = 100)
toc()

griddat = lapply(amph_comdat, "[[", 1)
griddat = rbindlist(griddat)

# alphabetize occurrence dataframe by species
occdat = lapply(amph_comdat, "[[", 2)
occdat = rbindlist(occdat, use.names = TRUE, fill = TRUE)
x.index = which(colnames(occdat) == "x")
y.index = which(colnames(occdat) == "y")
spnames = colnames(occdat)[-c(x.index, y.index)]
occdat = setcolorder(occdat, c('x', 'y', sort(spnames)))

# replace NAs with zeros
occdat[is.na(occdat)] <- 0

fwrite(griddat, "data/derivative_data/gridcell_data/env_forest/50_km/amph_comdat.csv", row.names = F)
fwrite(occdat, "data/derivative_data/gridcell_data/env_forest/50_km/amph_comdat_occr.csv", row.names = F)

# takes 13 min

# # error here because I didn't use forest env
# tic()
# amph_comdat = get_gridcell_vert_parallel(trait_dat = amph_spdat, occ = occr, env = env, rich_min = 5, ncore = 7, nsim = 100)
# write.csv(amph_comdat, "data/derivative_data/gridcell_data/amphibians_comdat/amph_comdat_parallel_forestsOnly.csv", row.names = F)
# toc()


# Mammals -----------------------------------------------------------------


# load mammal occurrence data (cover0 includes all species, even if ranges cover less than 1/3 grid cell)
occr = list.files("data/derivative_data/pres_abs_cover0/resolution_50km/mammals/", pattern = "pres_abs",
                  full.names = T)
occr = lapply(occr, fread)
xy = occr[[1]][,1:2]
occr = foreach(i = 1:length(occr), .combine = "cbind") %do% {
  o = occr[[i]][,3:ncol(occr[[i]])]
}

# find duplicate column names
cnames = colnames(occr)
cnames.dups = duplicated(cnames)
cnames.dups = cnames[which(cnames.dups)]

# remove duplicate column names by combining the duplicated species and replacing in the datatable
for(i in cnames.dups) {
  sp.dup = occr[,which(colnames(occr) == i), with = F]
  sp.dup = apply(sp.dup, 1, sum)
  sp.dup = ifelse(sp.dup == 1, 1, 0)
  occr = occr[, !which(colnames(occr) == i), with = F]
  occr = occr[, i := sp.dup]
}

# alphabatize columns by species
occr = setcolorder(occr, sort(colnames(occr)))

occr = cbind(xy, occr)

colnames(occr)[1:2] = c("x", "y")
colnames(occr) = gsub(" ", "_", colnames(occr))

# list of species from IUCN range maps
spp <- colnames(occr)[-1:-2]
spp <- data.frame(sciname = gsub(" ","_",spp))

# tree.sp = read.csv("data/original/phylotrees/mammals/taxonomy_mamPhy_5911species.csv")
# tree.sp = tree.sp %>% 
#   dplyr::select(Species_Name) %>% 
#   rename(sciname = Species_Name)

# get species and traits in common between iucn and trait dataset
traits_mammals = inner_join(spp, traits, by = "sciname")
# traits_mammals = inner_join(traits_mammals, tree.sp, by = "sciname")

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
mammals_spdat = get_sp_dat(trait_dat = traits_mammals, occ = occr, env = env.forest, traitbase = "moura")
write.csv(mammals_spdat, "data/derivative_data/species_data/mammals_spdat_moura_forestsOnly.csv", row.names = F)

mammals_spdat = read.csv("data/derivative_data/species_data/mammals_spdat_moura_forestsOnly.csv")


tic()
mammals_comdat = gridcell_dat(occ = occr, tr = mammals_spdat, env = env.forest, reps = 100)
toc()

griddat = lapply(mammals_comdat, "[[", 1)
griddat = rbindlist(griddat)

# alphabetize occurrence dataframe by species
occdat = lapply(mammals_comdat, "[[", 2)
occdat = rbindlist(occdat, use.names = TRUE, fill = TRUE)
x.index = which(colnames(occdat) == "x")
y.index = which(colnames(occdat) == "y")
spnames = colnames(occdat)[-c(x.index, y.index)]
occdat = setcolorder(occdat, c('x', 'y', sort(spnames)))

# replace NAs with zeros
occdat[is.na(occdat)] <- 0


fwrite(griddat, 
       "data/derivative_data/gridcell_data/env_forest/50_km/mammals_comdat.csv", row.names = F)
fwrite(occdat, 
       "data/derivative_data/gridcell_data/env_forest/50_km/mammals_comdat_occr.csv", row.names = F)


# tic()
# mammals_comdat = get_gridcell_vert_parallel(trait_dat = mammals_spdat, occ = occr, env = env.forest, rich_min = 5, ncore = 7, nsim = 100)
# write.csv(mammals_comdat, "data/derivative_data/gridcell_data/env_forest/50_km/mammals_comdat.csv", row.names = F)
# toc()


# # takes 1210 sec
# tic()
# mammals_comdat = get_gridcell_vert_parallel(trait_dat = mammals_spdat, occ = occr, env = env, rich_min = 5, ncore = 7, nsim = 100)
# write.csv(mammals_comdat, "data/derivative_data/gridcell_data/mammals_comdat/mammals_comdat_parallel_forestsOnly.csv", row.names = F)
# toc()

# Reptiles --------------------------------------------------------------

# load reptile occurrence data (cover0 includes all species, even if ranges cover less than 1/3 grid cell)
occr = list.files("data/derivative_data/pres_abs_cover0/resolution_50km/reptiles/", pattern = "pres_abs",
                  full.names = T)
occr = lapply(occr, fread)
xy = occr[[1]][,1:2]
occr = foreach(i = 1:length(occr), .combine = "cbind") %do% {
  o = occr[[i]][,3:ncol(occr[[i]])]
}
occr = cbind(xy, occr)

colnames(occr)[1:2] = c("x", "y")
colnames(occr) = gsub(" ", "_", colnames(occr))

# list of species from IUCN range maps
spp <- colnames(occr)[-1:-2]
spp <- data.frame(sciname = gsub(" ","_",spp))

# tree = read.tree("data/original/phylotrees/squamates/1-s2.0-S0006320716301483-mmc6/IBC06772-mmc6.tre")
# tree.sp = data.frame(sciname = tree$tip.label)

# get species and traits in common between iucn and trait dataset
traits_rept = inner_join(spp, traits, by = "sciname")
# traits_squam = inner_join(traits_rept, tree.sp, by = "sciname")

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
rept_spdat = get_sp_dat(trait_dat = traits_rept, occ = occr, env = env.forest, traitbase = "moura")
write.csv(rept_spdat, "data/derivative_data/species_data/reptiles_spdat_moura_forestsOnly.csv", row.names = F)

rept_spdat = read.csv("data/derivative_data/species_data/reptiles_spdat_moura_forestsOnly.csv")

# time for 50 km grid cells: 12 min
tic()
rept_comdat = gridcell_dat(occ = occr, tr = rept_spdat, env = env.forest, reps = 100)
toc()

griddat = lapply(rept_comdat, "[[", 1)
griddat = rbindlist(griddat)

# alphabetize occurrence dataframe by species
occdat = lapply(rept_comdat, "[[", 2)
occdat = rbindlist(occdat, use.names = TRUE, fill = TRUE)
x.index = which(colnames(occdat) == "x")
y.index = which(colnames(occdat) == "y")
spnames = colnames(occdat)[-c(x.index, y.index)]
occdat = setcolorder(occdat, c('x', 'y', sort(spnames)))

# replace NAs with zeros
occdat[is.na(occdat)] <- 0


fwrite(griddat, 
          "data/derivative_data/gridcell_data/env_forest/50_km/reptiles_comdat.csv", row.names = F)
fwrite(occdat, 
       "data/derivative_data/gridcell_data/env_forest/50_km/reptiles_comdat_occr.csv", row.names = F)


# takes 2825 sec
# tic()
# rept_comdat = get_gridcell_vert_parallel(trait_dat = rept_spdat, occ = occr, env = env, rich_min = 5, ncore = 7, nsim = 100)
# write.csv(rept_comdat, "data/derivative_data/gridcell_data/reptiles_comdat/rept_comdat_parallel_forestsOnly.csv", row.names = F)
# toc()



# Birds -------------------------------------------------------------------

# load bird (breeding + resident) occurrence data (cover0 includes all species, even if ranges cover less than 1/3 grid cell)
occr = list.files("data/derivative_data/pres_abs_cover0/resolution_50km/birds_breed_resident/", pattern = "pres_abs",
                  full.names = T)
occr = lapply(occr, fread)
xy = occr[[1]][,1:2]
occr = foreach(i = 1:length(occr), .combine = "cbind") %do% {
  o = occr[[i]][,3:ncol(occr[[i]])]
}

# find duplicate column names
cnames = colnames(occr)
cnames.dups = duplicated(cnames)
cnames.dups = unique(cnames[which(cnames.dups)])

# remove duplicate column names by combining the duplicated species and replacing in the datatable
for(i in cnames.dups) {
  sp.dup = occr[,which(colnames(occr) == i), with = F]
  sp.dup = apply(sp.dup, 1, sum)
  sp.dup = ifelse(sp.dup == 1, 1, 0)
  occr = occr[, !which(colnames(occr) == i), with = F]
  occr = occr[, c(i) := sp.dup]
}

# alphabatize columns by species
occr = setcolorder(occr, sort(colnames(occr)))

occr = cbind(xy, occr)

colnames(occr)[1:2] = c("x", "y")
colnames(occr) = gsub(" ", "_", colnames(occr))

# list of species from IUCN range maps
spp <- colnames(occr)[-1:-2]
spp <- data.frame(sciname = gsub(" ","_",spp))

# tree.sp = read.csv("data/original/phylotrees/mammals/taxonomy_mamPhy_5911species.csv")
# tree.sp = tree.sp %>% 
#   dplyr::select(Species_Name) %>% 
#   rename(sciname = Species_Name)

# Elton trait data

traits.birds <- read.csv("data/original/trait_data/elton_traits/BirdFuncDat.txt", sep = "\t")
names(traits.birds)

# Get traits of interest - ForStrat-Value, BodyMass-Value
traits.birds = traits.birds %>% 
  dplyr::select(PassNonPass:BLFamilyEnglish, Scientific:PelagicSpecialist, Nocturnal:BodyMass.Value) %>% 
  rename(body_size = BodyMass.Value,
         sciname = Scientific,
         family = BLFamilyLatin) %>% 
  mutate(sciname = gsub(" ","_", sciname)) %>% 
  drop_na()

names(traits.birds)

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
traits.birds = traits.birds %>% 
  left_join(traits.birds.moura, by = "sciname") %>% 
  filter(!(EcoMar == 1 & EcoTer == 0 & EcoFresh == 0)) %>% 
  filter(!(EcoMar == 1 & EcoTer == 0 & EcoFresh == 1))

# drop aerial only species
traits.birds = traits.birds %>% 
  filter(!(ForStrat.aerial == 1 & ForStrat.canopy == 0 & ForStrat.midhigh == 0 & ForStrat.understory == 0 & ForStrat.ground == 0 & ForStrat.wataroundsurf == 0 & ForStrat.watbelowsurf == 0))

# subset occurrence data to species that occur in both iucn and trait datasets
occr = occr %>% 
  dplyr::select(x,y,traits_birds$sciname)

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
traits.birds = traits.birds %>% 
  filter(sciname %in% colnames(occr))

# test to make sure species in trait and occr dataframes are in the same order
sum(traits.birds$sciname != colnames(occr[3:ncol(occr)]))


# * - species level data ------------------------------------------------------
bird_spdat = get_sp_dat(trait_dat = traits.birds, occ = occr, env = env, filter_col = NA, 
                        traitbase = "elton")
write.csv(bird_spdat, "data/derivative_data/species_data/birds_spdat_elton_forestsOnly.csv", row.names = F)

bird_spdat = read.csv("data/derivative_data/species_data/birds_spdat_elton_forestsOnly.csv")

# add arboreal, terrestrial, and fossorial classifications to species data for proportion calculations
bird_spdat = bird_spdat %>% 
  mutate(Arb = case_when(ForStrat.understory > 0 | ForStrat.midhigh > 0 |
                           ForStrat.canopy > 0 | ForStrat.aerial > 0 ~ 1,
                         .default = 0),
         Ter = case_when(ForStrat.ground > 0 | ForStrat.wataroundsurf > 0 |
                           ForStrat.watbelowsurf > 0 ~ 1, .default = 0),
         Fos = 0) # no fossorial birds

dn = traits %>% dplyr::select(sciname, Diu, Noc)
bird_spdat = bird_spdat %>% left_join(dn, by = "sciname")

tic()
birds_comdat = gridcell_dat(occ = occr, tr = bird_spdat, env = env.forest, reps = 100)
toc()

griddat = lapply(birds_comdat, "[[", 1)
griddat = rbindlist(griddat)

# alphabetize occurrence dataframe by species
occdat = lapply(birds_comdat, "[[", 2)
occdat = rbindlist(occdat, use.names = TRUE, fill = TRUE)
x.index = which(colnames(occdat) == "x")
y.index = which(colnames(occdat) == "y")
spnames = colnames(occdat)[-c(x.index, y.index)]
occdat = setcolorder(occdat, c('x', 'y', sort(spnames)))

# replace NAs with zeros
occdat[is.na(occdat)] <- 0


fwrite(griddat, 
       "data/derivative_data/gridcell_data/env_forest/50_km/birds_comdat.csv", row.names = F)
fwrite(occdat, 
       "data/derivative_data/gridcell_data/env_forest/50_km/birds_comdat_occr.csv", row.names = F)


# # takes 7051 sec sec
# tic()
# bird_comdat = get_gridcell_vert_parallel(trait_dat = bird_spdat, occ = occr, env = env, rich_min = 5, ncore = 7, nsim = 100)
# write.csv(bird_comdat, "data/derivative_data/gridcell_data/birds_comdat/birds_comdat_parallel_elton_forestsOnly.csv", row.names = F)
# toc()




