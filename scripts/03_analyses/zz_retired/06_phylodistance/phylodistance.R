library(picante)
library(foreach)

data(phylocom)
names(phylocom)

phy = phylocom$phylo
commtest = phylocom$sample
traits = phylocom$traits

phy
comm
traits

prunedphy = prune.sample(comm, phy)
prunedphy

par(mfrow = c(2,3))
for (i in row.names(comm)) {
  plot(prunedphy, show.tip.label = F, main = i)
  tiplabels(tip = which(prunedphy$tip.label %in% names(which(comm[i,]>0))), pch = 19, cex = 2)
}

par(mfrow = c(2,2))
for (i in row.names(comm)) {
  plot(phy, show.tip.label = F, main = i)
  tiplabels(pch=22, col = traits[,i], bg = traits[,i], cex = 1.5)
}

pd.result = pd(comm, phy, include.root = T)
pd.result

phydist = cophenetic(phy)
ses.mpd.result = ses.mpd(comm, phydist, null.model = "taxa.labels", abundance.weighted = F, runs = 99)
ses.mpd.result


# Amphibians --------------------------------------------------------------
tree = read.tree("data/original/phylotrees/amph_shl_new_Consensus_7238.tre")
tree

sp = read.csv("data/derivative_data/species_data/amph_spdat_moura_forestsOnly.csv")
tree.sp = tree$tip.label
head(sort(tree.sp))

# community data matrix for amphibians
comm = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/amph_comdat_occr.csv")

# fix hyphens in column names
colnames(comm)[which(colnames(comm) == "Pristimantis_w.nigrum")] = "Pristimantis_w-nigrum"
colnames(comm)[which(colnames(comm) == "Scinax_x.signatus")] = "Scinax_x-signatus"

#rownames(comm) = paste0(1:nrow(comm))
commdat = as.matrix(comm)
rownames(commdat) = paste0(commdat[,"x"], "_", commdat[,"y"])
commdat = commdat[,-c(1,2)]
griddat = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/amph_comdat.csv")

# moura and vertlife tree use same taxonomic backbone, so no need for phylomatcher
match = match.phylo.comm(phy = tree, comm = commdat)

phydist = cophenetic(match$phy)

# SES mpd by geographic realm
# takes 30 min
tic()
griddat.sub = griddat %>% filter(rich >= 5)
commdat = commdat[which(griddat$rich >= 5),]

ses.mpd = foreach(i = unique(griddat.sub$realm)) %do% {
  griddat.realm = griddat.sub[which(griddat.sub$realm==i),]
  comm.realm = commdat[which(griddat.sub$realm == i),]
  
  # remove species not in the realm
  comm.realm = comm.realm[,which(apply(comm.realm,2,sum)>0)]
  
  # prune tree to realm
  match = match.phylo.comm(phy = tree, comm = comm.realm)
  
  # distance
  phydist = cophenetic(match$phy)
  
  ses.mpd.i = ses.mpd(samp = match$comm, dis = phydist, null.model = "richness", runs = 100)
  ses.mpd.i = cbind(griddat.realm, ses.mpd.i)

}
toc()


head(ses.mpd[[1]])
nrow(ses.mpd[[1]])

library(data.table)
library(tidyverse)
library(colorspace)
ses.mpd = rbindlist(ses.mpd)

#griddat = left_join(griddat, ses.mpd, by = c("x", "y"))

ggplot(ses.mpd, aes(x = vert.mean.ses, y = mpd.obs, color = y)) +
  geom_point(size = 0.5) +
  scale_color_continuous_divergingx("spectral") +
  theme_classic()

ggplot(ses.mpd, aes(x = vert.mean, y = mpd.obs, color = y)) +
  geom_point(size = 0.5) +
  scale_color_continuous_divergingx("spectral") +
  theme_classic()

ggplot(ses.mpd, aes(x = vert.mean.ses, y = mpd.obs.z, color = y)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  ggtitle("Amphibians") +
  scale_color_continuous_divergingx("spectral") +
  theme_classic()

ggplot(ses.mpd, aes(x = vert.mean, y = mpd.obs.z, color = y)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  ggtitle("Amphibians") +
  scale_color_continuous_divergingx("spectral") +
  theme_classic()

ggplot(ses.mpd, aes(x = vert.mean.ses, y = mpd.obs.z, color = biome, fill = biome)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  ggtitle("Amphibians") +
  #scale_color_continuous_divergingx("spectral") +
  theme_classic()

ggplot(ses.mpd, aes(x = vert.mean, y = mpd.obs.z, color = biome, fill = biome)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  ggtitle("Amphibians") +
  #scale_color_continuous_divergingx("spectral") +
  theme_classic()

ggplot(ses.mpd %>% filter(grepl("Tropical", biome)), aes(x = vert.mean.ses, y = mpd.obs.z, color = biome, fill = biome)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  ggtitle("Amphibians: Tropical") +
  #scale_color_continuous_divergingx("spectral") +
  theme_classic()

ggplot(ses.mpd %>% filter(grepl("Temperate", biome)), aes(x = vert.mean.ses, y = mpd.obs.z, color = biome, fill = biome)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  ggtitle("Amphibians: Temperate") +
  #scale_color_continuous_divergingx("spectral") +
  theme_classic()

ggplot(ses.mpd %>% filter(grepl("Tropical", biome)), aes(x = vert.mean, y = mpd.obs.z, color = biome, fill = biome)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  ggtitle("Amphibians: Tropical") +
  #scale_color_continuous_divergingx("spectral") +
  theme_classic()

ggplot(ses.mpd %>% filter(grepl("Temperate", biome)), aes(x = vert.mean, y = mpd.obs.z, color = biome, fill = biome)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  ggtitle("Amphibians: Temperate") +
  #scale_color_continuous_divergingx("spectral") +
  theme_classic()


# note: need to rerun comdat minus the two species that were removed if using phylodistance
fwrite(ses.mpd, file = "results/phylodistance/amphibians_phylodistance.csv", row.names = F)

ses.mpd = read.csv("results/phylodistance/amphibians_phylodistance.csv")

wd = vect("data/original/rnaturalearth_world.shp") %>% 
  crop(ext(-180,180,-60,83.64)) %>% 
  project("+proj=cea +datum=WGS84")

library(tidyterra)
ggplot(ses.mpd) +
  geom_spatvector(data = wd, fill = "gray90", color = "black") +
  geom_tile(aes(x, y, fill = mpd.obs.z)) +
  geom_spatvector(data = wd, fill = NA, color = "black") +
  scale_fill_continuous_divergingx("spectral") +
  coord_sf(crs = "+proj=cea +datum=WGS84") +
  ggtitle("Amphibians") +
  theme_classic()


# Squamates --------------------------------------------------------------
tree = read.tree("data/original/phylotrees/squamates/1-s2.0-S0006320716301483-mmc6/IBC06772-mmc6.tre")
tree

sp = read.csv("data/derivative_data/species_data/squamates_spdat_moura_forestsOnly.csv")
tree.sp = tree$tip.label
head(sort(tree.sp))

# community data matrix for amphibians
comm = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/squamates_comdat_occr.csv")

#rownames(comm) = paste0(1:nrow(comm))
commdat = as.matrix(comm)
rownames(commdat) = paste0(commdat[,"x"], "_", commdat[,"y"])
commdat = commdat[,-c(1,2)]
griddat = read.csv("data/derivative_data/gridcell_data/env_forest/50_km/squamates_comdat.csv")

# moura and vertlife tree use same taxonomic backbone, so no need for phylomatcher
match = match.phylo.comm(phy = tree, comm = commdat)

phydist = cophenetic(match$phy)

# SES mpd by geographic realm
# takes 30 min
tic()
griddat.sub = griddat %>% filter(rich >= 5)
commdat = commdat[which(griddat$rich >= 5),]

ses.mpd = foreach(i = unique(griddat.sub$realm)) %do% {
  griddat.realm = griddat.sub[which(griddat.sub$realm==i),]
  comm.realm = commdat[which(griddat.sub$realm == i),]
  
  # remove species not in the realm
  comm.realm = comm.realm[,which(apply(comm.realm,2,sum)>0)]
  
  # prune tree to realm
  match = match.phylo.comm(phy = tree, comm = comm.realm)
  
  # distance
  phydist = cophenetic(match$phy)
  
  ses.mpd.i = ses.mpd(samp = match$comm, dis = phydist, null.model = "richness", runs = 100)
  ses.mpd.i = cbind(griddat.realm, ses.mpd.i)
  
}
toc()


head(ses.mpd[[1]])
nrow(ses.mpd[[1]])

library(data.table)
library(tidyverse)
library(colorspace)
ses.mpd = rbindlist(ses.mpd)

#griddat = left_join(griddat, ses.mpd, by = c("x", "y"))

ggplot(ses.mpd, aes(x = vert.mean.ses, y = mpd.obs, color = y)) +
  geom_point(size = 0.5) +
  scale_color_continuous_divergingx("spectral") +
  theme_classic()

ggplot(ses.mpd, aes(x = vert.mean, y = mpd.obs, color = y)) +
  geom_point(size = 0.5) +
  scale_color_continuous_divergingx("spectral") +
  theme_classic()

ggplot(ses.mpd, aes(x = vert.mean.ses, y = mpd.obs.z, color = y)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  ggtitle("Squamates") +
  scale_color_continuous_divergingx("spectral") +
  theme_classic()

ggplot(ses.mpd, aes(x = vert.mean, y = mpd.obs.z, color = y)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  ggtitle("Squamates") +
  scale_color_continuous_divergingx("spectral") +
  theme_classic()

ggplot(ses.mpd, aes(x = vert.mean.ses, y = mpd.obs.z, color = biome, fill = biome)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  ggtitle("Squamates") +
  #scale_color_continuous_divergingx("spectral") +
  theme_classic()

ggplot(ses.mpd, aes(x = vert.mean, y = mpd.obs.z, color = biome, fill = biome)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  ggtitle("Squamates") +
  #scale_color_continuous_divergingx("spectral") +
  theme_classic()

ggplot(ses.mpd %>% filter(grepl("Tropical", biome)), aes(x = vert.mean.ses, y = mpd.obs.z, color = biome, fill = biome)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  ggtitle("Squamates: Tropical") +
  #scale_color_continuous_divergingx("spectral") +
  theme_classic()

ggplot(ses.mpd %>% filter(grepl("Temperate", biome)), aes(x = vert.mean.ses, y = mpd.obs.z, color = biome, fill = biome)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  ggtitle("Squamates: Temperate") +
  #scale_color_continuous_divergingx("spectral") +
  theme_classic()

ggplot(ses.mpd %>% filter(grepl("Tropical", biome)), aes(x = vert.mean, y = mpd.obs.z, color = biome, fill = biome)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  ggtitle("Squamates: Tropical") +
  #scale_color_continuous_divergingx("spectral") +
  theme_classic()

ggplot(ses.mpd %>% filter(grepl("Temperate", biome)), aes(x = vert.mean, y = mpd.obs.z, color = biome, fill = biome)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  ggtitle("Squamates: Temperate") +
  #scale_color_continuous_divergingx("spectral") +
  theme_classic()


# note: need to rerun comdat minus the two species that were removed if using phylodistance
fwrite(ses.mpd, file = "results/phylodistance/squamates_phylodistance.csv", row.names = F)

wd = vect("data/original/rnaturalearth_world.shp") %>% 
  crop(ext(-180,180,-60,83.64)) %>% 
  project("+proj=cea +datum=WGS84")
  
library(tidyterra)
ggplot(ses.mpd) +
  geom_spatvector(data = wd, fill = "gray90", color = "black") +
  geom_tile(aes(x, y, fill = mpd.obs.z)) +
  geom_spatvector(data = wd, fill = NA, color = "black") +
  scale_fill_continuous_divergingx("spectral") +
  coord_sf(crs = "+proj=cea +datum=WGS84") +
  ggtitle("Squamates") +
  theme_classic()
