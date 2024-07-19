# get gridded presence absence data at 50 km resolution from polygon range maps
# Lydia Soifer

library(terra)
library(sf)
library(letsR)
library(data.table)
library(tidyverse)
library(tidyterra)


# Mammals -----------------------------------------------------------------

# data = st_read("./../data/ranges/MAMMALS_TERRESTRIAL_ONLY/MAMMALS_TERRESTRIAL_ONLY.shp")
# 
# # test behaviour of letsR when there are multiple polygons for one species
# u_arctos = data %>% filter(sci_name == "Ursus arctos")
# u_arctos_pa = lets.presab(u_arctos, resol = 1, xmn = -190, xmx = 190, ymn = 18, ymx = 80)
# View(u_arctos_pa$Presence_and_Absence_Matrix) # combines all polygons into 1 for the species
# plot(u_arctos_pa)
# 
# unique(u_arctos$seasonal)
# 
# u_arctos_nonbreed = lets.presab(u_arctos, resol = 1, xmn = -190, xmx = 190, ymn = 18, ymx = 80, seasonal = 3)
# plot(u_arctos_nonbreed)
# 
# u_arctos_breed = lets.presab(u_arctos, resol = 1, xmn = -190, xmx = 190, ymn = 18, ymx = 80, seasonal = 1)
# plot(u_arctos_breed)
# 
# rm(list = ls())
# gc()

# verticality data
vert = read.csv("data/original/Moura_arboreality_vertebrates.csv")

# mammal data
data = vect("data/original/ranges/MAMMALS_TERRESTRIAL_ONLY/MAMMALS_TERRESTRIAL_ONLY.shp")

vert_mammal = vert %>% dplyr::filter(Class == "Mammalia") %>% 
  rename(sci_name = "Scientific.Name")

# reduce data to only species for which we have verticality data
data = data %>% inner_join(vert_mammal, by = "sci_name")

# check if all polygons are valid
valid = is.valid(data)
for (n in which(valid == FALSE)) {print(n)}

# reproject first because reprojecting after making valid undoes valid polygons
data_reproj = project(data, "+proj=cea +datum=WGS84")
data_reproj_valid = makeValid(data_reproj)

# recheck polygon validity
valid_check = is.valid(data_reproj_valid)
which(valid_check == FALSE)
rm(data, data_reproj)
gc()

library(numform)
for (i in 1:ceiling(nrow(data_reproj_valid) / 1000)) {
  st = 1 + 1000*(i-1)
  en = st + 999
  en = ifelse(en > nrow(data_reproj_valid), nrow(data_reproj_valid), en)
  
  d = data_reproj_valid[st:en]
  
  writeVector(d, paste0("tempfiles/mammals_", f_pad_zero(i,2), ".shp"))
}

# split range file for processing
library(doParallel)
library(foreach)
cl = makeCluster(5)
registerDoParallel(cl)
foreach(i = 1:ceiling(nrow(data_reproj_valid) / 1000), .packages = c("letsR", "numform")) %dopar% {
  
  f = list.files(path = "tempfiles/", pattern = ".shp", full.names = T)
  f = f[[i]]
  f = vect(f)
  
  # st = 1 + 2000*(i-1)
  # en = st + 1999
  # en = ifelse(en > nrow(data_reproj_valid), nrow(data_reproj_valid), en)
  # 
  # d = data_reproj_valid[st:en]
  
  # presence: extant or probably extant
  # origin: native only
  # seasonal: all resident, breeding, non-breeding
  # keep cells in matrix where richness = 0
  presab = lets.presab(f, resol = 50000, 
                            crs = "+proj=cea +datum=WGS84", crs.grid = "+proj=cea +datum=WGS84",
                            count = T, cover = 0,
                            xmn = -20592508, xmx = 20588492, ymn = -5743602, ymx = 6573398,
                            remove.cells = FALSE)
  
  lets.save(presab, file = paste0("data/derivative_data/pres_abs_cover0/resolution_50km/mammals/letsR_mammals_", f_pad_zero(i,2), ".Rdata"))
  
  # write occr
  write.csv(presab$Presence_and_Absence_Matrix, 
            paste0("data/derivative_data/pres_abs_cover0/resolution_50km/mammals/pres_abs_mammals_", f_pad_zero(i,2), ".csv"), row.names = F)
  # write species list
  write.csv(data.frame(species=presab$Species_name), 
            paste0("data/derivative_data/pres_abs_cover0/resolution_50km/mammals/sp_list_mammals_", f_pad_zero(i,2), ".csv"), row.names = F)
}
stopCluster(cl)

# remove temp files
f.rm = list.files(path = "tempfiles/", full.names = T)
file.remove(f.rm)


# presence: extant or probably extant
# origin: native only
# seasonal: all resident, breeding, non-breeding
# keep cells in matrix where richness = 0
# presab_mammals = lets.presab(data_reproj_valid, resol = 111000, 
#                              crs = "+proj=cea +datum=WGS84", crs.grid = "+proj=cea +datum=WGS84",
#                              count = T, cover = 0.1,
#                              xmn = -20592508, xmx = 20588492, ymn = -5743602, ymx = 6573398,
#                              presence = c(1,2), origin = 1, seasonal = 1:5,
#                              remove.cells = FALSE)
# 
# lets.save(presab_mammals, file = "data/pres_abs_cover0/letsR_mammals.Rdata")
# 
# # write occr
# write.csv(presab_mammals$Presence_and_Absence_Matrix, "data/pres_abs_cover0/pres_abs_mammals.csv",row.names = F)
# # write species list
# write.csv(data.frame(species=presab_mammals$Species_name), "data/pres_abs_cover0/sp_list_mammals.csv",row.names = F)
# 
# splist = data.frame(sci_name = presab_mammals$Species_name)
# vert_mammal = vert_mammal %>% inner_join(splist, by = "sci_name")
# 
# # add verticality data to ranges
# mammals_vertmap = lets.maplizer(presab_mammals,
#                                 y = vert_mammal$Verticality,
#                                 z = vert_mammal$sci_name,
#                                 func = mean, ras = T)
# 
# plot(mammals_vertmap$Raster)
# mammals_vertmap$Matrix[1:5, 1:3]
# 
# mat = mammals_vertmap$Matrix
# colnames(mat) = c("x", "y", "vert_mean_mammals")
# 
# writeRaster(mammals_vertmap$Raster, filename = "data/verticality_rasters_cover0/mammals_verticality.tif")
# write.csv(mat, "data/verticality_mat_cover0.csv", row.names = F)
# 
# rm(data_reproj_valid, mammals_vertmap, n, presab_mammals, splist, valid, valid_check, vert_mammal)
# gc()


# Amphibians --------------------------------------------------------------

data = vect("data/original/ranges/AMPHIBIANS/AMPHIBIANS.shp")
gc()

vert_amph = vert %>% dplyr::filter(Class == "Amphibia") %>% 
  rename(sci_name = "Scientific.Name")

# reduce data to only species for which we have verticality data
data = data %>% inner_join(vert_amph, by = "sci_name") %>% 
  filter(presence == 1 | presence == 2) %>% 
  filter(origin == 1)
gc()

# check if all polygons are valid
#valid = is.valid(data)
#for (n in which(valid == FALSE)) {print(n)}

# reproject first because reprojecting after making valid undoes valid polygons
data_reproj = project(data, "+proj=cea +datum=WGS84")
rm(data)
gc()
data_reproj_valid = makeValid(data_reproj)

# recheck polygon validity
valid_check = is.valid(data_reproj_valid)
which(valid_check == FALSE)
rm(data, data_reproj)
gc()

# presence: extant or probably extant
# origin: native only
# seasonal: all resident, breeding, non-breeding
# keep cells in matrix where richness = 0
# presab_amph = lets.presab(data_reproj_valid, resol = 111000, 
#                           crs = "+proj=cea +datum=WGS84", crs.grid = "+proj=cea +datum=WGS84",
#                           count = T, cover = 0.05,
#                           xmn = -20592508, xmx = 20588492, ymn = -5743602, ymx = 6573398,
#                           remove.cells = FALSE)
# 
# lets.save(presab_amph, file = "data/derivative_data/pres_abs_cover0.05/letsR_amph.Rdata")

# # write occr
# write.csv(presab_amph$Presence_and_Absence_Matrix, "data/derivative_data/pres_abs_cover0.05/pres_abs_amph.csv",row.names = F)
# # write species list
# write.csv(data.frame(species=presab_amph$Species_name), "data/derivative_data/pres_abs_cover0.05/sp_list_amph.csv",row.names = F)

presab_amph = lets.presab(data_reproj_valid, resol = 50000, 
                          crs = "+proj=cea +datum=WGS84", crs.grid = "+proj=cea +datum=WGS84",
                          count = T, cover = 0,
                          xmn = -20592508, xmx = 20588492, ymn = -5743602, ymx = 6573398,
                          remove.cells = FALSE)

lets.save(presab_amph, file = "data/derivative_data/pres_abs_cover0/resolution_50km/letsR_amph.Rdata")

# write occr
fwrite(presab_amph$Presence_and_Absence_Matrix, "data/derivative_data/pres_abs_cover0/resolution_50km/pres_abs_amph.csv",row.names = F)
# write species list
write.csv(data.frame(species=presab_amph$Species_name), "data/derivative_data/pres_abs_cover0/resolution_50km/sp_list_amph.csv",row.names = F)

# add amphibian verticality data to verticality data matrix
vert_mat = read.csv("data/verticality_mat_cover0.csv") %>% 
  as.data.frame()

splist = data.frame(sci_name = presab_amph$Species_name)

amph_df = as.data.frame(data_reproj_valid)
length(unique(amph_df$sci_name))
sort(unique(amph_df$sci_name))[1:100]
sort(splist$sci_name[1:100])

# check why there is such a drastic reduction in num species after running lets.presab
# amph_dat = st_read("./../data/ranges/AMPHIBIANS/AMPHIBIANS.shp")
# a_hylonomos = amph_dat %>% filter(sci_name == "Adelastes hylonomos")
# plot(wrld_simpl$geometry)
# a_hylonomos
# plot(a_hylonomos$geometry, color = "red", add = T)
# 
# # range is super tiny in Venezuela
# # how does lets.preab deal with super tiny range
# 
# a_hylonomos_cea = st_transform(a_hylonomos, "+proj=cea +datum=WGS84")
# 
# presab_a_hylonomous = lets.presab(a_hylonomos_cea, resol = 111000, 
#                           crs = "+proj=cea +datum=WGS84", crs.grid = "+proj=cea +datum=WGS84",
#                           count = T, cover = 0.3,
#                           xmn = -20592508, xmx = 20588492, ymn = -5743602, ymx = 6573398,
#                           remove.cells = T)
# 
# presab_a_hylonomous = lets.presab(a_hylonomos_cea, resol = 111000, 
#                                   crs = "+proj=cea +datum=WGS84", crs.grid = "+proj=cea +datum=WGS84",
#                                   count = T, cover = 0,
#                                   xmn = -20592508, xmx = 20588492, ymn = -5743602, ymx = 6573398,
#                                   remove.cells = T)
# 
# a_mar = amph_dat %>% filter(sci_name == "Adelophryne maranguapensis")
# a_mar

# so demanding the range covers 30% of the cell removes a bunch of species with very tiny ranges

vert_amph = vert_amph %>% inner_join(splist, by = "sci_name")

# add verticality data to ranges
amph_vertmap = lets.maplizer(presab_amph,
                             y = vert_amph$Verticality,
                             z = vert_amph$sci_name,
                             func = mean, ras = T)

plot(amph_vertmap$Raster)
amph_vertmap$Matrix[1:5, 1:3]

df_amph = amph_vertmap$Matrix %>% as.data.frame()
colnames(df_amph) = c("x", "y", "vert_mean_amphibian")

# add amphibian verticality to verticality data frame
vert_mat = left_join(vert_mat, df_amph, by = c("x", "y"))
head(vert_mat)

writeRaster(amph_vertmap$Raster, filename = "data/verticality_rasters_cover0/amphibian_verticality.tif",
            overwrite = T)
write.csv(vert_mat, "data/verticality_mat_cover0.csv", row.names = F)

rm(data_reproj_valid, mammals_vertmap, n, presab_mammals, splist, valid, valid_check, vert_amph)
rm(a_hylonomos, a_hylonomos_cea, a_mar, amph_dat, amph_df, amph_vertmap, df_amph, presab_a_hylonomous, presab_amph, wrld_simpl)
gc()



# Reptiles ----------------------------------------------------------------

data = vect("data/original/ranges/reptiles/Gard_1_7_ranges.shp")

vert_rept = vert %>% dplyr::filter(Class == "Reptilia") %>% 
  rename(sci_name = "Scientific.Name")

# reduce data to only species for which we have verticality data
data = data %>% inner_join(vert_rept, by = c("binomial" = "sci_name"))
gc()

# check if all polygons are valid
#valid = is.valid(data)
#for (n in which(valid == FALSE)) {print(n)}

# reproject first because reprojecting after making valid undoes valid polygons
data_reproj = project(data, "+proj=cea +datum=WGS84")
rm(data)
gc()
data_reproj_valid = makeValid(data_reproj)

# recheck polygon validity
#valid_check = is.valid(data_reproj_valid)
#which(valid_check == FALSE)
rm(data_reproj)
gc()

library(numform)
for (i in 1:ceiling(nrow(data_reproj_valid) / 1000)) {
  st = 1 + 1000*(i-1)
  en = st + 999
  en = ifelse(en > nrow(data_reproj_valid), nrow(data_reproj_valid), en)
  
  d = data_reproj_valid[st:en]
  
  writeVector(d, paste0("tempfiles/rept_", f_pad_zero(i,2), ".shp"))
}

# split range file for processing
library(doParallel)
library(foreach)
cl = makeCluster(5)
registerDoParallel(cl)
foreach(i = 1:ceiling(nrow(data_reproj_valid) / 1000), .packages = c("letsR", "numform")) %dopar% {
  
  f = list.files(path = "tempfiles/", pattern = ".shp", full.names = T)
  f = f[[i]]
  f = vect(f)
  
  # st = 1 + 2000*(i-1)
  # en = st + 1999
  # en = ifelse(en > nrow(data_reproj_valid), nrow(data_reproj_valid), en)
  # 
  # d = data_reproj_valid[st:en]
  
  # presence: extant or probably extant
  # origin: native only
  # seasonal: all resident, breeding, non-breeding
  # keep cells in matrix where richness = 0
  presab_rept = lets.presab(f, resol = 50000, 
                            crs = "+proj=cea +datum=WGS84", crs.grid = "+proj=cea +datum=WGS84",
                            count = T, cover = 0,
                            xmn = -20592508, xmx = 20588492, ymn = -5743602, ymx = 6573398,
                            remove.cells = FALSE)
  
  lets.save(presab_rept, file = paste0("data/derivative_data/pres_abs_cover0/resolution_50km/reptiles/letsR_rept_", f_pad_zero(i,2), ".Rdata"))
  
  # write occr
  write.csv(presab_rept$Presence_and_Absence_Matrix, 
            paste0("data/derivative_data/pres_abs_cover0/resolution_50km/reptiles/pres_abs_rept_", f_pad_zero(i,2), ".csv"), row.names = F)
  # write species list
  write.csv(data.frame(species=presab_rept$Species_name), 
            paste0("data/derivative_data/pres_abs_cover0/resolution_50km/reptiles/sp_list_rept_", f_pad_zero(i,2), ".csv"), row.names = F)
}
stopCluster(cl)

# remove temp files
f.rm = list.files(path = "tempfiles/", full.names = T)
file.remove(f.rm)

# presence: extant or probably extant
# origin: native only
# seasonal: all resident, breeding, non-breeding
# keep cells in matrix where richness = 0
# presab_rept = lets.presab(data_reproj_valid, resol = 50000, 
#                           crs = "+proj=cea +datum=WGS84", crs.grid = "+proj=cea +datum=WGS84",
#                           count = T, cover = 0,
#                           xmn = -20592508, xmx = 20588492, ymn = -5743602, ymx = 6573398,
#                           remove.cells = FALSE)

# #lets.save(presab_rept, file = "data/pres_abs_cover0/letsR_rept.Rdata")
# lets.save(presab_rept, file = "data/derivative_data/pres_abs_cover0/resolution_50km/letsR_rept.Rdata")
# 
# # write occr
# write.csv(presab_rept$Presence_and_Absence_Matrix, "data/derivative_data/pres_abs_cover0/resolution_50km/pres_abs_rept.csv",row.names = F)
# # write species list
# write.csv(data.frame(species=presab_rept$Species_name), "data/derivative_data/pres_abs_cover0/resolution_50km/sp_list_rept.csv",row.names = F)

# # add reptile verticality data to verticality data matrix
# vert_mat = read.csv("data/verticality_mat_cover0.csv") %>% 
#   as.data.frame()
# 
# splist = data.frame(sci_name = presab_rept$Species_name)
# 
# rept_df = as.data.frame(data_reproj_valid)
# length(unique(rept_df$binomial))
# sort(unique(rept_df$binomial))[1:100]
# sort(splist$sci_name[1:100])
# 
# # vert_sf = st_read("./../data/ranges/reptiles/Gard_1_7_ranges.shp")
# # a_anz = vert_sf %>% filter(binomial == "Abronia anzuetoi")
# # plot(a_anz$geometry)
# 
# vert_rept = vert_rept %>% inner_join(splist, by = "sci_name")
# 
# # add verticality data to ranges
# rept_vertmap = lets.maplizer(presab_rept,
#                              y = vert_rept$Verticality,
#                              z = vert_rept$sci_name,
#                              func = mean, ras = T)
# 
# plot(rept_vertmap$Raster)
# rept_vertmap$Matrix[1:5, 1:3]
# 
# df_rept = rept_vertmap$Matrix %>% as.data.frame()
# colnames(df_rept) = c("x", "y", "vert_mean_reptile")
# 
# # add amphibian verticality to verticality data frame
# vert_mat = left_join(vert_mat, df_rept, by = c("x", "y"))
# head(vert_mat)
# 
# writeRaster(rept_vertmap$Raster, filename = "data/verticality_rasters_cover0/reptile_verticality.tif")
# write.csv(vert_mat, "data/verticality_mat_cover0.csv", row.names = F)
# 
# gc()


# Birds Breeding and resident -------------------------------------------------------------------

vert_birds = vert %>% dplyr::filter(Class == "Aves") %>% 
  rename(sci_name = "Scientific.Name")

data_layers = st_layers("data/original/ranges/BOTW_2022.2/BOTW.gdb")

data_proxy = vect("data/original/ranges/BOTW_2022.2/BOTW.gdb", layer = "All_Species", proxy = T)

data = query(data_proxy, start = 1, vars = c("sisid", "sci_name", "presence", "origin", "seasonal", "Shape_Length", "Shape_Area"),
             where = "origin = 1 AND (presence = 1 OR presence = 2) AND (seasonal = 1 OR seasonal = 2)")

# reduce data to only species for which we have verticality data
data = data %>% inner_join(vert_birds, by = "sci_name")
gc()

# check if all polygons are valid
#valid = is.valid(data)
#for (n in which(valid == FALSE)) {print(n)}

# reproject first because reprojecting after making valid undoes valid polygons
data_reproj = project(data, "+proj=cea +datum=WGS84")
rm(data)
gc()
data_reproj_valid = makeValid(data_reproj)

# recheck polygon validity
#valid_check = is.valid(data_reproj_valid)
#which(valid_check == FALSE)
rm(data_reproj)
gc()

library(numform)
for (i in 1:ceiling(nrow(data_reproj_valid) / 1000)) {
  st = 1 + 1000*(i-1)
  en = st + 999
  en = ifelse(en > nrow(data_reproj_valid), nrow(data_reproj_valid), en)
  
  d = data_reproj_valid[st:en]
  
  writeVector(d, paste0("tempfiles/birds_", f_pad_zero(i,2), ".shp"))
}

# split range file for processing
library(doParallel)
library(foreach)
cl = makeCluster(2)
registerDoParallel(cl)
foreach(i = c(5,8), .packages = c("letsR", "numform")) %dopar% {
  
  f = list.files(path = "tempfiles/", pattern = ".shp", full.names = T)
  f = f[[i]]
  f = vect(f)
  
  # st = 1 + 2000*(i-1)
  # en = st + 1999
  # en = ifelse(en > nrow(data_reproj_valid), nrow(data_reproj_valid), en)
  # 
  # d = data_reproj_valid[st:en]
  
  # presence: extant or probably extant
  # origin: native only
  # seasonal: all resident, breeding, non-breeding
  # keep cells in matrix where richness = 0
  presab = lets.presab(f, resol = 50000, 
                       crs = "+proj=cea +datum=WGS84", crs.grid = "+proj=cea +datum=WGS84",
                       count = T, cover = 0,
                       xmn = -20592508, xmx = 20588492, ymn = -5743602, ymx = 6573398,
                       remove.cells = FALSE)
  
  lets.save(presab, file = paste0("data/derivative_data/pres_abs_cover0/resolution_50km/birds_breed_resident/letsR_birds_", f_pad_zero(i,2), ".Rdata"))
  
  # write occr
  write.csv(presab$Presence_and_Absence_Matrix, 
            paste0("data/derivative_data/pres_abs_cover0/resolution_50km/birds_breed_resident/pres_abs_birds_", f_pad_zero(i,2), ".csv"), row.names = F)
  # write species list
  write.csv(data.frame(species=presab$Species_name), 
            paste0("data/derivative_data/pres_abs_cover0/resolution_50km/birds_breed_resident/sp_list_birds_", f_pad_zero(i,2), ".csv"), row.names = F)
}
stopCluster(cl)

# remove temp files
f.rm = list.files(path = "tempfiles/", full.names = T)
file.remove(f.rm)


# # presence: extant or probably extant
# # origin: native only
# # seasonal: all resident, breeding, non-breeding
# # keep cells in matrix where richness = 0
# presab_birdsbreed = lets.presab(data_reproj_valid, resol = 111000, 
#                           crs = "+proj=cea +datum=WGS84", crs.grid = "+proj=cea +datum=WGS84",
#                           count = T, cover = 0,
#                           xmn = -20592508, xmx = 20588492, ymn = -5743602, ymx = 6573398,
#                           remove.cells = FALSE)
# 
# lets.save(presab_birdsbreed, file = "data/pres_abs_cover0/letsR_birds_breed.Rdata")
# 
# # write occr
# fwrite(presab_birdsbreed$Presence_and_Absence_Matrix, "data/pres_abs_cover0/pres_abs_birds_breeding_resident.csv",row.names = F)
# # write species list
# write.csv(data.frame(species=presab_birdsbreed$Species_name), "data/pres_abs_cover0/sp_list_breeding_resident.csv",row.names = F)
# 
# # add bird verticality data to verticality data matrix
# vert_mat = read.csv("data/verticality_mat_cover0.csv") %>% 
#   as.data.frame()
# 
# splist = data.frame(sci_name = presab_birdsbreed$Species_name)
# 
# birdsbreed_df = as.data.frame(data_reproj_valid)
# length(unique(birdsbreed_df$sci_name))
# sort(unique(birdsbreed_df$sci_name))[1:100]
# sort(splist$sci_name[1:100])
# 
# a_kat = query(data_proxy, where = "sci_name = 'Acanthiza katherina'")
# a_kat
# plot(a_kat)
# # range is long and narrow, so is not included
# 
# vert_birdsbreed = vert_birds %>% inner_join(splist, by = "sci_name")
# 
# # add verticality data to ranges
# birdsbreed_vertmap = lets.maplizer(presab_birdsbreed,
#                              y = vert_birdsbreed$Verticality,
#                              z = vert_birdsbreed$sci_name,
#                              func = mean, ras = T)
# 
# plot(birdsbreed_vertmap$Raster)
# birdsbreed_vertmap$Matrix[1:5, 1:3]
# 
# df_birdsbreed = birdsbreed_vertmap$Matrix %>% as.data.frame()
# colnames(df_birdsbreed) = c("x", "y", "vert_mean_birds_breeding")
# 
# # add bird verticality to verticality data frame
# vert_mat = left_join(vert_mat, df_birdsbreed, by = c("x", "y"))
# head(vert_mat)
# 
# writeRaster(birdsbreed_vertmap$Raster, filename = "data/verticality_rasters_cover0/birds_breeding_resident_verticality.tif")
# write.csv(vert_mat, "data/verticality_mat_cover0.csv", row.names = F)
# 
# gc()



# Birds non-Breeding and resident -------------------------------------------------------------------

vert_birds = vert %>% dplyr::filter(Class == "Aves") %>% 
  rename(sci_name = "Scientific.Name")

data_layers = st_layers("data/original/ranges/BOTW_2022.2/BOTW.gdb")

data_proxy = vect("data/original/ranges/BOTW_2022.2/BOTW.gdb", layer = "All_Species", proxy = T)

data = query(data_proxy, start = 1, vars = c("sisid", "sci_name", "presence", "origin", "seasonal", "Shape_Length", "Shape_Area"),
             where = "origin = 1 AND (presence = 1 OR presence = 2) AND (seasonal = 1 OR seasonal = 3)")

df = as.data.frame(data)
rm(df)
gc()

# reduce data to only species for which we have verticality data
data = data %>% inner_join(vert_birds, by = "sci_name")
gc()

# check if all polygons are valid
#valid = is.valid(data)
#for (n in which(valid == FALSE)) {print(n)}

# reproject first because reprojecting after making valid undoes valid polygons
data_reproj = project(data, "+proj=cea +datum=WGS84")
rm(data)
gc()
data_reproj_valid = makeValid(data_reproj)

# recheck polygon validity
#valid_check = is.valid(data_reproj_valid)
#which(valid_check == FALSE)
rm(data_reproj)
gc()

# presence: extant or probably extant
# origin: native only
# seasonal: all resident, breeding, non-breeding
# keep cells in matrix where richness = 0
presab_birds_nonbreed = lets.presab(data_reproj_valid, resol = 111000, 
                                crs = "+proj=cea +datum=WGS84", crs.grid = "+proj=cea +datum=WGS84",
                                count = T, cover = 0,
                                xmn = -20592508, xmx = 20588492, ymn = -5743602, ymx = 6573398,
                                remove.cells = FALSE)

lets.save(presab_birds_nonbreed, file = "data/pres_abs_cover0/letsR_birds_nonbreed.Rdata")

load("data/pres_abs/letsR_birds_nonbreed.Rdata")
pam$Richness_Raster = terra::unwrap(pam$Richness_Raster)

presab_birds_nonbreed = pam
rm(pam)

# write occr
fwrite(presab_birds_nonbreed$Presence_and_Absence_Matrix, "data/pres_abs_cover0/pres_abs_birds_nonbreeding_resident.csv",row.names = F)
# write species list
write.csv(data.frame(species=presab_birds_nonbreed$Species_name), "data/pres_abs_cover0/sp_list_nonbreeding_resident.csv",row.names = F)

# add bird verticality data to verticality data matrix
vert_mat = read.csv("data/verticality_mat_cover0.csv") %>% 
  as.data.frame()

splist = data.frame(sci_name = presab_birds_nonbreed$Species_name)

birds_nonbreed_df = as.data.frame(data_reproj_valid)
length(unique(birds_nonbreed_df$sci_name))
# sort(unique(birds_nonbreed_df$sci_name))[1:100]
# sort(splist$sci_name[1:100])

vert_birds_nonbreed = vert_birds %>% inner_join(splist, by = "sci_name")

# add verticality data to ranges
birds_nonbreed_vertmap = lets.maplizer(presab_birds_nonbreed,
                                   y = vert_birds_nonbreed$Verticality,
                                   z = vert_birds_nonbreed$sci_name,
                                   func = mean, ras = T)

plot(birds_nonbreed_vertmap$Raster)
birds_nonbreed_vertmap$Matrix[1:5, 1:3]

df_birds_nonbreed = birds_nonbreed_vertmap$Matrix %>% as.data.frame()
colnames(df_birds_nonbreed) = c("x", "y", "vert_mean_birds_nonbreeding")

# add bird verticality to verticality data frame
vert_mat = left_join(vert_mat, df_birds_nonbreed, by = c("x", "y"))
head(vert_mat)

writeRaster(birds_nonbreed_vertmap$Raster, filename = "data/verticality_rasters_cover0/birds_nonbreeding_resident_verticality.tif")
write.csv(vert_mat, "data/verticality_mat_cover0.csv", row.names = F)

gc()










