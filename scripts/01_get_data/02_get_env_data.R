# This code loads all environmental variables, projects to cylindrical equal area and 50 km resolution
# Saves to a csv file (each row = 1 cell represented by x,y coords)

library(terra)
library(sf)
library(fasterize)
library(doParallel)
library(cluster)
library(foreach)

# change tempdir so that I can easily delete temp files that take up too much space
# tempdir(tmpdir = "tempfiles/")

# template raster used to get grid level pres/abs data at 111km
mapa = rast(extent = c(-20592508, 20588492, -5743602, 6573398), crs = "+proj=cea +datum=WGS84")
#res(mapa) = 111000
res(mapa) = 50000

## Current climate conditions (CHELSA v2.1)

## Chelsa climate data
## Downloaded from https://chelsa-climate.org/downloads/
bio = list.files(path = "data/original/env_data/chelsa/", pattern = ".tif", full.names = T)
#bio = lapply(bio, rast)
bio = bio[c(1:5,9:12)]
bio18 = bio[grepl("bio18", bio)] # precip warmest quarter

# Resmaple chelsa climate data to 50 km resolution and project to cylindrical equal area projection
cl = makeCluster(7)
registerDoParallel(cl)
bio = foreach(i = 1:length(bio), .packages = c("terra")) %dopar% {
  mapa = rast(extent = c(-20592508, 20588492, -5743602, 6573398), crs = "+proj=cea +datum=WGS84")
  res(mapa) = 50000
  r = rast(bio[[i]])
  r = project(r, "+proj=cea +datum=WGS84")
  r = resample(r, mapa, threads = T)
  writeRaster(r, paste0("data/derivative_data/resampled_env_rasters_50km/", names(r), ".tif"))
}
stopCluster(cl)

# bio18
mapa = rast(extent = c(-20592508, 20588492, -5743602, 6573398), crs = "+proj=cea +datum=WGS84")
res(mapa) = 50000
r = rast(bio18)
r = project(r, "+proj=cea +datum=WGS84")
r = resample(r, mapa, threads = T)
writeRaster(r, paste0("data/derivative_data/resampled_env_rasters_50km/", names(r), ".tif"))


# future climate ----------------------------------------------------------
# Load future climate data files and make raster stacks for each scenario
# gfdl-esm4-ssp585
bio.f.gfdl = list.files(path = "data/original/env_data/chelsa_future/2071-2100/gfdl-esm4-ssp585/", pattern = ".tif", full.names = T)
bio.f.gfdl = lapply(bio.f.gfdl, rast)
bio.f.gfdl = rast(bio.f.gfdl)

bio.f.ipsl = list.files(path = "data/original/env_data/chelsa_future/2071-2100/ipsl-cm61-lr_ssp585/", pattern = ".tif", full.names = T)
bio.f.ipsl = lapply(bio.f.ipsl, rast)
bio.f.ipsl = rast(bio.f.ipsl)

bio.f.mpi = list.files(path = "data/original/env_data/chelsa_future/2071-2100/mpi-esm1-2-hr_ssp585/", pattern = ".tif", full.names = T)
bio.f.mpi = lapply(bio.f.mpi, rast)
bio.f.mpi = rast(bio.f.mpi)

bio.f.mri = list.files(path = "data/original/env_data/chelsa_future/2071-2100/mri-esm2-0_ssp585/", pattern = ".tif", full.names = T)
bio.f.mri = lapply(bio.f.mri, rast)
bio.f.mri = rast(bio.f.mri)

bio.f.ukesm = list.files(path = "data/original/env_data/chelsa_future/2071-2100/ukesm1-0-II_ssp585/", pattern = ".tif", full.names = T)
bio.f.ukesm = lapply(bio.f.ukesm, rast)
bio.f.ukesm = rast(bio.f.ukesm)

bio.f = list(bio.f.gfdl, bio.f.ipsl, bio.f.mpi, bio.f.mri, bio.f.ukesm)
bio.f = lapply(bio.f, crop, ext(-180, 180, -64.7937, 90)) 

writeRaster(bio.f[[1]], "tempfiles/bio_future_gfdl.tif")
writeRaster(bio.f[[2]], "tempfiles/bio_future_ipsl.tif")
writeRaster(bio.f[[3]], "tempfiles/bio_future_mpi.tif")
writeRaster(bio.f[[4]], "tempfiles/bio_future_mri.tif")
writeRaster(bio.f[[5]], "tempfiles/bio_future_ukesm.tif")

# crop bio18
bio18.1 = bio.f.gfdl[[grepl("bio18", names(bio.f.gfdl))]]
bio18.2 = bio.f.ipsl[[grepl("bio18", names(bio.f.ipsl))]]
bio18.3 = bio.f.mpi[[grepl("bio18", names(bio.f.mpi))]]
bio18.4 = bio.f.mri[[grepl("bio18", names(bio.f.mri))]]
bio18.5 = bio.f.ukesm[[grepl("bio18", names(bio.f.ukesm))]]

bio18 = c(bio18.1, bio18.2, bio18.3, bio18.4, bio18.5)
bio18 = crop(bio18, ext(-180, 180, -64.7937, 90))


# Get raster stacks for each bio layer (5 scenarios per bio layer)
bio.f = list.files(path = "tempfiles/", pattern = ".tif", full.names = T)
bio.f = lapply(bio.f, rast)


mat.f = rast(lapply(bio.f, "[[", 1))
precip_ann.f = rast(lapply(bio.f, "[[", 2))
precip_wet.f = rast(lapply(bio.f, "[[", 3))
precip_dry.f = rast(lapply(bio.f, "[[", 4))
precip_sea.f = rast(lapply(bio.f, "[[", 5))
temp_diu.f = rast(lapply(bio.f, "[[", 6))
temp_sea.f = rast(lapply(bio.f, "[[", 7))
maxtemp_warm.f = rast(lapply(bio.f, "[[", 8))
min_temp_cold.f = rast(lapply(bio.f, "[[", 9))

# Take mean across the scenarios and resample to 50 km resolution and cylindrical equal area projection

mat.f.mean = mean(mat.f)
mat.f.proj = project(mat.f.mean, "+proj=cea +datum=WGS84")
mat.f.resamp = resample(mat.f.proj, mapa, threads = T)
names(mat.f.resamp) = "mat"
writeRaster(mat.f.resamp, "data/derivative_data/resampled_env_rasters_50km/chelsa_future/2071_2100/ensemble/chelsa_future_bio_01.tif")

maxtemp_warm.f = mean(maxtemp_warm.f)
maxtemp_warm.f = project(maxtemp_warm.f, "+proj=cea +datum=WGS84")
maxtemp_warm.f = resample(maxtemp_warm.f, mapa, threads = T)
names(maxtemp_warm.f) = "max_temp_warm"
writeRaster(maxtemp_warm.f, "data/derivative_data/resampled_env_rasters_50km/chelsa_future/2071_2100/ensemble/chelsa_bio_05.tif")

min_temp_cold.f = mean(min_temp_cold.f)
min_temp_cold.f = project(min_temp_cold.f, "+proj=cea +datum=WGS84")
min_temp_cold.f = resample(min_temp_cold.f, mapa, threads = T)
names(min_temp_cold.f) = "min_temp_cold"
writeRaster(min_temp_cold.f, "data/derivative_data/resampled_env_rasters_50km/chelsa_future/2071_2100/ensemble/chelsa_bio_06.tif")

# temp_sea.f = mean(temp_sea.f)
# temp_sea.f = project(temp_sea.f, "+proj=cea +datum=WGS84")
# temp_sea.f = resample(temp_sea.f, mapa, threads = T)
# names(temp_sea.f) = "temp_sea"
# writeRaster(temp_sea.f, "data/derivative_data/resampled_env_rasters_50km/chelsa_future/2071_2100/ensemble/chelsa_bio_04.tif")

temp_diu.f = mean(temp_diu.f)
temp_diu.f = project(temp_diu.f, "+proj=cea +datum=WGS84")
temp_diu.f = resample(temp_diu.f, mapa, threads = T)
names(temp_diu.f) = "temp_diu"
writeRaster(temp_diu.f, "data/derivative_data/resampled_env_rasters_50km/chelsa_future/2071_2100/ensemble/chelsa_bio_02.tif")

precip_wet.f = mean(precip_wet.f)
precip_wet.f = project(precip_wet.f, "+proj=cea +datum=WGS84")
precip_wet.f = resample(precip_wet.f, mapa, threads = T)
names(precip_wet.f) = "precip_wet"
writeRaster(precip_wet.f, "data/derivative_data/resampled_env_rasters_50km/chelsa_future/2071_2100/ensemble/chelsa_bio_13.tif")

precip_dry.f = mean(precip_dry.f)
precip_dry.f = project(precip_dry.f, "+proj=cea +datum=WGS84")
precip_dry.f = resample(precip_dry.f, mapa, threads = T)
names(precip_dry.f) = "precip_dry"
writeRaster(precip_dry.f, "data/derivative_data/resampled_env_rasters_50km/chelsa_future/2071_2100/ensemble/chelsa_bio_14.tif")

# precip_ann.f = mean(precip_ann.f)
# precip_ann.f = project(precip_ann.f, "+proj=cea +datum=WGS84")
# precip_ann.f = resample(precip_ann.f, mapa, threads = T)
# writeRaster(precip_ann.f, "data/derivative_data/resampled_env_rasters/chelsa_future/2071_2100/ensemble/chelsa_bio_12.tif")

precip_sea.f = mean(precip_sea.f)
precip_sea.f = project(precip_sea.f, "+proj=cea +datum=WGS84")
precip_sea.f = resample(precip_sea.f, mapa, threads = T)
names(precip_sea.f) = "precip_sea"
writeRaster(precip_sea.f, "data/derivative_data/resampled_env_rasters_50km/chelsa_future/2071_2100/ensemble/chelsa_bio_15.tif")

# add bio18: average monthly precip during warmest quarter
precip_warm.f = mean(bio18)
precip_warm.f = project(precip_warm.f, "+proj=cea +datum=WGS84")
precip_warm.f = resample(precip_warm.f, mapa, threads = T)
names(precip_warm.f) = "precip_warm"
writeRaster(precip_warm.f, "data/derivative_data/resampled_env_rasters_50km/chelsa_future/2071_2100/ensemble/chelsa_bio_18.tif")

# canopy height
# canopy height from Lang et al. 2023 https://doi.org/10.1038/s41559-023-02206-6
# This layer represents canopy height after masking out non-woody and human-altered habitats
# based on Hansen et al. 2022 land cover.
# then aggregating and averaging to coarser resolution
canopy_height = rast("data/derivative_data/resampled_env_rasters_50km/canopy_height_lang2023/canopy_height_agg_ForestOnly.tif")
canopy_height = crop(canopy_height, ext(-180, 180, -64.7937, 90))
canopy_height = project(canopy_height, "+proj=cea +datum=WGS84")
canopy_height = resample(canopy_height, mapa)
names(canopy_height) = "canopy_height"
writeRaster(canopy_height, "data/derivative_data/resampled_env_rasters_50km/canopy_height_lang2023/canopy_height_agg_ForestOnly_cea.tif", overwrite = T)

# vegetation density
# Crowther et al. 2015 https://doi.org/10.1038/nature14967
veg_den = rast("data/original/vegetation_density/Crowther_Nature_Files_Revision_01_WGS84_GeoTiff/Crowther_Nature_Biome_Revision_01_WGS84_GeoTiff.tif")
veg_den = crop(veg_den, ext(-180, 180, -64.7937, 90))
veg_den = project(veg_den, "+proj=cea +datum=WGS84")
veg_den = resample(veg_den, mapa)
names(veg_den) = "veg_den"
writeRaster(veg_den, "data/derivative_data/resampled_env_rasters_50km/veg_den.tif")

# elevation
# SRTM30+ Global 1-km Digital Elevation Model (DEM): Version 11: Land Surface
# downloaded from ERDDAP
elev = rast("data/original/srtm30plus/srtm30plus.tif")
elev = crop(elev, ext(-180, 180, -64.7937, 90))
elev = project(elev, "+proj=cea +datum=WGS84")
elev = resample(elev, mapa)
names(elev) = "elev"
writeRaster(elev, "data/derivative_data/resampled_env_rasters_50km/elevation.tif")

## Climate Velocity
## Downloaded from [Dryad](http://datadryad.org/resource/doi:10.5061/dryad.b13j1).
## Citation: Sandel, B. et al. (2011) The influence of Late Quaternary climate-change
## velocity on species endemism. Science.

clim_velocity = rast("data/original/env_data/doi_10_5061_dryad_b13j1__v20111101/Dryad Archive Files/Velocity.tif")
clim_velocity = project(clim_velocity, mapa)
names(clim_velocity) = "clim_velocity"
writeRaster(clim_velocity, "data/derivative_data/resampled_env_rasters_50km/clim_velocity.tif", overwrite = T)



# add biorealms -----------------------------------------------------------

# realms
# from Dinnerstein et al. 2017 https://doi.org/10.1093/biosci/bix014
biorealm = vect("data/original/env_data/ecoregions/Ecoregions2017.shp")
  
biorealm = project(biorealm, "+proj=cea datum=WGS84")
realm = rasterize(biorealm, mapa, field = "REALM")
realm = subst(realm, "N/A", NA)
names(realm) = "realm"
#writeRaster(realm, "data/derivative_data/resampled_env_rasters_50km/realm.tif")

# biomes: also from Dinnerstein
biome = rasterize(biorealm, mapa, field = "BIOME_NAME")
biome = subst(biome, "N/A", NA)
names(biome) = "biome"
#writeRaster(biome, "data/derivative_data/resampled_env_rasters_50km/biome.tif")

# ecoregion: also from Dinnerstein
ecoregion = rasterize(biorealm, mapa, field = "ECO_NAME")
ecoregion = subst(ecoregion, "N/A", NA)
names(ecoregion) = "ecoregion"
#writeRaster(ecoregion, "data/derivative_data/resampled_env_rasters_50km/ecoregion.tif")

# add canopy height calculated without human modified land or bare ground
# (masked by land cover and aggregated with median)
canopy_height2 = rast("data/derivative_data/resampled_env_rasters_50km/canopy_height_lang2023/canopy_height_agg_ForestOnly_cea.tif")

# add foliage height diversity (masked by landcover same as canopy height and aggregated with median value)
fhd = rast("data/derivative_data/resampled_env_rasters_50km/gedi_fhd/fhd_mask_aggregated.tif")
#fhd = rast("data/original/gedi_fhd/gediv002_fhd-pai-1m-a0_vf_20190417_20230316_12000m.tif")
fhd = project(fhd$mean, "+proj=cea +datum=WGS84")
fhd= resample(fhd, mapa)
names(fhd) = "fhd"

# stack rasters
r = rast(list.files(path = "data/derivative_data/resampled_env_rasters_50km/", pattern = ".tif", full.names = T))
r$biome = biome
r$realm = realm
r$ecoregion = ecoregion
bionames = names(r)[3:12]
bionames = sapply(bionames, strsplit, split = "_")
bionames = sapply(bionames, "[[", 2)
#bionames[10:13] = c(hurs_mean", "hurs_range", "vpd_mean", "vpd_range")
names(r)[3:12] =  bionames
#r = c(r, realm, biome)
r = c(r, canopy_height2)
names(r)[18] = "canopy_height2"
r = c(r, fhd)

env_df = as.data.frame(r, xy = T)

# check units
# summary(env_df)

# rescale env data to proper units
env_df = env_df %>% 
  dplyr::mutate(bio4 = bio4/100) %>% 
  dplyr::rename(mat = bio1, temp_diu = bio2, temp_sea = bio4, tmax_warm = bio5, tmin_cold = bio6,
         precip_ann = bio12, precip_wet = bio13, precip_dry = bio14, precip_sea = bio15, precip_warm = bio18) %>% 
  dplyr::mutate(veg_complexity = canopy_height * veg_den)

write.csv(env_df, "data/derivative_data/env_data_50km.csv", row.names = F)


# bio1 = annual mean temp
# bio5 = mean daily maximum air temperature of the warmest month
# bio6 = mean daily minimum air temperature of the coldest month
# bio4 = temp sesaonality
# bio2 = mean diurnal air temp range
# bio13 = precip amount of the wettest month
# bio14 = precip amount of the driest month
# bio12 = annual precip amount
# bio15 = precip seasonality
# bio18 = precip of the warmest quarter
# hurs_mean = Mean monthly near-surface relative humidity
# hurs_range = Annual range of monthly near-surface relative humidity
# vpd_mean = Mean monthly vapor pressure deficit
# vpd_range = annual range of monthly vapor pressure deficity


