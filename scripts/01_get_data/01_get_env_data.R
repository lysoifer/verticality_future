library(terra)
library(sf)
library(fasterize)

# change tempdir so that I can easily delete temp files that take up too much space
tempdir(tmpdir = "tempfiles/")

# template raster used to get grid level pres/abs data at 111km
mapa = rast(extent = c(-20592508, 20588492, -5743602, 6573398), crs = "+proj=cea +datum=WGS84")
res(mapa) = 111000

## Current climate conditions (CHELSA v2.1)

## Chelsa climate data
## Downloaded from https://chelsa-climate.org/downloads/
bio = list.files(path = "data/original/env_data/chelsa/", pattern = ".tif", full.names = T)
bio = lapply(bio, rast)
bio = rast(bio)
bio = crop(bio, ext(-180, 180, -64.7937, 90))

tempsea = rast("data/original/env_data/chelsa/CHELSA_bio4_1981-2010_V.2.1.tif")

mat = project(bio$`CHELSA_bio1_1981-2010_V.2.1`, "+proj=cea +datum=WGS84")
mat = resample(mat, mapa, threads = T)
writeRaster(mat, "data/derivative_data/resampled_env_rasters/chelsa_bio_01.tif")

tic()
max_temp_warm = project(bio$`CHELSA_bio5_1981-2010_V.2.1`, "+proj=cea +datum=WGS84")
max_temp_warm = resample(max_temp_warm, mapa, threads = T)
writeRaster(max_temp_warm, "data/derivative_data/resampled_env_rasters/chelsa_bio_05.tif")
toc()

tic()
min_temp_cold = project(bio$`CHELSA_bio6_1981-2010_V.2.1`, "+proj=cea +datum=WGS84")
toc()
tic()
min_temp_cold = resample(min_temp_cold, mapa, threads = T)
toc()
writeRaster(min_temp_cold, "data/derivative_data/resampled_env_rasters/chelsa_bio_06.tif")


temp_sea = project(bio$`CHELSA_bio4_1981-2010_V.2.1`, "+proj=cea +datum=WGS84")
temp_sea = resample(temp_sea, mapa, threads = T)
writeRaster(temp_sea, "data/derivative_data/resampled_env_rasters/chelsa_bio_04.tif")

temp_diu = project(bio$`CHELSA_bio2_1981-2010_V.2.1`, "+proj=cea +datum=WGS84")
temp_diu = resample(temp_diu, mapa, threads = T)
writeRaster(temp_diu, "data/derivative_data/resampled_env_rasters/chelsa_bio_02.tif")

precip_wet = project(bio$`CHELSA_bio13_1981-2010_V.2.1`, "+proj=cea +datum=WGS84")
precip_wet = resample(precip_wet, mapa, threads = T)
writeRaster(precip_wet, "data/derivative_data/resampled_env_rasters/chelsa_bio_13.tif")

precip_dry = project(bio$`CHELSA_bio14_1981-2010_V.2.1`, "+proj=cea +datum=WGS84")
precip_dry = resample(precip_dry, mapa, threads = T)
writeRaster(precip_dry, "data/derivative_data/resampled_env_rasters/chelsa_bio_14.tif")

precip_ann = project(bio$`CHELSA_bio12_1981-2010_V.2.1`, "+proj=cea +datum=WGS84")
precip_ann = resample(precip_ann, mapa, threads = T)
writeRaster(precip_ann, "data/derivative_data/resampled_env_rasters/chelsa_bio_12.tif")

precip_sea = project(bio$`CHELSA_bio15_1981-2010_V.2.1`, "+proj=cea +datum=WGS84")
precip_sea = resample(precip_sea, mapa, threads = T)
writeRaster(precip_sea, "data/derivative_data/resampled_env_rasters/chelsa_bio_15.tif")

hurs_mean = project(bio$`CHELSA_hurs_mean_1981-2010_V.2.1`, "+proj=cea +datum=WGS84")
hurs_mean = resample(hurs_mean, mapa, threads = T)
writeRaster(hurs_mean, "data/derivative_data/resampled_env_rasters/chelsa_hurs_mean.tif")

hurs_range = project(bio$`CHELSA_hurs_range_1981-2010_V.2.1`, "+proj=cea +datum=WGS84")
hurs_range = resample(hurs_range, mapa, threads = T)
writeRaster(hurs_range, "data/derivative_data/resampled_env_rasters/chelsa_hurs_range.tif")

vpd_mean = project(bio$`CHELSA_vpd_mean_1981-2010_V.2.1`, "+proj=cea +datum=WGS84")
vpd_mean = resample(vpd_mean, mapa, threads = T)
writeRaster(vpd_mean, "data/derivative_data/resampled_env_rasters/chelsa_vpd_mean.tif")

vpd_range = project(bio$`CHELSA_vpd_range_1981-2010_V.2.1`, "+proj=cea +datum=WGS84")
vpd_range = resample(vpd_range, mapa, threads = T)
writeRaster(vpd_range, "data/derivative_data/resampled_env_rasters/chelsa_vpd_range.tif")
tmpFiles(remove=T, old = T)


# future climate ----------------------------------------------------------

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


mat.f.mean = mean(mat.f)
mat.f.proj = project(mat.f.mean, "+proj=cea +datum=WGS84")
mat.f.resamp = resample(mat.f.proj, mapa, threads = T)
names(mat.f.resamp) = "mat"
writeRaster(mat.f.resamp, "data/derivative_data/resampled_env_rasters/chelsa_future/2071_2100/ensemble/chelsa_future_bio_01.tif")

maxtemp_warm.f = mean(maxtemp_warm.f)
maxtemp_warm.f = project(maxtemp_warm.f, "+proj=cea +datum=WGS84")
maxtemp_warm.f = resample(maxtemp_warm.f, mapa, threads = T)
names(maxtemp_warm.f) = "max_temp_warm"
writeRaster(maxtemp_warm.f, "data/derivative_data/resampled_env_rasters/chelsa_future/2071_2100/ensemble/chelsa_bio_05.tif")

min_temp_cold.f = mean(min_temp_cold.f)
min_temp_cold.f = project(min_temp_cold.f, "+proj=cea +datum=WGS84")
min_temp_cold.f = resample(min_temp_cold.f, mapa, threads = T)
names(min_temp_cold.f) = "min_temp_cold"
writeRaster(min_temp_cold.f, "data/derivative_data/resampled_env_rasters/chelsa_future/2071_2100/ensemble/chelsa_bio_06.tif")

temp_sea.f = mean(temp_sea.f)
temp_sea.f = project(temp_sea.f, "+proj=cea +datum=WGS84")
temp_sea.f = resample(temp_sea.f, mapa, threads = T)
names(temp_sea.f) = "temp_sea"
writeRaster(temp_sea.f, "data/derivative_data/resampled_env_rasters/chelsa_future/2071_2100/ensemble/chelsa_bio_04.tif")

temp_diu.f = mean(temp_diu.f)
temp_diu.f = project(temp_diu.f, "+proj=cea +datum=WGS84")
temp_diu.f = resample(temp_diu.f, mapa, threads = T)
names(temp_diu.f) = "temp_diu"
writeRaster(temp_diu.f, "data/derivative_data/resampled_env_rasters/chelsa_future/2071_2100/ensemble/chelsa_bio_02.tif")

precip_wet.f = mean(precip_wet.f)
precip_wet.f = project(precip_wet.f, "+proj=cea +datum=WGS84")
precip_wet.f = resample(precip_wet.f, mapa, threads = T)
names(precip_wet.f) = "precip_wet"
writeRaster(precip_wet.f, "data/derivative_data/resampled_env_rasters/chelsa_future/2071_2100/ensemble/chelsa_bio_13.tif")

precip_dry.f = mean(precip_dry.f)
precip_dry.f = project(precip_dry.f, "+proj=cea +datum=WGS84")
precip_dry.f = resample(precip_dry.f, mapa, threads = T)
names(precip_dry.f) = "precip_dry"
writeRaster(precip_dry.f, "data/derivative_data/resampled_env_rasters/chelsa_future/2071_2100/ensemble/chelsa_bio_14.tif")

precip_ann.f = mean(precip_ann.f)
precip_ann.f = project(precip_ann.f, "+proj=cea +datum=WGS84")
precip_ann.f = resample(precip_ann.f, mapa, threads = T)
writeRaster(precip_ann.f, "data/derivative_data/resampled_env_rasters/chelsa_future/2071_2100/ensemble/chelsa_bio_12.tif")

precip_sea.f = mean(precip_sea.f)
precip_sea.f = project(precip_sea.f, "+proj=cea +datum=WGS84")
precip_sea.f = resample(precip_sea.f, mapa, threads = T)
names(precip_sea.f) = "precip_sea"
writeRaster(precip_sea.f, "data/derivative_data/resampled_env_rasters/chelsa_future/2071_2100/ensemble/chelsa_bio_15.tif")



# canopy height
# canopy height from Lang et al. 2023 https://doi.org/10.1038/s41559-023-02206-6
canopy_height = rast("data/original/canopy_height/coarsened/canopy_height_reduceRes_mosaic.tif")
crs(canopy_height) = "+init=epsg:4326" # set crs because it was not being read correctly leading to warnings when projecting to cea
canopy_height = crop(canopy_height, ext(-180, 180, -64.7937, 90))
canopy_height = project(canopy_height, "+proj=cea +datum=WGS84")
canopy_height = resample(canopy_height, mapa)
names(canopy_height) = "canopy_height"
writeRaster(canopy_height, "data/derivative_data/resampled_env_rasters/canopy_height.tif")

# vegetation density
# Crowther et al. 2015 https://doi.org/10.1038/nature14967
veg_den = rast("data/original/vegetation_density/Crowther_Nature_Files_Revision_01_WGS84_GeoTiff/Crowther_Nature_Biome_Revision_01_WGS84_GeoTiff.tif")
veg_den = crop(veg_den, ext(-180, 180, -64.7937, 90))
veg_den = project(veg_den, "+proj=cea +datum=WGS84")
veg_den = resample(veg_den, mapa)
names(veg_den) = "veg_den"
writeRaster(veg_den, "data/derivative_data/resampled_env_rasters/veg_den.tif")

# elevation
# SRTM30+ Global 1-km Digital Elevation Model (DEM): Version 11: Land Surface
# downloaded from ERDDAP
elev = rast("data/original/srtm30plus/srtm30plus.tif")
elev = crop(elev, ext(-180, 180, -64.7937, 90))
elev = project(elev, "+proj=cea +datum=WGS84")
elev = resample(elev, mapa)
names(elev) = "elev"
writeRaster(elev, "data/derivative_data/resampled_env_rasters/elevation.tif")

## Climate Velocity
## Downloaded from [Dryad](http://datadryad.org/resource/doi:10.5061/dryad.b13j1).
## Citation: Sandel, B. et al. (2011) The influence of Late Quaternary climate-change
## velocity on species endemism. Science.

clim_velocity = rast("data/original/env_data/doi_10_5061_dryad_b13j1__v20111101/Dryad Archive Files/Velocity.tif")
clim_velocity = project(clim_velocity, mapa)
names(clim_velocity) = "clim_velocity"
writeRaster(clim_velocity, "data/derivative_data/resampled_env_rasters/clim_velocity.tif", overwrite = T)

# realms
# from Dinnerstein et al. 2017 https://doi.org/10.1093/biosci/bix014
biorealm = vect("data/original/env_data/ecoregions/Ecoregions2017.shp")
biorealm = project(biorealm, "+proj=cea datum=WGS84")
realm = rasterize(biorealm, mapa, field = "REALM")
realm = subst(realm, "N/A", NA)
names(realm) = "realm"
#writeRaster(realm, "data/derivative_data/resampled_env_rasters/realm.tif")

# biomes: also from Dinnerstein
biome = rasterize(biorealm, mapa, field = "BIOME_NAME")
biome = subst(biome, "N/A", NA)
names(biome) = "biome"
#writeRaster(biome, "data/derivative_data/resampled_env_rasters/biome.tif")

# stack rasters
r = rast(list.files(path = "data/derivative_data/resampled_env_rasters/", pattern = ".tif", full.names = T))
bionames = names(r)[2:14]
bionames = sapply(bionames, strsplit, split = "_")
bionames = sapply(bionames, "[[", 2)
bionames[10:13] = c("hurs_mean", "hurs_range", "vpd_mean", "vpd_range")
names(r)[2:14] =  bionames
r = c(r, realm, biome)

env_df = as.data.frame(r, xy = T)

# check units
# summary(env_df)

# rescale env data to proper units
env_df = env_df %>% 
  dplyr::mutate(bio4 = bio4/100) %>% 
  dplyr::rename(mat = bio1, temp_diu = bio2, temp_sea = bio4, tmax_warm = bio5, tmin_cold = bio6,
         precip_ann = bio12, precip_wet = bio13, precip_dry = bio14, precip_sea = bio15) %>% 
  dplyr::mutate(veg_complexity = canopy_height * veg_den)

write.csv(env_df, "data/derivative_data/env_data.csv", row.names = F)

# bio1 = annual mean temp
# bio5 = mean daily maximum air temperature of the warmest month
# bio6 = mean daily minimum air temperature of the coldest month
# bio4 = temp sesaonality
# bio2 = mean diurnal air temp range
# bio13 = precip amount of the wettest month
# bio14 = precip amount of the driest month
# bio12 = annual precip amount
# bio15 = precip seasonality
# hurs_mean = Mean monthly near-surface relative humidity
# hurs_range = Annual range of monthly near-surface relative humidity
# vpd_mean = Mean monthly vapor pressure deficit
# vpd_range = annual range of monthly vapor pressure deficity


