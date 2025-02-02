# This script downloads CHELSA climate data for present and future time periods, canopy height, and elevation
# For future climate data, the five ssp585 scenarios were downloaded (I later average these to obtain a future climate layer)

# download climate data
library(terra)
library(numform)
library(rerddap)
library(tidyverse)

# REFERENCE PROJECTION MAP
mapa <- rast(xmin = -20592508, xmax = 20588492, ymin = -5743602, ymax = 6573398,
             crs = "+proj=cea +datum=WGS84")
res(mapa) <- 111000
values(mapa) = 1

## CHELSA V2.1 bioclimate data
bio_layers = c(1,2,4,5,6,12:18)
chelsa_bio_paths = paste0("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio",bio_layers ,"_1981-2010_V.2.1.tif")
dest = paste0("data/original/env_data/chelsa/CHELSA_bio",bio_layers ,"_1981-2010_V.2.1.tif")
for(i in 11:length(chelsa_bio_paths)){
  curl::curl_download(chelsa_bio_paths[i], dest[i])
}

bio = list.files(path = "data/original/env_data/chelsa/", pattern = ".tif", full.names = T)
bio = lapply(bio, rast)
bio = rast(bio)
bio = project(bio, mapa)

chelsa = rast("data/original/env_data/chelsa/CHELSA_bio1_1981-2010_V.2.1.tif")

bio_layers2 = c("hurs_mean", "hurs_range", "vpd_mean", "vpd_range", "gdd0")
chelsa_bio_paths2 = paste0("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_",bio_layers2 ,"_1981-2010_V.2.1.tif")

dest2 = paste0("data/original/env_data/chelsa/CHELSA_",bio_layers2 ,"_1981-2010_V.2.1.tif")
for(i in 1:length(chelsa_bio_paths2)){
  curl::curl_download(chelsa_bio_paths2[i], dest2[i])
}


## CHELSA future climate projections CMIP6
# 2071-2100_gfdl-esm4_ssp585_V.2.1
f = read.delim("data/original/env_data/chelsa_future/2071-2100/gfdl-esm4-ssp585/envidatS3paths.txt", header = F)[,1]
f = str_squish(f)
fname = str_split(f, "/")
fname = sapply(fname, "[[", 13)

for(i in 1:length(f)) {
  dest = paste0("data/original/env_data/chelsa_future/2071-2100/gfdl-esm4-ssp585/", fname[i])
  curl::curl_download(f[i], dest)
}

# add bio18
f = "https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/2071-2100/GFDL-ESM4/ssp585/bio/CHELSA_bio18_2071-2100_gfdl-esm4_ssp585_V.2.1.tif"
f = str_squish(f)
fname = str_split(f, "/")
fname = sapply(fname, "[[", 11)
dest = paste0("data/original/env_data/chelsa_future/2071-2100/gfdl-esm4-ssp585/", fname)
curl::curl_download(f, dest)

# 2071-2100_ispl-cm61-lr_ssp585
f = read.delim("data/original/env_data/chelsa_future/2071-2100/ipsl-cm61-lr_ssp585/envidatS3paths.txt", header = F)[,1]
f = str_squish(f)
fname = str_split(f, "/")
fname = sapply(fname, "[[", 13)

for(i in 1:length(f)) {
  dest = paste0("data/original/env_data/chelsa_future/2071-2100/ipsl-cm61-lr_ssp585/", fname[i])
  curl::curl_download(f[i], dest)
}

# add bio18
f = "https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/2071-2100/IPSL-CM6A-LR/ssp585/bio/CHELSA_bio18_2071-2100_ipsl-cm6a-lr_ssp585_V.2.1.tif"
f = str_squish(f)
fname = str_split(f, "/")
fname = sapply(fname, "[[", 11)
dest = paste0("data/original/env_data/chelsa_future/2071-2100/ipsl-cm61-lr_ssp585/", fname)
curl::curl_download(f, dest)


# 2071-2100_mpi-esm-1-2-hr_ssp585
f = read.delim("data/original/env_data/chelsa_future/2071-2100/mpi-esm1-2-hr_ssp585/envidatS3paths.txt", header = F)[,1]
f = str_squish(f)
fname = str_split(f, "/")
fname = sapply(fname, "[[", 13)

for(i in 1:length(f)) {
  dest = paste0("data/original/env_data/chelsa_future/2071-2100/mpi-esm1-2-hr_ssp585/", fname[i])
  curl::curl_download(f[i], dest)
}

# add bio18
f = "https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/2071-2100/MPI-ESM1-2-HR/ssp585/bio/CHELSA_bio18_2071-2100_mpi-esm1-2-hr_ssp585_V.2.1.tif"
f = str_squish(f)
fname = str_split(f, "/")
fname = sapply(fname, "[[", 11)
dest = paste0("data/original/env_data/chelsa_future/2071-2100/mpi-esm1-2-hr_ssp585/", fname)
curl::curl_download(f, dest)

# mri-esm2-0_ssp585
f = read.delim("data/original/env_data/chelsa_future/2071-2100/mri-esm2-0_ssp585/envidatS3paths.txt", header = F)[,1]
f = str_squish(f)
fname = str_split(f, "/")
fname = sapply(fname, "[[", 13)

for(i in 1:length(f)) {
  dest = paste0("data/original/env_data/chelsa_future/2071-2100/mri-esm2-0_ssp585/", fname[i])
  curl::curl_download(f[i], dest)
}

# add bio18
f = "https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/2071-2100/MRI-ESM2-0/ssp585/bio/CHELSA_bio18_2071-2100_mri-esm2-0_ssp585_V.2.1.tif"
f = str_squish(f)
fname = str_split(f, "/")
fname = sapply(fname, "[[", 11)
dest = paste0("data/original/env_data/chelsa_future/2071-2100/mri-esm2-0_ssp585/", fname)
curl::curl_download(f, dest)

# ukesm1-0-II_ssp585
f = read.delim("data/original/env_data/chelsa_future/2071-2100/ukesm1-0-II_ssp585/envidatS3paths.txt", header = F)[,1]
f = str_squish(f)
fname = str_split(f, "/")
fname = sapply(fname, "[[", 13)

for(i in 1:length(f)) {
  dest = paste0("data/original/env_data/chelsa_future/2071-2100/ukesm1-0-II_ssp585/", fname[i])
  curl::curl_download(f[i], dest)
}

# add bio18
f = "https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/2071-2100/UKESM1-0-LL/ssp585/bio/CHELSA_bio18_2071-2100_ukesm1-0-ll_ssp585_V.2.1.tif"
f = str_squish(f)
fname = str_split(f, "/")
fname = sapply(fname, "[[", 11)
dest = paste0("data/original/env_data/chelsa_future/2071-2100/ukesm1-0-II_ssp585/", fname)
curl::curl_download(f, dest)




## canopy height
## Lang, N., Jetz, W., Schindler, K. et al. A high-resolution canopy height model of the Earth. 
## Nat Ecol Evol 7, 1778–1789 (2023). https://doi.org/10.1038/s41559-023-02206-6
## Downloaded from google earth engine where I aggregated data to 1km resolution

canopy_rr1 = rast("data/original/canopy_height/coarsened/canopy_height_reduceRes_1km_1.tif")
canopy_rr2 = rast("data/original/canopy_height/coarsened/canopy_height_reduceRes_1km_2.tif")
canopy_height = mosaic(canopy_rr1, canopy_rr2, filename = "data/original/canopy_height/coarsened/canopy_height_reduceRes_mosaic.tif")

# elevation
browse("srtm30plus_v11_land")
srtm30_info = info("srtm30plus_v11_land")
lat = c(-89.99583333333334, 89.99583333333334)
lon = c(-179.99583333333334, 0)
srtm30 = griddap(srtm30_info, latitude = lat, longitude = lon, fields = "elev",
                 store = disk(path = "data/original/elev/"), read = F)

srtm_1 = rast("data/original/elev/8f929b1d9881beb1cb9a7a6f2893e27f.nc")

lon = c(0,179.99583333333334)
srtm30_2 = griddap(srtm30_info, latitude = lat, longitude = lon, fields = "elev",
                   store = disk(path = "data/original/elev2/"), read = F)

srtm_2 = rast("data/original/elev2/f4ec53517f84b8ee8920455d92104e6c.nc")

srtm = mosaic(srtm_1, srtm_2, filename = "data/original/srtm30plus/srtm30plus.tif")
