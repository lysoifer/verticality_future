# Get data for modelling vertical microclimates in boreal, temperate, and tropical forests
# last updated June 12th, 2025
# run on hipergator

library(terra)
library(microclim3D) # installed from lysoifer/microclim3D
library(data.table)
library(sf)
library(gfcanalysis)
library(dplyr)
library(tidyterra)
library(fasterize)
library(ggplot2)


# FUNCTIONS ---------------------------------------------------------------

# find name of Hansen tile
find_gfw_tile = function(tiles, output_folder, 
                         images = c("treecover2000", "lossyear", "gain", "datamask"),
                         dataset = "GFC-2022-v1.10") 
{
  stopifnot(all(images %in% c("treecover2000", "lossyear", 
                              "gain", "datamask", "first", "last")))
  if (!file_test("-d", output_folder)) {
    stop("output_folder does not exist")
  }
  message(paste(dim(tiles)[1], "tiles to download/check."))
  successes <- 0
  failures <- 0
  skips <- 0
  for (n in 1:dim(tiles)[1]) {
    gfc_tile <- tiles[n, ]
    min_x <- st_bbox(gfc_tile)[1]
    max_y <- st_bbox(gfc_tile)[4]
    if (min_x < 0) {
      min_x <- paste0(sprintf("%03i", abs(min_x)), "W")
    }
    else {
      min_x <- paste0(sprintf("%03i", min_x), "E")
    }
    if (max_y < 0) {
      max_y <- paste0(sprintf("%02i", abs(max_y)), "S")
    }
    else {
      max_y <- paste0(sprintf("%02i", max_y), "N")
    }
    file_root <- paste0("Hansen_", dataset, "_")
    file_suffix <- paste0("_", max_y, "_", min_x, ".tif")
    filenames <- paste0(file_root, images, file_suffix)
    filenames = paste0(output_folder, filenames)
    return(filenames)
  }}


# Copernicus closed forest raster -----------------------------------------

# run on hipergator
# Select forested land cover
# slow to run

# lc = rast("/orange/scheffers/lydia.soifer/big_data/land_use/copernicus/PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif")
# mask landcovers that are not closed forest
# notforest = c(0,121,123,122,124,125,126,20,30,90,100,60,40,50,70,80,200)
# lc.closedforest = subst(lc, from =notforest, to = rep(NA, 17), filename = "data/derivative_data/land_cover/copernicus_closedForest.tif")
#writeRaster(lc, "data/derivative_data/land_cover/copernicus_closedForest.tif")


# SELECT POINTS -----------------------------------------------------------
# select points in boreal, temperate, and tropical forests and download gedi tiles

pts = data.frame(lon = c(-89.54, -91.08,-75.57, 
                         137.77, 10.50, 117.69),
                 lat = c(50.75,37.69,-5.24, 
                         50.8, 48.33, 4.95),
                 region = c("canada_boreal", "mark_twain_nf", "pacaya_samiria",
                            "russia_boreal", "europe_temperate", "danum_valley"))
pts = vect(pts, geom = c("lon", "lat"), crs = "epsg:4326")
pts = buffer(pts, 100000)


# download gedi data - done on desktop then transfered files to hipergator
# get gedi tiles for polygons of interest

gedi = list()
for(p in 1:nrow(pts)) {
  g = get_gedi_tiles(e = pts[p], fpath = "D:/data/gedi_L2b/tiles/", unzip = T)
  gedi[[p]] = g
}
names(gedi) = pts$region


# FILTER GEDI POINTS ---------------------------
# Extract ecoregion, land cover, and forestloss

# Copernicus closed forest landcover (to remove areas that are not forest)
# not closed forest revalued to NA
landcover = rast("data/derivative_data/land_cover/copernicus_closedForest.tif")

# Biomes
temp = rast("/orange/scheffers/lydia.soifer/big_data/climate/chelsa/monthly/tas/CHELSA_tas_01_2000_V.2.1.tif") # template for rasterizing biomes
biomes = vect("data/original/ecoregions/Ecoregions2017.shp")
biomes_key = biomes %>% 
  tidyterra::select(BIOME_NUM, BIOME_NAME) %>% 
  as.data.frame() %>% 
  distinct() %>% 
  filter(BIOME_NAME != "N/A")

biomes = biomes %>% 
  tidyterra::select(BIOME_NUM)

biomes.r = rasterize(biomes, temp, field = "BIOME_NUM")

# identify GEDI tiles
gedi = list()
for(p in 1:nrow(pts)) {
  g = get_gedi_tiles(e = pts[p], fpath = "/orange/scheffers/lydia.soifer/big_data/vegetation/gedi_L2b/tiles/", unzip = F, download = F)
  gedi[[p]] = g
}
names(gedi) = pts$region


# remove points with forest loss between 2001 and 2024, anything but closed forest land cover, anything above 500 m elevation

for(i in 1:length(gedi)) {
  g = fread(gedi[[i]])
  # extract land cover
  g[,landcover := terra::extract(landcover, .SD)[,2], .SDcols = c("lon_lm_a0", "lat_lm_a0")]
  
  # extract forestloss year (0 = no forest loss)
  e = st_as_sf(data.frame(lat = range(g$lat_lm_a0), 
                          lon = range(g$lon_lm_a0)),
               coords = c("lon", "lat"), crs = 4326)
  tiles = calc_gfc_tiles(e)
  gfw_path = find_gfw_tile(tiles,
                           output_folder = "/orange/scheffers/lydia.soifer/big_data/vegetation/forestloss_hansen/",
                           images = c("lossyear"),
                           dataset = "GFC-2024-v1.12")
  forestloss = rast(gfw_path)
  g[,forestloss := terra::extract(forestloss, .SD)[,2], .SDcols = c("lon_lm_a0", "lat_lm_a0")]
  
  # extract biome
  g[,biome := terra::extract(biomes.r, .SD)[,2], .SDcols = c("lon_lm_a0", "lat_lm_a0")]
  g = left_join(g, biomes_key, by = c("biome" = "BIOME_NUM"))
  
  # remove areas that have lost forest
  g = g[forestloss==0,]
  
  # remove areas that do not have closed forest landcover
  g = g[!is.na(landcover),]
  
  # remove elevation above 500 m (we will model lowland forests only for consistency)
  g = g[elev_lm_a0<500, ]
  
  # name of region
  g[, region_name := names(gedi)[i]]
  
  fname = paste0("data/derivative_data/gedi_filtered/gedi_", names(gedi)[i], ".csv")
  fwrite(g, file = fname, row.names = F)
}



# crop to the coldest month of the year
# coldest_month = c(1, 1, 7)
# for(p in 1:length(gedi.crop)) {
#   gedi.crop[[p]] = gedi.crop[[p]] %>% 
#     mutate(month = month(as.Date(date))) %>% 
#     filter(month == coldest_month[p])
# }
# 
# test = gedi.crop[[3]]
# 