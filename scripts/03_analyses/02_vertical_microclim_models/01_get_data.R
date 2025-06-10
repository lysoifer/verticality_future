# Get data for modelling vertical microclimates in boreal, temperate, and tropical forests

library(terra)
source("C:/Users/lydia.soifer/Dropbox (UFL)/projects/3d_sdms/scripts/package_devel/R/downloads/get_gedi_tiles.R")
source("C:/Users/lydia.soifer/Dropbox (UFL)/projects/3d_sdms/scripts/package_devel/R/vegetation/gedi_import.R")

# Model microclimates within protected areas (from UNEP-WCMC)
pas2 = vect("D:/data/protected_areas/unep-wcmc_wdpa/WDPA_Jun2025_Public_shp/WDPA_Jun2025_Public_shp_2/WDPA_Jun2025_Public_shp-polygons.shp",
            proxy = T)

# path of least resistance - just select some points that look forested from google earth imagery

pts = data.frame(lon = c(-89.54, -91.08,-75.57),
                 lat = c(50.75,37.69,-5.24),
                 region = c("canada", "mark_twain_nf", "pacaya_samiria"))
pts = vect(pts, geom = c("lon", "lat"), crs = "epsg:4326")
pts = buffer(pts, 10000)


# get gedi tiles for polygons of interest
gedi = list()
for(p in 1:nrow(pts)) {
  g = get_gedi_tile(e = pts[p], fpath = "D:/data/gedi_L2b/tiles/", unzip = T)
  gedi[[p]] = g
}
names(gedi) = pts$region

# crop gedi tiles to polygons of interest
# will need to adjust of I have multiple tiles per area
gedi.crop = list()
tiles = list()
for(p in 1:nrow(pts)) {
  tiles[[p]] = fread(gedi[[p]])
  e = ext(pts[p])
  # gedi.crop[[p]] = g %>% 
  #   filter(lon_lm_a0 > e[1] &
  #            lon_lm_a0 < e[2] &
  #            lat_lm_a0 > e[3] &
  #            lat_lm_a0 < e[4])
}

# crop to the coldest month of the year
coldest_month = c(1, 1, 7)
for(p in 1:length(gedi.crop)) {
  gedi.crop[[p]] = gedi.crop[[p]] %>% 
    mutate(month = month(as.Date(date))) %>% 
    filter(month == coldest_month[p])
}

test = gedi.crop[[3]]













