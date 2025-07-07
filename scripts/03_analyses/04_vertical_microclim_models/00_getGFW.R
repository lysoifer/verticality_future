# get global forest watch data for the tiles I am interested 
# run on hipergator
# last updated: June 12, 2025

library(gfcanalysis)
library(sf)

pts = data.frame(lon = c(-89.54, -91.08,-75.57, 137.77, 10.50, 117.69),
                 lat = c(50.75,37.69,-5.24, 50.8, 48.33, 4.95),
                 region = c("canada_boreal", "mark_twain_nf", "pacaya_samiria",
                            "russia_boreal", "europe_temperate", "danum_valley"))
pts = vect(pts, geom = c("lon", "lat"), crs = "epsg:4326")
pts = buffer(pts, 100000)
pts = st_as_sf(pts)

for(i in 1:nrow(pts)) {
  tiles = calc_gfc_tiles(pts[i,])
  download_tiles(tiles, output_folder = "/orange/scheffers/lydia.soifer/big_data/vegetation/forestloss_hansen/",
                 images = c("lossyear"),
                 dataset = "GFC-2024-v1.12")
}
