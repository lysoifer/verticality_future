# Download MODIS landcover data and soil data for regional extents
# The extents cover 1x1 degree gedi tiles
# RUN ON HIPERGATOR

library(microclim3D) # on github lysoifer/microclim3D
library(data.table)
library(terra)

start = "2019-01-01"
end = "2022-12-31"
outpath = "/orange/scheffers/lydia.soifer/big_data/land_use/landcover_modis/regions2/"
usr = "lysoifer"
pwd = "Timp4mEDa."

# These are the processed gedi data, filtering out points with forest loss, etc.
gedi = list.files("data/derivative_data/gedi_filtered/", full.names = T)

# Get extents of gedi data
e = list()
for(i in 1:length(gedi)) {
  g = fread(gedi[i], select = c("lat_lm_a0", "lon_lm_a0"))
  e[[i]] = ext(min(g$lon_lm_a0), max(g$lon_lm_a0), min(g$lat_lm_a0), max(g$lat_lm_a0))
}
regions = gsub(".csv", "", basename(gedi))
regions = gsub("gedi_", "", regions)
names(e) = regions

# Download modis landcover using modisfast wrapper
microclim3D::get_modis_lc(e, start, end, usr, pwd, collection = "MCD12Q1.061",
                          variables = c("LC_Type1", "QC"), regions = names(e),
                          outpath = outpath)



# Download soil data ------------------------------------------------------

for(i in 1:length(e)) {
  out = paste0("/orange/scheffers/lydia.soifer/big_data/soil/soilgrid/regions/", names(e)[i], "/")
  if(!dir.exists(out)) dir.create(out)
  e2 = extend(e[[i]], 1) # using this for danum because I am getting NA soil values at a couple of points
  roi = vect(e2, crs = "epsg:4326")
  get_soil(roi, outdir = paste0(out, "soilgrids.tif"), overwrite = T)
}





