library(terra)
hgt = rast("data/derivative_data/resampled_env_rasters_50km/canopy_height.tif")
fhd = rast("data/original/gedi_fhd/gediv002_fhd-pai-1m-a0_vf_20190417_20230316_12000m.tif")

fhd = project(fhd, hgt)
fhd = resample(fhd, hgt)
fhd = fhd$mean
names(fhd) = "fhd"

r = c(hgt, fhd)
df = as.data.frame(r) %>% 
  drop_na()

ggplot(df, aes(canopy_height, fhd)) +
  geom_point(pch = ".", color = "gray80") +
  geom_smooth(method = "lm") +
  theme_classic()
