# run micropoint models for 100 points in each forested area
library(terra)
library(mcera5)
library(ncdf4)
library(data.table)
library(microclim3D)
library(microclimf)
library(foreach)
library(numform)

# read in gedi data
# danum gedi points are all labeled as mangroves - resselect point

### HARD CODES ###
#gedi = fread("data/derivative_data/gedi_filtered/gedi_pacaya_samiria.csv")

era = paste0("/orange/scheffers/dklinges/soiltemp/SoilTemp-big-data/era5/mcera5_",2019:2022, "_time.nc")
yrs = 2019:2022
land_sea_mask = rast("/blue/scheffers/dklinges/soiltemp/SoilTemp-daily-offset/data/enviro_parameters/land_cover/landprop_30km.tif")

# For all the regions to run in a loop
# landcover for each region

gedifiles = list.files("data/derivative_data/gedi_filtered/", full.names = T)

lcpath = list.dirs("/orange/scheffers/lydia.soifer/big_data/land_use/landcover_modis/regions2/data/", full.names = T, recursive = T)
lcpath = grep("MCD", lcpath, value = T)

# soil for each region
spath = list.files("/orange/scheffers/lydia.soifer/big_data/soil/soilgrid/regions/", recursive = T, full.names = T)

# global elevation and elevation of era5
elev_path = "/orange/scheffers/lydia.soifer/big_data/topography/elevation_geodataR/elevation/wc2.1_30s/wc2.1_30s_elev.tif"
dtmc_path = "/orange/scheffers/lydia.soifer/big_data/topography/dtmc_era5/demera5.tif"

# file paths to save output for each region
regions = sapply(spath, function(x)stringr::str_split_i(x, "/", 10))
op = paste0("/blue/scheffers/lydia.soifer/verticality_future/data/derivative_data/micropoint_prep/prep_micropoint_", regions, ".rds")



#########

## Functions ##

# prep micropoint is from the microclim3D package that I am working on 
# Here is the version used in this analysis
# As the package is under active development, the function may change overtime
prep_micropoint = function(e, start, end, era_path = NA, clima = NA, landcover_path, soilpath, elev_path, dtmc_path) {
  # extend requested extent (so distance weighting works)
  e_extend = extend(e, 1)
  # extract climate data
  # Run if file path to .nc file is provided (file path from running get_era5)
  if(!is.na(era_path)) {
    climr = extract_clima(era_path, 
                          start_time = as.POSIXlt(start, tz = "UTC"), 
                          end_time = as.POSIXlt(end, tz = "UTC"),
                          e_extend[1], e_extend[2], e_extend[3], e_extend[4],
                          format = "micropoint")
    climr = lapply(climr, terra::unwrap)
  }
  
  # Run this if climate rasters are provided from running gutted and sped
  # up version of extract_clima
  # right now, I am not going to unwrap clima because I am exporting to desktop
  # if(is(clima, "list")) {
  #   climr = lapply(clima, terra::unwrap)
  # }
  if(is(clima, "list")) {
    climr = clima
  }
  
  
  # VEGETATION PARAMETERS
  # landcover
  lc = modisfast::mf_import_data(landcover_path, collection = "MCD12Q1.061")
  lc = terra::project(lc, "epsg:4326", method = "mode")
  #lc = terra::crop(lc, clima[[2]])
  # reclassify landcover to categories used in microclimf
  lc = microclim3D:::reclass_landcover(lc[[grep("LC", names(lc), value = T)]], type = "modis_igbp")
  names(lc) = as.character(time(lc))
  # reproject landcover to match projection of gedi data ("epsg:4326)
  #lc = .reproj_bigr(lc, "epsg:4326", e, method = "near")
  
  # need to deal with PAI
  # will have to estimate monthly change in lai based on modis or landsat and apply offset to gedi
  
  # EXTENT SHOULD BE PER 1 X 1 DEGREE TILE
  # get habitat-based vegetation data
  #load("scripts/package_devel/data/globclim.rda") # global climate variables (see documentation in microclimf)
  # extract time
  # if climr is a list of climate over multiple years
  if(is(climr[[1]], "list")) {
    tme = climr[[1]][[1]]
  } else { # climr only represents one year
    tme = climr$obs_time
  }
  
  #globclim = microclimf::globclim
  #use landcover in the most recent year that will be modeled
  # this should be sufficient because modis documentation says that the data should not be used to evaluate landcover change over time
  vegp = microclimf::vegpfromhab(
    habitats = lc[[nlyr(lc)]], # this selects 2022 landcover
    clump0 = F, # issue with indexing of pai in clumpestimate
    tme = as.POSIXlt(tme),
    lat = mean(e[3], e[4]),
    long = mean(e[1], e[2])
  )
  #vegp$pai = rast(vegp$pai, crs = crs(vegp$hgt), extent = ext(vegp$hgt))
  #vegp = lapply(vegp, unwrap)
  #vegp$clump = rast(vegp$clump, crs = crs(vegp$hgt), extent = ext(vegp$hgt))
  
  # get leaf reflectance - can work on this for improving the model??
  # albedo = rast(list.files(albedopath, pattern = ".tif$", full.names = T), lyrs = seq(1, 23, 2))
  # eproj_alb = project(e, "epsg:4326", crs(albedo))
  # albedo = crop(albedo, eproj_alb)
  # albedo = project(albedo, "epsg:4326")
  # albedo = resample(albedo, vegp$hgt)
  # albedo = mean(albedo, na.rm = T)
  
  # can the model take monthly leafr values?
  # or use same months for pai and albedo - see example in package - errors here
  #leafr = leafrfromalb(pai = pai$m_1, x = vegp$x, alb = albedo)
  
  #vegp$leafr = leafr
  
  # GROUND PARAMETERS
  # Soil data is from soilgrids
  # base file names must be as follows in soilnames. These are the defaults when downloading data
  # from soilgrids.
  #soil = list.files(soilpath, full.names = T)
  #soilnames = c("clay_0-5cm_mean_30s.tif", "sand_0-5cm_mean_30s.tif", "silt_0-5cm_mean_30s.tif")
  #soil = soil[grepl(paste(c(soilnames[1], soilnames[2], soilnames [3]), collapse = "|"), soil)]
  soil = rast(soilpath)
  soil = microclim3D:::.reproj_bigr(soil, "epsg:4326", e, method = "bilinear")
  soil = crop(soil, e)
  #soil = resample(soil, vegp$hgt)
  soil = soil/10
  
  # extend elevation to remove edge effects when calculating aspect and slope
  elev = rast(elev_path)
  #eproj = project(extend(e, c(1,1,1,1)), "epsg:4326", crs(elev))
  #elev = crop(elev, eproj)
  elev = crop(elev, e_extend)
  #elev = project(elev, "epsg:4326")
  #elev = project(elev, vegp$hgt, align_only = T)
  #elev = resample(elev, vegp$hgt)
  
  asp = terra::terrain(elev, "aspect", unit = "degrees")
  slp = terra::terrain(elev, "slope", unit = "degrees")
  
  # era5 elevation for altitudinal correction
  dtmc = rast(dtmc_path)
  crs(dtmc) = "epsg:4326"
  dtmc = crop(dtmc, e_extend)
  
  return(list(climr = climr, vegp = vegp, soil = wrap(soil), elev = wrap(elev), asp = wrap(asp), slp = wrap(slp), dtmc = wrap(dtmc)))
}


gedi = gedifiles[4]
landcover_path = lcpath[4]
soilpath = spath[4]
outpath = op[4]


# This function preps the data
# First it extracts the climate data from era5 for 2019-2022
# Next, select gedi points in months with min cold and max warm temps
# Finally, run prep_micropoint, which serves as input for the microclimate models
# save the output so I can run the micropoint model on the desktop, because linux is returning NA values in rcpp function call within micropoint
model_prep = function(gedi, era, yrs, land_sea_mask, landcover_path, soilpath, elev_path, dtmc_path, outpath) {
  # Extract era5 climate data for each year that gedi data is available (and that we have downloaded era5 data) - 2019-2022
  # Doing this outspide of prep micropoint because I need to calculate coldest and warmest months
  # and I want to do it for 2019-2022
  # I provide the land sea mask separately so that I can use the global era5 data that Dave has already downloaded
  # takes about 3 min per year
  
  gedi = fread(gedi)
  
  # I want to do this for each year from 2019 to 2022
  # find extent of gedi and extend by one cell in each direction so distance weighting works on edges
  e = ext(min(gedi$lon_lm_a0), max(gedi$lon_lm_a0), min(gedi$lat_lm_a0), max(gedi$lat_lm_a0))
  e_extend = extend(e, 2) # extending so distance weighting works on edges
  
  clima = list()
  for(i in 1:length(yrs)) {
    st = paste0(yrs[i], "-01-01")
    en = paste0(yrs[i], "-12-31")
    c = microclim3D:::extract_clima_2(era[i], long_min = e_extend[1], long_max = e_extend[2], lat_min = e_extend[3], lat_max = e_extend[4],
                                      start_time = as.POSIXlt(st, tz = "UTC"),
                                      end_time = as.POSIXlt(en, tz = "UTC"),
                                      land_sea_mask = land_sea_mask,
                                      dtr_cor = T,
                                      format = "micropoint")
    clima[[i]] = c
  }
  names(clima) = as.character(2019:2022)
  
  # calculate coldest and warmest month for each year
  mcold = list()
  mhot = list()
  for(i in 1:length(yrs)) {
    climai = lapply(clima[[i]], terra::unwrap)
    
    # summarize daily temperatures
    temp = climai$temp
    
    # min and max temp per day
    days = rep(1:(nlyr(temp)/24), each = 24)
    temp.minday = tapp(temp, days, fun = min)
    temp.maxday = tapp(temp, days, fun = max)
    
    # mean min and max per month
    if(length(days)/24 == 365){
      mon = c(rep(1,31),rep(2,28),rep(3,31), rep(4,30), rep(5,31),
              rep(6,30), rep(7,31), rep(8,31), rep(9,30), rep(10,31),
              rep(11,30), rep(12,31))
    }
    if(length(days)/24==366){
      mon = c(rep(1,31),rep(2,28),rep(3,31), rep(4,30), rep(5,31),
              rep(6,30), rep(7,31), rep(8,31), rep(9,30), rep(10,31),
              rep(11,30), rep(12,31))
    }
    
    temp.minday = tapp(temp.minday, mon, fun = mean)
    temp.maxday = tapp(temp.maxday, mon, fun = mean)
    
    # find coldest (lowest min temps on average) and hottest (hottest max temps on average) month
    mcold[[i]] = which.min(temp.minday)
    mhot[[i]] = which.max(temp.maxday)
    
  }
  names(mcold) = as.character(yrs)
  names(mhot) = as.character(yrs)
  
  mcold = rast(mcold) # month with coldest min temps on average
  mhot = rast(mhot) # month with hottest max temps on average
  
  # Prepare micropoint inputs -----------------------------------------------
  
  # First select the gedi points that we will use to model
  
  # subset gedi to points occurring in coldest and hottest months for each year
  gedi[,month := month(as.Date(date))]
  gedi[,year := year(as.Date(date))]
  gedi[,yearmonth := format(date, "%Y-%m")]
  
  yrs.mcold = paste0(yrs, "-", f_pad_zero(values(mcold),2))
  yrs.mhot = paste0(yrs, "-", f_pad_zero(values(mhot),2))
  
  gedisub = gedi[(yearmonth %in% yrs.mcold) | (yearmonth %in% yrs.mhot),]
  
  # remove any points in water / not in forest as designated by modis
  lc = modisfast::mf_import_data(landcover_path, collection = "MCD12Q1.061")
  lc = terra::project(lc, "epsg:4326", method = "mode")
  lc = microclim3D:::reclass_landcover(lc[[grep("LC", names(lc), value = T)]], type = "modis_igbp")
  names(lc) = as.character(time(lc))
  
  # landcover selection based on 2022 landcover (this should be fine as modis documentation states that 
  # there is too much error to examine changes in landcover over time) Also conservative using latest year
  gedisub[, lc := terra::extract(lc[["2022-01-01"]], .SD)[,2], .SDcols = c("lon_lm_a0", "lat_lm_a0")]
  gedisub = gedisub[lc <=5,] # restricts to points that fall in a forest
  
  # require height > 5m
  gedisub = gedisub[rh_100_a0 > 5,]
  
  # restrict to power beams
  gedisub = gedisub[beam == "BEAM0101" | beam == "BEAM0110" | beam == "BEAM1000" | beam == "BEAM1011",]
  
  # need to identify if the location of the point matches months of tmin or tmax
  # First identify the month of coldest min temp based on location
  gedisub[, c("tmin_2019", "tmin_2020", "tmin_2021", "tmin_2022") := terra::extract(mcold, .SD)[,2:5], .SDcols = c("lon_lm_a0", "lat_lm_a0")]
  gedisub[, tmin_2019 := paste0("2019-", f_pad_zero(tmin_2019,2))]
  gedisub[, tmin_2020 := paste0("2020-", f_pad_zero(tmin_2020,2))]
  gedisub[, tmin_2021 := paste0("2021-", f_pad_zero(tmin_2021,2))]
  gedisub[, tmin_2022 := paste0("2022-", f_pad_zero(tmin_2022,2))]
  
  # match year of observation to month of coldest mintemp
  gedisamp.mcold = gedisub[yearmonth == tmin_2019 | yearmonth == tmin_2020 | yearmonth == tmin_2021 | yearmonth == tmin_2022,]
  
  # randomly sample 100 gedi points if we have more than 100 points in the subset
  if(nrow(gedisamp.mcold)>100) {
    gedisamp.mcold = gedisamp.mcold[sample(1:nrow(gedisamp.mcold), 100, replace = F),] 
  }
  
  # same for max temps
  gedisub[, c("tmax_2019", "tmax_2020", "tmax_2021", "tmax_2022") := terra::extract(mhot, .SD)[,2:5], .SDcols = c("lon_lm_a0", "lat_lm_a0")]
  gedisub[, tmax_2019 := paste0("2019-", f_pad_zero(tmax_2019,2))]
  gedisub[, tmax_2020 := paste0("2020-", f_pad_zero(tmax_2020,2))]
  gedisub[, tmax_2021 := paste0("2021-", f_pad_zero(tmax_2021,2))]
  gedisub[, tmax_2022 := paste0("2022-", f_pad_zero(tmax_2022,2))]
  
  gedisamp.mhot = gedisub[yearmonth == tmax_2019 | yearmonth == tmax_2020 | yearmonth == tmax_2021 | yearmonth == tmax_2022,]
  
  # randomly sample 100 gedi points if we have more than 100 points in the subset
  if(nrow(gedisamp.mhot)>100) {
    gedisamp.mhot = gedisamp.mhot[sample(1:nrow(gedisamp.mhot), 100, replace = F),] 
  }
  
  # don't need start and end because I extracted clim data already
  microin = prep_micropoint(e = e, start = NA, end = NA, era_path = NA, clima = clima, landcover_path, soilpath, elev_path, dtmc_path)
  
  # make sure all rasters are wrapped
  # files to export to run on local desktop
  # microin, gedi
  out = list(gedisamp.mcold, gedisamp.mhot, microin)
  names(out) = c("gedisamp.mcold", "gedisamp.mhot", "microin")
  saveRDS(out, outpath)
  
}

# rerun europe temperate with height > 5m
# temperate forests for months with min cold temps don't have 100 points because there
# aren't that many points with leaf-on conditions, which is needed for accurate vertical profile measurements
for(i in 1:6) {
  print(i)
  model_prep(gedi = gedifiles[i],
             era = era,
             yrs = yrs,
             land_sea_mask = land_sea_mask,
             landcover_path = lcpath[i],
             soilpath = spath[i],
             elev_path = elev_path,
             dtmc_path = dtmc_path,
             outpath = op[i])
}







# run micropoint on desktop
# run_micropoint(tme = microin$climr$obs_time, 
#                gedi = gedi.m12[1:2], 
#                climr = microin$climr[2:10], 
#                vegp = microin$vegp, 
#                soil = microin$soil,
#                elev = microin$elev,
#                asp = microin$asp,
#                slp = microin$slp,
#                dtmc = microin$dtmc,
#                method = "temporal_month",
#                plotout = NA,
#                n,
#                maxiter=100,
#                fout, 
#                modis_path = NA,
#                reqhgts = NA)

