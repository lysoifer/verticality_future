# RUN MICROPOINT MODELS
# MUST RUN ON DESKTOP FOR NOW
# micropoint prep run on hipergator because large data is stored there
# I exported the inputed needed to run micropoint, including
# the env variables (climate, veg, soil params), and sampled gedi
# points for min temp of coldest month and max temp of warmest month.
# Months with min and max temps were determined using era5 climate grids
# For example, I calculated min temp of each day, then averaged min temps
# by month, then found the month with the coldest min temp. This was done
# for each year.
# In general months are more consistent in temperate/boreal than in tropics

library(microclim3D)
library(foreach)
library(data.table)
library(soiltexture)
library(dplyr)
library(tidyr)
soilparams = micropoint::soilparams

# Function to run micropoint - part of microclim3D package, but package currently isn't functioning
run_micropoint = function(tme, gedi, climr, vegp, soil, elev, asp, slp, dtmc,
                          method, plotout, n, maxiter, fout, modis_path, reqhgts = NA, vertpai_method = "pai") {
  # if(method == "vertprof") {
  #   vertmat = matrix(ncol = 7, nrow = 0)
  #   colnames(vertmat) = c("z", "tair", "canopy_height", "elev", "lon", "lat", "doy")
  #   write.table(vertmat, fout)
  # }
  
  # Unpack spatrasters
  # if we have multiple years of climate data
  if(is(climr[[1]], "list")) {
    for(a in 1:length(climr)) {
      climr[[a]] = lapply(climr[[a]], unwrap)
    }
  } else { # if there is only one year of climate data
    climr = unwrap(climr)
  }
  
  # unpack vegp
  vegp = lapply(vegp, unwrap)
  soil = unwrap(soil)
  elev = unwrap(elev)
  asp = unwrap(asp)
  slp = unwrap(slp)
  dtmc = unwrap(dtmc)
  
  # Loop through GEDI points and calculate vertical microclimate gradient for each
  v = foreach(i = 1:nrow(gedi), .combine = rbind) %do% {
    print(i)
    
    # If we are modelling for only the month in which the GEDI point was observed
    # Then we need to pull out the year and month of the observation
    # and select the climate data for the correct year
    if(method == "temporal_month") {
      #gedi[, year:= year(as.Date(date))]
      pyear = gedi[i, year]
      pmonth = gedi[i, month]
      climri = climr[[as.character(pyear)]]
      tmei = climri$obs_time
      climri = climri[2:10]
    }
    
    
    # extract climate data from raster at requested point
    # and apply distance weighting correction
    
    clim_point = microclim3D:::clim_distweight(lon = as.numeric(gedi[i, "lon_lm_a0"]), lat = as.numeric(gedi[i,"lat_lm_a0"]),
                                               climri, tmei) %>%
      dplyr::select(!timezone) %>%
      mutate(obs_time = format(obs_time, format = "%Y-%m-%d %H:%M:%S")) %>%
      as.data.frame()
    
    # extract dtmc and elevation at point
    dtmc_p = terra::extract(dtmc, gedi[i,c("lon_lm_a0", "lat_lm_a0")], ID = F)[1,1]
    elev_p = as.numeric(gedi[i, "elev_lm_a0"])
    
    # apply altitudinal correction
    altcor = microclim3D:::altcorrect_point(clim_point, dtmc_p, elev_p, altcor = 2)
    clim_point$temp = altcor[["tc"]]
    clim_point$pres = altcor[["pk"]]
    
    # vegetation parameters
    h = as.numeric(gedi[i, "rh_100_a0"]) # height
    pai = as.numeric(gedi[i, "pai_a0"]) #pai
    
    # get vertical pai profile
    pai_z = as.numeric(gedi[i, pai_a0:pai_l21])
    pavd = as.numeric(gedi[i,pavd_0_5:pavd_95_100]) # pavd
    
    paii = microclim3D:::.pai_vertprofile(pai_z, h, vertpai_method)
    
    # If we want to model climate over an entire year, we need monthly estimates of pai
    # MODIS provides temporally resolved PAI, but GEDI is a single time point
    # The way I have approached it here is to model temporal variation in GEDI PAI
    # assuming a constant offset with MODIS
    # I additionally assume that the temporal variation is constant across the vertical profile
    # which is probably not a good assumption as understory greenup is typically offset from
    # canopy greenup
    # I also should be estimating temporal variability in clumpiness, which I am not currently doing
    
    # get vertical profile by season
    if(method == "temporal_year") {
      modis = rast(modis_path)
      pai_z_season = pai_seasonality(gedi[i,], xy_crs = "epsg:4326", modis)
      pai = pai_z_season[[1]]
      paii = pai_z_season[[2]]
    }
    
    # test if height is greater than 5m (errors when there is only one paiz value - not sure why)
    #print(pai_z)
    print(h)
    
    # extract values from vegp raster to get veg parameters at point
    # pai and height are derived from GEDI
    # remaining is derived from auto-generated params based on habitat type and lat/lon
    vegp_p = microclim3D:::vegp_point(vegp, lon = as.numeric(gedi[i, "lon_lm_a0"]),
                                      lat = as.numeric(gedi[i, "lat_lm_a0"]),
                                      pai, h)
    
    # soil parameters
    # get soil parameters
    asp_p = terra::extract(asp$aspect, gedi[i, c("lon_lm_a0", "lat_lm_a0")])[1,2]
    slp_p = terra::extract(slp$slope, gedi[i, c("lon_lm_a0", "lat_lm_a0")])[1,2]
    
    # get ground parameters based on soilgrids and average values from micropoint::soilparams
    ground_p = microclim3D:::.ground_point(soil, asp_p, slp_p,
                                           gedi[i, lon_lm_a0],
                                           gedi[i, lat_lm_a0])
    
    # apply dtr_correct
    clim_point_correct = micropoint:::dtr_correct(clim_point, zref = 2,
                                                  lat = gedi[i, lat_lm_a0],
                                                  long = gedi[i, lon_lm_a0])
    
    # model hourly microclimate over one year
    # takes into account temporal variation in PAI by modeling by month and changing PAI with each month based on changes in MODIS PAI
    if(method == "temporal_year") {
      # NEEDS SOME WORK, SEE JUNE 9 2025 IN LAB NOTEBOOK
      mout = foreach(m=1:12, .combine = "rbind") %do% {
        vegp_p$pai = pai[m]
        hrs = which(month(clim_point$obs_time)==m)
        out = micropoint::runpointmodel(
          clim_point_correct[hrs,],
          reqhgt = reqhgt,
          vegp = vegp_p,
          paii = paii[,m],
          groundp = ground_p,
          lat = as.numeric(gedi[i, "lat_lm_a0"]),
          long = as.numeric(gedi[i, "lon_lm_a0"])
        )
      }
      
    }
    
    # model microclimates at specified reqhgts for the month of the gedi shot
    if(method == "temporal_month") {
      #m = month(as.Date(gedi[i,date]))
      hrs = which(month(clim_point$obs_time)==pmonth)
      req = reqhgts
      if(is.na(sum(reqhgts))) {req = c(0.15, 2, seq(5,h - h%%5, 5), h-0.1)}
      
      # test of max canopy height is in req, and if yes, don't model at top of canopy
      if(h %in% req) {req = req[req!=h]}
      # pai is currently constant
      mout = foreach(hi = req, .combine = rbind) %do% {
        out = micropoint::runpointmodel(
          clim_point_correct[hrs,],
          reqhgt = hi,
          vegp = vegp_p,
          paii = paii,
          groundp = ground_p,
          lat = as.numeric(gedi[i, "lat_lm_a0"]),
          long = as.numeric(gedi[i, "lon_lm_a0"]),
          maxiter = 100
        )
        out$canopy_hgt = vegp_p$h
        out$reqhgt = hi
        out
      }
    }
    
    pai_test = pai_z[pai_z>0]
    if(method == "vertprof" & h >= 6 & length(pai_test)>1 & min(pai > 0.2)) {
      # NEEDS A BIT OF WORK SEE JUNE9, 2025 IN LAB NOTEBOOK
      hour_yr = gedi[i, doy]*24-12 # noon on day of year that gedi shot was taken
      tair = c()
      height = c()
      hour = c()
      for(i in seq(hour_yr,hour_yr + 1,1)) {
        vertprof <- micropoint::plotprofile(clim_point_correct, hr = hour_yr, plotout = plotout,
                                            vegp_p, paii = paii, ground_p,
                                            lat = as.numeric(gedi[i, "lat_lm_a0"]),
                                            long= as.numeric(gedi[i, "lon_lm_a0"]),
                                            maxiter = maxiter, n = length(paii))
        tair = c(tair, vertprof$var)
        height = c(height, vertprof$z)
        hour = c(hour, rep(i, length(vertprof$z)))
      }
      
      mat = cbind(z = vertprof$z, var = vertprof$var,
                  canopy_height = vegp_p$h,
                  elev = elev_p,
                  lon = gedi[i,lon_lm_a0],
                  lat = gedi[i, lat_lm_a0],
                  doy = gedi[i, doy])
      colnames(mat)[2] = plotout
      
      # mat
      
      #vertmat = rbind(vertmat, mat)
      #write.csv(vertmat, file = fout, row.names = F)
    }
    #toc()
    mout$shot_num = gedi[i, shot_num]
    mout
  }
  #return(vertmat)
  if(!is.na(fout)) {
    if(!dir.exists(dirname(fout))) {dir.create(dirname(fout))}
    write.csv(v, file = fout, row.names = F)
  }
  return(v)
}



# Run vertical microclimate profiles --------------------------------------

# Run mintemp vertical profiles
indat = list.files("D:/projects/verticality_future/data/derivative_data/micropoint_prep/", full.names = T)
# run vertical profiles in the coldest month
for(i in 1:6) {
  dat = readRDS(indat[i])
  gedi.tmin = dat$gedisamp.mcold
  microin = dat$microin
  
  nm = basename(indat[i])
  nm = gsub("prep_micropoint_", "", nm)
  nm = gsub(".rds", "", nm)
  
  print(i)
  print(nm)
  print(nrow(gedi.tmin))
  print(unique(gedi.tmin$BIOME_NAME))
  
  
  run = run_micropoint(
    tme = microin$climr$`2019`$obs_time,
    gedi = gedi.tmin,
    climr = microin$climr,
    vegp = microin$vegp,
    soil = microin$soil,
    elev = microin$elev,
    asp = microin$asp,
    slp = microin$slp,
    dtmc = microin$dtmc,
    method = "temporal_month",
    reqhgts = NA,
    vertpai_method = "pai",
    fout = paste0("results/microclim_models/",nm,"/mintemp_cold.csv"))
  
}

# run vertical profiles in the hottest month - need to rerun danum with the larger soil grid (currently prepping micropoint inupt on hipergator)
for(i in 1:6) {
  dat = readRDS(indat[i])
  gedi.tmax = dat$gedisamp.mhot
  microin = dat$microin
  
  nm = basename(indat[i])
  nm = gsub("prep_micropoint_", "", nm)
  nm = gsub(".rds", "", nm)
  
  print(i)
  print(nm)
  print(nrow(gedi.tmax))
  print(unique(gedi.tmax$BIOME_NAME))
  
  
  run = run_micropoint(
    tme = microin$climr$`2019`$obs_time,
    gedi = gedi.tmax,
    climr = microin$climr,
    vegp = microin$vegp,
    soil = microin$soil,
    elev = microin$elev,
    asp = microin$asp,
    slp = microin$slp,
    dtmc = microin$dtmc,
    method = "temporal_month",
    reqhgts = NA,
    vertpai_method = "pai",
    fout = paste0("results/microclim_models/",nm,"/tmax_warm.csv"))
  
}



# VISUALIZATION -----------------------------------------------------------

# Load in modeling data
micromods = list.files("results/microclim_models/", recursive = T, full.names = T)
micromods.tmin = grep("mintemp", micromods, value = T) # points for the coldest month of the year
micromods.tmax = grep("tmax", micromods, value = T) # points for the warmest month of the year


# Make dataframes -----------------------------------------------------------

# dataframe for month with coldest min temp
mm.tmin = foreach(f = 1:length(micromods.tmin), .combine = rbind) %do% {
  fi = fread(micromods.tmin[f])
  nm = stringr::str_split_i(dirname(micromods.tmin[f]), "/", 3)
  fi$region = nm
  fi
}

# dataframe for month with hottest max temp
mm.tmax = foreach(f = 1:length(micromods.tmax), .combine = rbind) %do% {
  fi = fread(micromods.tmax[f])
  nm = stringr::str_split_i(dirname(micromods.tmax[f]), "/", 3)
  fi$region = nm
  fi
}

# dataframe for PAD
pad = foreach(f = 1:length(indat), .combine = rbind) %do% {
  dat = readRDS(indat[f])
  gedi.tmin = dat$gedisamp.mcold
  gedi.tmax = dat$gedisamp.mhot
  gedi = rbind(gedi.tmin, gedi.tmax, fill = T)
  cnames = grep("^pavd", names(gedi), value = T)
  cnames = c("shot_num", "rh_100_a0", grep("[0-9]$", cnames, value = T))
  gedi = gedi[, ..cnames]
  gedi[gedi==0] <- NA # convert zeros to NAs
  gedi = gedi %>% 
    pivot_longer(cols = 3:ncol(gedi), names_to = "height", values_to = "pavd") %>% 
    drop_na()
  
  # average by relhgt categories
  
  nm = basename(indat[f])
  nm = gsub("prep_micropoint_", "", nm)
  nm = gsub(".rds", "", nm)
  gedi$region = nm
  gedi
}

# summarize by relhgt
pad2 = pad %>% 
  mutate(h = case_when(height == "pavd_0_5" ~ 0.15,
                       height == "pavd_5_10" ~ 5,
                       height == "pavd_10_15" ~ 10,
                       height == "pavd_15_20" ~ 15,
                       height == "pavd_20_25" ~ 20,
                       height == "pavd_25_30" ~ 25,
                       height == "pavd_30_35" ~ 30,
                       height == "pavd_35_40" ~ 35,
                       height == "pavd_40_45" ~ 40),
         relhgt = h/rh_100_a0, # calculate relative height
         frelhgt = cut(relhgt, breaks = seq(0,1,0.1), include_lowest = T, right = F)) %>%  # break points
  group_by(frelhgt, region, shot_num) %>% 
  summarize(pavd = mean(pavd, na.rm = T)) %>% # mean pavd per relhgt 
  separate(frelhgt, into = c("minbin", "maxbin"), remove = F, sep = ",") %>% 
  mutate(minbin = as.numeric(gsub("\\[", "", minbin)),
         maxbin = as.numeric(gsub(")", "", maxbin)),
         region2 = case_when(region == "canada_boreal" ~ "Boreal Canada",
                             region == "russia_boreal" ~ "Boreal Russia",
                             region == "europe_temperate" ~ "Temperate Europe",
                             region == "mark_twain_nf" ~ "Temperate USA",
                             region == "pacaya_samiria" ~ "Tropical S America",
                             region == "danum_valley" ~ "Tropical Borneo"))
        
boreal = ggplot(pad2 %>% filter(region2 == "Boreal Canada")) +
  geom_line(aes(x = minbin, y = pavd, group = shot_num)) +
  stat_summary(aes(x = minbin, y = pavd), color = "red") +
  scale_y_continuous(limits = c(0,0.7)) +
  coord_flip() +
  theme_classic()

temperate = ggplot(pad2 %>% filter(region2 == "Temperate USA")) +
  geom_line(aes(x = minbin, y = pavd, group = shot_num)) +
  stat_summary(aes(x = minbin, y = pavd), color = "red") +
  scale_y_continuous(limits = c(0,0.7)) +
  coord_flip() +
  theme_classic()

tropical = ggplot(pad2 %>% filter(region2 == "Tropical S America")) +
  geom_line(aes(x = minbin, y = pavd, group = shot_num)) +
  stat_summary(aes(x = minbin, y = pavd), color = "red") +
  scale_y_continuous(limits = c(0,0.7)) +
  coord_flip() +
  theme_classic()

pavdplt = boreal / temperate / tropical  

# ground to canopy difs ---------------------------------------------------

# Calculate average difference between ground climate and climates at each height in the canopy

# for month with coldest min temp
df2.tmin = mm.tmin %>% 
  mutate(hr = hour(obs_time), yday = yday(obs_time)) %>% 
  group_by(yday, reqhgt, canopy_hgt, shot_num, region) %>% 
  summarise(tair.min = min(tair),
            tair.max = max(tair),
            tleaf.min = min(tleaf),
            tleaf.max = max(tleaf),
            rh.min = min(relhum),
            windspeed.mean = mean(windspeed)) %>% 
  group_by(reqhgt, canopy_hgt, shot_num, region) %>% 
  summarise(tair.min = mean(tair.min),
            tair.max = mean(tair.max),
            tleaf.min = mean(tleaf.min),
            tleaf.max = mean(tleaf.max),
            rh.min = mean(rh.min),
            windspeed.mean = mean(windspeed.mean)) %>% 
  group_by(shot_num, region) %>% 
  mutate(tair.min.ground = tair.min[reqhgt==0.15],
         tair.max.ground = tair.max[reqhgt == 0.15],
         tleaf.min.ground = tleaf.min[reqhgt == 0.15],
         tleaf.max.ground = tleaf.max[reqhgt == 0.15],
         rh.min.ground = rh.min[reqhgt == 0.15],
         windspeed.mean.ground = windspeed.mean[reqhgt==0.15],
         tair.min.dif = tair.min - tair.min.ground,
         tair.max.dif = tair.max - tair.max.ground,
         tleaf.min.dif = tleaf.min - tleaf.min.ground,
         tleaf.max.dif = tleaf.max - tleaf.max.ground,
         rh.min.dif = rh.min - rh.min.ground,
         windspeed.mean.dif = windspeed.mean - windspeed.mean.ground,
         relhgt = reqhgt/canopy_hgt,
         frelhgt = cut(relhgt, breaks = seq(0,1,0.1), include_lowest = T, right = F)) %>% 
  group_by(frelhgt, region) %>% 
  summarise(tair.min.dif = median(tair.min.dif),
            tair.min = median(tair.min),
            tair.max.dif = median(tair.max.dif),
            tair.max = median(tair.max),
            tleaf.min.dif = median(tleaf.min.dif),
            tleaf.min = median(tleaf.min),
            tleaf.max.dif = median(tleaf.max.dif),
            tleaf.max = median(tleaf.max),
            rh.min.dif = median(rh.min.dif),
            rh.min = median(rh.min),
            windspeed.mean.dif = median(windspeed.mean.dif),
            windspeed.mean = median(windspeed.mean)) %>% 
  separate(frelhgt, into = c("minbin", "maxbin"), remove = F, sep = ",") %>% 
  mutate(minbin = as.numeric(gsub("\\[", "", minbin)),
         maxbin = as.numeric(gsub(")", "", maxbin)),
         lat = case_when(grepl("boreal", region) ~ "boreal",
                         region == "europe_temperate" ~ "temperate",
                         region == "mark_twain_nf" ~ "temperate",
                         region == "danum_valley" ~ "tropical",
                         region == "pacaya_samiria" ~ "tropical"),
         region2 = case_when(region == "canada_boreal" ~ "Boreal Canada",
                             region == "russia_boreal" ~ "Boreal Russia",
                             region == "europe_temperate" ~ "Temperate Europe",
                             region == "mark_twain_nf" ~ "Temperate USA",
                             region == "pacaya_samiria" ~ "Tropical S America",
                             region == "danum_valley" ~ "Tropical Borneo"))

# for month with warmest max temps
df2.tmax = mm.tmax %>% 
  mutate(hr = hour(obs_time), yday = yday(obs_time)) %>% 
  group_by(yday, reqhgt, canopy_hgt, shot_num, region) %>% 
  summarise(tair.min = min(tair),
            tair.max = max(tair),
            tleaf.min = min(tleaf),
            tleaf.max = max(tleaf),
            rh.min = min(relhum),
            windspeed.mean = mean(windspeed)) %>% 
  group_by(reqhgt, canopy_hgt, shot_num, region) %>% 
  summarise(tair.min = mean(tair.min),
            tair.max = mean(tair.max),
            tleaf.min = mean(tleaf.min),
            tleaf.max = mean(tleaf.max),
            rh.min = mean(rh.min),
            windspeed.mean = mean(windspeed.mean)) %>% 
  group_by(shot_num, region) %>% 
  mutate(tair.min.ground = tair.min[reqhgt==0.15],
         tair.max.ground = tair.max[reqhgt == 0.15],
         tleaf.min.ground = tleaf.min[reqhgt == 0.15],
         tleaf.max.ground = tleaf.max[reqhgt == 0.15],
         rh.min.ground = rh.min[reqhgt == 0.15],
         windspeed.mean.ground = windspeed.mean[reqhgt==0.15],
         tair.min.dif = tair.min - tair.min.ground,
         tair.max.dif = tair.max - tair.max.ground,
         tleaf.min.dif = tleaf.min - tleaf.min.ground,
         tleaf.max.dif = tleaf.max - tleaf.max.ground,
         rh.min.dif = rh.min - rh.min.ground,
         windspeed.mean.dif = windspeed.mean - windspeed.mean.ground,
         relhgt = reqhgt/canopy_hgt,
         frelhgt = cut(relhgt, breaks = seq(0,1,0.1), include_lowest = T, right = F)) %>% 
  group_by(frelhgt, region) %>% 
  summarise(tair.min.dif = median(tair.min.dif),
            tair.min = median(tair.min),
            tair.max.dif = median(tair.max.dif),
            tair.max = median(tair.max),
            tleaf.min.dif = median(tleaf.min.dif),
            tleaf.min = median(tleaf.min),
            tleaf.max.dif = median(tleaf.max.dif),
            tleaf.max = median(tleaf.max),
            rh.min.dif = median(rh.min.dif),
            rh.min = median(rh.min),
            windspeed.mean.dif = median(windspeed.mean.dif),
            windspeed.mean = median(windspeed.mean)) %>% 
  separate(frelhgt, into = c("minbin", "maxbin"), remove = F, sep = ",") %>% 
  mutate(minbin = as.numeric(gsub("\\[", "", minbin)),
         maxbin = as.numeric(gsub(")", "", maxbin)),
         lat = case_when(grepl("boreal", region) ~ "boreal",
                         region == "europe_temperate" ~ "temperate",
                         region == "mark_twain_nf" ~ "temperate",
                         region == "danum_valley" ~ "tropical",
                         region == "pacaya_samiria" ~ "tropical"),
         region2 = case_when(region == "canada_boreal" ~ "Boreal Canada",
                             region == "russia_boreal" ~ "Boreal Russia",
                             region == "europe_temperate" ~ "Temperate Europe",
                             region == "mark_twain_nf" ~ "Temperate USA",
                             region == "pacaya_samiria" ~ "Tropical S America",
                             region == "danum_valley" ~ "Tropical Borneo"))



# climate on the ground ---------------------------------------------------

# calculate the average climate on the ground
df.ground.tmin = mm.tmin %>% 
  filter(reqhgt == 0.15) %>% 
  mutate(hr = hour(obs_time), yday = yday(obs_time)) %>% 
  group_by(yday, reqhgt, canopy_hgt, shot_num, region) %>% 
  summarise(tair.min = min(tair),
            tair.max = max(tair),
            tleaf.min = min(tleaf),
            tleaf.max = max(tleaf),
            rh.min = min(relhum),
            windspeed.mean = mean(windspeed)) %>% 
  group_by(reqhgt, canopy_hgt, shot_num, region) %>% 
  summarise(tair.min = mean(tair.min),
            tair.max = mean(tair.max),
            tleaf.min = mean(tleaf.min),
            tleaf.max = mean(tleaf.max),
            rh.min = mean(rh.min),
            windspeed.mean = mean(windspeed.mean)) %>% 
  group_by(reqhgt, region) %>% 
  summarise(tair.min = mean(tair.min),
            tair.max = mean(tair.max),
            tleaf.min = mean(tleaf.min),
            tleaf.max = mean(tleaf.max),
            rh.min = mean(rh.min),
            windspeed.mean = mean(windspeed.mean)) %>% 
  mutate(region2 = case_when(region == "canada_boreal" ~ "Boreal Canada",
                              region == "russia_boreal" ~ "Boreal Russia",
                              region == "europe_temperate" ~ "Temperate Europe",
                              region == "mark_twain_nf" ~ "Temperate USA",
                              region == "pacaya_samiria" ~ "Tropical S America",
                              region == "danum_valley" ~ "Tropical Borneo"))

df.ground.tmax = mm.tmax %>% 
  filter(reqhgt == 0.15) %>% 
  mutate(hr = hour(obs_time), yday = yday(obs_time)) %>% 
  group_by(yday, reqhgt, canopy_hgt, shot_num, region) %>% 
  summarise(tair.min = min(tair),
            tair.max = max(tair),
            tleaf.min = min(tleaf),
            tleaf.max = max(tleaf),
            rh.min = min(relhum),
            windspeed.mean = mean(windspeed)) %>% 
  group_by(reqhgt, canopy_hgt, shot_num, region) %>% 
  summarise(tair.min = mean(tair.min),
            tair.max = mean(tair.max),
            tleaf.min = mean(tleaf.min),
            tleaf.max = mean(tleaf.max),
            rh.min = mean(rh.min),
            windspeed.mean = mean(windspeed.mean)) %>% 
  group_by(reqhgt, region) %>% 
  summarise(tair.min = mean(tair.min),
            tair.max = mean(tair.max),
            tleaf.min = mean(tleaf.min),
            tleaf.max = mean(tleaf.max),
            rh.min = mean(rh.min),
            windspeed.mean = mean(windspeed.mean)) %>% 
  mutate(region2 = case_when(region == "canada_boreal" ~ "Boreal Canada",
                             region == "russia_boreal" ~ "Boreal Russia",
                             region == "europe_temperate" ~ "Temperate Europe",
                             region == "mark_twain_nf" ~ "Temperate USA",
                             region == "pacaya_samiria" ~ "Tropical S America",
                             region == "danum_valley" ~ "Tropical Borneo"))



# * - plot latitudinal gradient in climate at different heights -----------

# ggplot(df2.tmax) +
#   geom_line(aes(x = region2, tair.max, color = region2, group = region2)) +
#   theme_bw()
# 
# ggplot(df2.tmin) +
#   geom_line(aes(x = region2, tair.min, color = minbin, group = region2), linewidth = 1) +
#   theme_bw()


# *- plot regions together and facet by variables -------------------------

region_order <- c("Boreal Canada", "Boreal Russia", "Temperate Europe", "Temperate USA", "Tropical Borneo", "Tropical S America")
spacing <- c(1, 1.2, 2, 2.2, 3, 3.2) 

df2.tmin %>% 
  arrange(region2) %>% 
  mutate(region2_num = spacing[match(region2, region_order)]) %>% 
  ggplot() +
  geom_line(aes(minbin, region2_num, color = tair.min.dif, group = region2_num), linewidth = 2) +
  coord_flip() +
  scale_color_continuous_divergingx(palette = "RdBu", rev = T, limits = c(-1,2.1),
                                    guide = guide_colorbar(title = "Tdif (\u00b0C)")) +
  scale_y_continuous(breaks = spacing, labels = region_order) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(color = "black", fill = "gray80"))
  

# plot individual variables -----------------------------------------------

# min during the coldest month
library(ggnewscale)
ggplot() +
  geom_line(data = df2.tmin, aes(minbin, tair.min.dif, color = tair.min.dif), linewidth = 2) +
  coord_flip() +
  scale_color_continuous_divergingx(palette = "RdBu", rev = T,
                                    guide = guide_colorbar(title = "Tdif (\u00b0C)")) +
  ggnewscale::new_scale_color() +
  geom_vline(data = df.ground.tmin, aes(xintercept = 0,color = tair.min), linewidth = 5) +
  scale_color_continuous_divergingx(palette = "Spectral", rev = T,
                                    guide = guide_colorbar(theme = theme(legend.title = element_blank(),
                                                                         legend.key.width = unit(60, "mm")),
                                                           position = "bottom")) +
  scale_x_continuous("Relative Height", breaks = seq(0,1,0.2), limits = c(0,0.95), expand = c(0,0)) +
  scale_y_continuous("Minimum temperature difference to the ground (\u00b0C)") +
  facet_wrap(~region2, nrow = 1) +
  theme_bw() +
  theme(panel.grid = element_blank())

# plot max temp of warmest month
ggplot() +
  geom_line(data = df2.tmax, aes(minbin, tair.max.dif, color = tair.max.dif), linewidth = 2) +
  coord_flip() +
  scale_color_continuous_divergingx(palette = "RdBu", rev = T,
                                    guide = guide_colorbar(title = "Tdif (\u00b0C)")) +
  ggnewscale::new_scale_color() +
  geom_vline(data = df.ground.tmax, aes(xintercept = 0,color = tair.min), linewidth = 5) +
  scale_color_continuous_divergingx(palette = "Spectral", rev = T,
                                    guide = guide_colorbar(theme = theme(legend.title = element_blank(),
                                                                         legend.key.width = unit(60, "mm")),
                                                           position = "bottom"), ) +
  scale_x_continuous("Relative Height", breaks = seq(0,1,0.2), limits = c(0,0.95), expand = c(0,0)) +
  scale_y_continuous("Maximum temperature difference to the ground (\u00b0C)") +
  facet_wrap(~region2, nrow = 1) +
  theme_bw() +
  theme(panel.grid = element_blank())

# plot min RH of warmest month
ggplot() +
  geom_line(data = df2.tmax, aes(minbin, rh.min.dif, color = rh.min.dif), linewidth = 2) +
  coord_flip() +
  scale_color_continuous_divergingx(palette = "RdBu", rev = T,
                                    guide = guide_colorbar(title = "Tdif (\u00b0C)")) +
  ggnewscale::new_scale_color() +
  geom_vline(data = df.ground.tmax, aes(xintercept = 0,color = rh.min), linewidth = 5) +
  scale_color_continuous_divergingx(palette = "Spectral", rev = T,
                                    guide = guide_colorbar(theme = theme(legend.title = element_blank(),
                                                                         legend.key.width = unit(60, "mm")),
                                                           position = "bottom"), ) +
  scale_x_continuous("Relative Height", breaks = seq(0,1,0.2), limits = c(0,0.95), expand = c(0,0)) +
  scale_y_continuous("Relative humidity: difference to the ground (%)") +
  facet_wrap(~region2, nrow = 1) +
  theme_bw() +
  theme(panel.grid = element_blank())


# plot all variables together ---------------------------------------------

# plot everything together
ggplot() +
  
  geom_line(data = pad, aes(factor("PAVD"), minbin, color = pavd), linewidth = 3) +
  scale_color_continuous_sequential(palette = "Greens", rev = T,
                                    guide = guide_colorbar(title = "PAVD (m2/m3)")) +

  ggnewscale::new_scale_color() +
  geom_line(data = df2.tmin, aes(factor("TMIN"), minbin, color = tair.min.dif), linewidth = 3) +
  scale_color_continuous_sequential(palette = "PuBu", rev = F,
                                    guide = guide_colorbar(title = "TMINdif (\u00b0C)")) +
  ggnewscale::new_scale_color() +
  geom_point(data = df.ground.tmin, aes(x = factor("TMIN"), y = 0, color = tair.min), shape = "\u2014", size = 10) +
  scale_color_continuous_sequential(palette = "PuBu", rev = F,
                                    guide = guide_colorbar(title = "Ground TMIN (\u00b0C)", 
                                                           theme = theme(legend.key.width = unit(60, "mm"),
                                                                         legend.title.position = "top"),
                                                           position = "bottom")) +
  ggnewscale::new_scale_color() +
  geom_line(data = df2.tmax, aes(factor("TMAX"), minbin, color = tair.max.dif), linewidth = 3) +
  scale_color_continuous_sequential(palette = "OrRd", rev = T,
                                    guide = guide_colorbar(title = "TMAXdif (%)")) +
  ggnewscale::new_scale_color() +
  geom_point(data = df.ground.tmax, aes(x = factor("TMAX"), y = 0, color = tair.max), shape = "\u2014", size = 10) +
  scale_color_continuous_sequential(palette = "OrRd", rev = T,
                                    guide = guide_colorbar(title = "Ground TMAX (\u00b0C)", 
                                                           theme = theme(legend.key.width = unit(60, "mm"),
                                                                         legend.title.position = "top"),
                                                           position = "bottom")) +

  ggnewscale::new_scale_color() +
  geom_line(data = df2.tmax, aes(factor("RH"), minbin, color = rh.min.dif), linewidth = 3) +
  scale_color_continuous_sequential(palette = "BluYl", rev = T,
                                    guide = guide_colorbar(title = "RHdif (%)")) +
  ggnewscale::new_scale_color() +
  #geom_vline(data = df.ground.tmax, aes(xintercept = 0,color = rh.min), linewidth = 5) +
  geom_point(data = df.ground.tmax, aes(x = factor("RH"), y = 0, color = rh.min), shape = "\u2014", size = 10) +
  scale_color_continuous_sequential(palette = "BluYl", rev = T,
                                    guide = guide_colorbar(title = "Ground RH (%)", 
                                                           theme = theme(legend.key.width = unit(60, "mm"),
                                                                         legend.title.position = "top"),
                                                           position = "bottom")) +
  
  scale_y_continuous("Relative Height", breaks = seq(0,1,0.2), limits = c(0,0.95), expand = c(0.01,0)) +
  #scale_x_discrete("Relative humidity: difference to the ground (%)") +
  facet_wrap(~region2, nrow = 1) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        panel.background = element_blank())
ggsave("figures/main_figs/vertical_gradients/vertical_gradients.svg",
       height = 200, width = 300, units = "mm")
ggsave("figures/main_figs/vertical_gradients/vertical_gradients.png",
       height = 200, width = 300, units = "mm")



# plot lines per shot number ----------------------------------------------



df2 = mm.tmin %>% 
  mutate(doy = yday(obs_time)) %>% 
  group_by(shot_num, region, canopy_hgt, reqhgt, doy) %>% 
  summarise(tair.min = min(tair)) %>% 
  group_by(shot_num, region, canopy_hgt, reqhgt) %>% 
  summarise(tair.min = mean(tair.min)) 

df2max = mm.tmax %>% 
  mutate(doy = yday(obs_time)) %>% 
  group_by(shot_num, region, canopy_hgt, reqhgt, doy) %>% 
  summarise(tair.max = max(tair), rhmin = min(relhum)) %>% 
  group_by(shot_num, region, canopy_hgt, reqhgt) %>% 
  summarise(tair.min = mean(tair.max), rhmin = mean(rhmin)) 

df2 %>% 
  filter(region == "canada_boreal") %>% 
  ggplot(aes(reqhgt, tair.min, group = shot_num)) +
  geom_line() +
  coord_flip() +
  theme_classic()

df2 %>% 
  filter(region == "russia_boreal") %>% 
  ggplot(aes(reqhgt, tair.min, group = shot_num)) +
  geom_line() +
  coord_flip() +
  theme_classic()

df2 %>% 
  filter(region ==  "mark_twain_nf") %>% 
  ggplot(aes(reqhgt, tair.min, group = shot_num)) +
  geom_line() +
  coord_flip() +
  theme_classic()

tropics = df2max %>% 
  filter(region ==  "pacaya_samiria") %>% 
  ggplot(aes(reqhgt, rhmin, group = shot_num)) +
  geom_line() +
  scale_y_continuous(limits = c(20,90)) +
  scale_x_continuous(limits = c(0,50)) +
  coord_flip() +
  theme_classic()

temperate = df2max %>% 
  filter(region ==  "mark_twain_nf") %>% 
  ggplot(aes(reqhgt, rhmin, group = shot_num)) +
  geom_line() +
  scale_y_continuous(limits = c(20,90)) +
  scale_x_continuous(limits = c(0,50)) +
  coord_flip() +
  theme_classic()


boreal = df2max %>% 
  filter(region ==  "canada_boreal") %>% 
  ggplot(aes(reqhgt, rhmin, group = shot_num)) +
  geom_line() +
  scale_y_continuous(limits = c(20,90)) +
  scale_x_continuous(limits = c(0,50)) +
  coord_flip() +
  theme_classic()

tropics2 = df2max %>% 
  filter(region ==  "danum_valley") %>% 
  ggplot(aes(reqhgt, rhmin, group = shot_num)) +
  geom_line() +
  scale_y_continuous(limits = c(20,90)) +
  scale_x_continuous(limits = c(0,50)) +
  coord_flip() +
  theme_classic()

temperate2 = df2max %>% 
  filter(region ==  "europe_temperate") %>% 
  ggplot(aes(reqhgt, rhmin, group = shot_num)) +
  geom_line() +
  scale_y_continuous(limits = c(20,90)) +
  scale_x_continuous(limits = c(0,50)) +
  coord_flip() +
  theme_classic()


boreal2 = df2max %>% 
  filter(region ==   "russia_boreal") %>% 
  ggplot(aes(reqhgt, rhmin, group = shot_num)) +
  geom_line() +
  scale_y_continuous(limits = c(20,90)) +
  scale_x_continuous(limits = c(0,50)) +
  coord_flip() +
  theme_classic()


library(patchwork)
rhplt = boreal / temperate /tropics

design = "
adg
beh
cfi"

pavdplt[[1]] + pavdplt[[2]] + pavdplt [[3]] + 
  rhplt[[1]] + rhplt[[2]] + rhplt[[3]] +
  boreal2 + temperate2 + tropics2 + plot_layout(design = design, byrow = T)

# max temp during the coldest month
df2 %>% 
  ggplot(aes(minbin, tair.max.dif, color = tair.max.dif)) +
  geom_line(linewidth = 2) +
  coord_flip() +
  scale_color_continuous_divergingx(palette = "Spectral", rev = T) +
  facet_wrap(~region, nrow = 1)

# min rh during the coldest month
df2 %>% 
  ggplot(aes(minbin, rh.min.dif, color = rh.min.dif)) +
  geom_line(linewidth = 2) +
  coord_flip() +
  scale_color_continuous_divergingx(palette = "Spectral", rev = T) +
  facet_wrap(~region, nrow = 1)
 
# min leaf temp during the coldest month
df2 %>% 
  ggplot(aes(minbin, tleaf.min.dif, color = tleaf.min.dif)) +
  geom_line(linewidth = 2) +
  coord_flip() +
  scale_color_continuous_divergingx(palette = "Spectral", rev = T) +
  facet_wrap(~region, nrow = 1)

# max leaf temp during the coldest month
df2 %>% 
  ggplot(aes(minbin, tleaf.max.dif, color = tleaf.max.dif)) +
  geom_line(linewidth = 2) +
  coord_flip() +
  scale_color_continuous_divergingx(palette = "Spectral", rev = T) +
  facet_wrap(~region, nrow = 1)

# average windspeed during the coldest month
# min leaf temp during the coldest month
df2 %>% 
  ggplot(aes(minbin, windspeed.mean.dif, color = windspeed.mean.dif)) +
  geom_line(linewidth = 2) +
  coord_flip() +
  scale_color_continuous_divergingx(palette = "Spectral", rev = T) +
  facet_wrap(~region, nrow = 1)









df %>% filter(region == "europe_temperate") %>% 
  ggplot(aes(region, minbin, fill = tair.min)) +
  scale_fill_continuous_sequential(palette = "Blues") +
  geom_tile()

df %>% filter(region == "mark_twain_nf") %>% 
  ggplot(aes(region, minbin, fill = tair.min)) +
  scale_fill_continuous_sequential(palette = "Blues") +
  geom_tile()

df %>% filter(region == "mark_twain_nf") %>%
  ggplot(aes(minbin, tair.min, color = tair.min)) +
  scale_fill_continuous_sequential(palette = "Blues") +
  geom_line(linewidth = 2) +
  coord_flip()

df %>% filter(region == "europe_temperate") %>%
  ggplot(aes(minbin, tair.min, color = tair.min)) +
  scale_fill_continuous_sequential(palette = "Blues") +
  geom_line(linewidth = 2) +
  coord_flip()

df %>% filter(region == "canada_boreal") %>%
  ggplot(aes(minbin, tair.min, color = tair.min)) +
  scale_fill_continuous_sequential(palette = "Blues") +
  geom_line(linewidth = 2) +
  coord_flip()

df %>% filter(region == "russia_boreal") %>%
  ggplot(aes(minbin, tair.min, color = tair.min)) +
  scale_fill_continuous_sequential(palette = "Blues") +
  geom_line(linewidth = 2) +
  coord_flip()

df %>% filter(region == "pacaya_samiria") %>%
  ggplot(aes(minbin, tair.min, color = tair.min)) +
  scale_fill_continuous_sequential(palette = "Blues") +
  geom_line(linewidth = 2) +
  coord_flip()

df %>% filter(region == "danum_valley") %>%
  ggplot(aes(minbin, tair.min, color = tair.min)) +
  scale_fill_continuous_sequential(palette = "Blues") +
  geom_line(linewidth = 2) +
  coord_flip()

# Load in data for Pacaya Samiria gedi points
indat = readRDS("D:/projects/verticality_future/data/derivative_data/micropoint_prep/prep_micropoint_pacaya_samiria.rds")

ps.geditmin = indat$gedisamp.mcold
ps.geditmax = indat$gedisamp.mhot
ps.microin = indat$microin

# takes about 15 min to run 100 points
tic()
tropical.ps.tmin = run_micropoint(
  tme = ps.microin$climr$`2019`$obs_time,
  gedi = ps.geditmin,
  climr = ps.microin$climr,
  vegp = ps.microin$vegp,
  soil = ps.microin$soil,
  elev = ps.microin$elev,
  asp = ps.microin$asp,
  slp = ps.microin$slp,
  dtmc = ps.microin$dtmc,
  method = "temporal_month",
  reqhgts = NA,
  vertpai_method = "pai",
  fout = "results/microclim_models/pacaya_samiria/mintemp_cold.csv"
)
toc()

df = tropical.ps.tmin %>%
  mutate(mday = mday(obs_time)) %>%
  group_by(shot_num, reqhgt, mday) %>%
  summarise(tair.minday = min(tair)) %>%
  group_by(shot_num, reqhgt) %>%
  summarise(tmin_cold = mean(tair.minday))
  # group_by(reqhgt) %>%
  # summarise(tmin_cold.mean = mean(tmin_cold), tmin_cold.sd = sd(tmin_cold))
  #

ggplot(df, aes(tmin_cold, reqhgt, group = as.factor(shot_num))) +
  geom_line() +
  theme_bw()

# Boreal: Canada
indat = readRDS("D:/projects/verticality_future/data/derivative_data/micropoint_prep/prep_micropoint_canada_boreal.rds")

geditmin = indat$gedisamp.mcold
geditmax = indat$gedisamp.mhot
microin = indat$microin

tic()
borealCanada.tmin = run_micropoint(
  tme = microin$climr$`2019`$obs_time,
  gedi = geditmin,
  climr = microin$climr,
  vegp = microin$vegp,
  soil = microin$soil,
  elev = microin$elev,
  asp = microin$asp,
  slp = microin$slp,
  dtmc = microin$dtmc,
  method = "temporal_month",
  reqhgts = NA,
  vertpai_method = "pai",
  fout = "results/microclim_models/canada_boreal/mintemp_cold.csv"
)
toc()

df = borealCanada.tmin %>%
  mutate(mday = mday(obs_time)) %>%
  group_by(shot_num, reqhgt, mday, canopy_hgt) %>%
  summarise(tair.minday = min(tair)) %>%
  group_by(shot_num, reqhgt, canopy_hgt) %>%
  summarise(tmin_cold = median(tair.minday))
# group_by(reqhgt) %>%
# summarise(tmin_cold.mean = mean(tmin_cold), tmin_cold.sd = sd(tmin_cold))
#

dfleaf = borealCanada.tmin %>%
  mutate(mday = mday(obs_time)) %>%
  group_by(shot_num, reqhgt, mday, canopy_hgt) %>%
  summarise(tleaf.minday = min(tleaf)) %>%
  group_by(shot_num, reqhgt, canopy_hgt) %>%
  summarise(tmin_cold = median(tleaf.minday))

df %>% 
  filter(tmin_cold < -35) %>% 
  mutate(relhgt = reqhgt/canopy_hgt) %>% 
  ggplot(aes(tmin_cold, reqhgt, group = as.factor(shot_num), color = tmin_cold)) +
  geom_line() +
  theme_bw()

dfleaf %>% 
  mutate(relhgt = reqhgt/canopy_hgt) %>% 
  ggplot(aes(tmin_cold, reqhgt, group = as.factor(shot_num), color = tmin_cold)) +
  geom_line() +
  theme_bw()

df %>% 
  filter(tmin_cold < -35) %>% 
  ungroup() %>% 
  mutate(relhgt = reqhgt/canopy_hgt,
         hgroup = cut_width(relhgt, width = 0.1, boundary = 0.1)) %>% 
  group_by(hgroup) %>%
  summarise(tmin_cold = median(tmin_cold)) %>%
  mutate(relhgt = seq(0,0.9,0.1),
         x = factor(1)) %>% 
  ggplot(aes(x, relhgt, fill = tmin_cold)) +
  geom_tile() +
  scale_fill_continuous_sequential("Blues") +
  theme_bw()

df2 = borealCanada.tmin %>% 
  mutate(mday = mday(obs_time)) %>% 
  group_by(shot_num, reqhgt, mday, canopy_hgt) %>% 
  summarise(tair.minday = min(tair))

df2 %>% 
  filter(reqhgt == 0.15) %>% 
  mutate(h = ifelse(reqhgt == 0.15, "ground", "canopy")) %>% 
  ggplot(aes(mday, tair.minday, color = h, group = shot_num)) +
  geom_line()




