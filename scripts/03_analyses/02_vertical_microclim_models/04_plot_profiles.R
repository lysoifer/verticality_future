# Plot boxplots of vertical microclimate profiles
# Then consider running an NMDS to look at bioclim at different latitudes and iin canopy and ground

library(ggthemes)
library(dplyr)
library(ggplot2)
library(foreach)
library(data.table)
library(tidyr)

# Load in modeling data
micromods = list.files("results/microclim_models/", recursive = T, full.names = T)
micromods.tmin = grep("mintemp", micromods, value = T) # points for the coldest month of the year
micromods.tmax = grep("tmax", micromods, value = T) # points for the warmest month of the year
indat = list.files("D:/projects/verticality_future/data/derivative_data/micropoint_prep/", full.names = T) # gedi points used for models

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
  cnames = c("shot_num", "rh_100_a0", "lat_lm_a0", "lon_lm_a0", grep("[0-9]$", cnames, value = T))
  gedi = gedi[, ..cnames]
  gedi[gedi==0] <- NA # convert zeros to NAs
  gedi = gedi %>% 
    pivot_longer(cols = 5:ncol(gedi), names_to = "height", values_to = "pavd") %>% 
    drop_na()
  
  nm = basename(indat[f])
  nm = gsub("prep_micropoint_", "", nm)
  nm = gsub(".rds", "", nm)
  gedi$region = nm
  gedi
}

# summarize PAVD by relhgt
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
         frelhgt = cut(relhgt, breaks = seq(0,1,0.1), include_lowest = T, right = T)) %>%  # break points
  filter(relhgt <= 1) %>% 
  group_by(frelhgt, region, shot_num, lat_lm_a0, lon_lm_a0, h) %>% 
  summarize(pavd = median(pavd, na.rm = T)) %>% # mean pavd per relhgt 
  separate(frelhgt, into = c("minbin", "maxbin"), remove = F, sep = ",") %>% 
  mutate(minbin = as.numeric(gsub("\\(", "", minbin)),
         maxbin = as.numeric(gsub("\\]", "", maxbin)),
         region2 = case_when(region == "canada_boreal" ~ "Boreal Canada",
                             region == "russia_boreal" ~ "Boreal Russia",
                             region == "europe_temperate" ~ "Temperate Europe",
                             region == "mark_twain_nf" ~ "Temperate USA",
                             region == "pacaya_samiria" ~ "Tropical S America",
                             region == "danum_valley" ~ "Tropical Borneo"))


# data frame for month with minimum cold temperatures

# Average min temp per day in the coldest month
# calculate min temp per day for each point
# take mean of min temps across days for each point
df2.tmin = mm.tmin %>% 
  mutate(doy = yday(obs_time), 
         relhgt = reqhgt/canopy_hgt,
         frelhgt = cut(relhgt, breaks = seq(0,1,0.1), include_lowest = T, right = T)) %>% 
  group_by(shot_num, region, canopy_hgt, frelhgt, relhgt, reqhgt, doy) %>% 
  summarise(tair.min = min(tair)) %>% 
  group_by(shot_num, region, frelhgt, relhgt, canopy_hgt, reqhgt) %>% 
  summarise(tair.min = mean(tair.min)) %>% 
  mutate(f_lat = case_when(grepl("boreal", region)~"Boreal",
                           grepl("mark_twain|europe", region)~"Temperate",
                           grepl("danum|pacaya", region)~"Tropical"))

# Average max temp and minRH per day in the warmest month
# calculate max temp and min rh per day for each point
# take mean of max temps and min rh across days for each point
df2.tmax = mm.tmax %>% 
  mutate(doy = yday(obs_time), 
         relhgt = reqhgt/canopy_hgt,
         frelhgt = cut(relhgt, breaks = seq(0,1,0.1), include_lowest = T, right = T)) %>% 
  group_by(shot_num, region, canopy_hgt, frelhgt, relhgt, reqhgt, doy) %>% 
  summarise(tair.max = max(tair), rhmin = min(relhum)) %>% 
  group_by(shot_num, region, frelhgt, relhgt, canopy_hgt, reqhgt) %>% 
  summarise(tair.max = mean(tair.max), rhmin = mean(rhmin)) %>% 
  mutate(f_lat = case_when(grepl("boreal", region)~"Boreal",
                           grepl("mark_twain|europe", region)~"Temperate",
                           grepl("danum|pacaya", region)~"Tropical"))

# summary dataframes
df2.tmax.summ = df2.tmax %>% 
  ungroup() %>% 
  group_by(f_lat, frelhgt) %>% 
  summarise(rhmin.med = median(rhmin), tair.max.med = median(tair.max),
            rhmin.sd = sd(rhmin), tair.max.sd = sd(tair.max))

df2.tmin.summ = df2.tmin %>% 
  ungroup() %>% 
  group_by(f_lat, frelhgt) %>% 
  summarise(tair.min.med = median(tair.min), tair.min.sd = sd(tair.min))

pad.summ = pad2 %>% 
  mutate(f_lat = case_when(grepl("boreal", region)~"Boreal",
                           grepl("mark_twain|europe", region)~"Temperate",
                           grepl("danum|pacaya", region)~"Tropical")) %>% 
  group_by(frelhgt, f_lat) %>% 
  summarize(val.median = median(pavd), val.sd = sd(pavd)) %>% 
  mutate(var = "PAVD (m2/m3)")


# VERTICAL PROFILES ----------------------------------------------------------------
pal = c(Boreal = "#56B4E9", Temperate = "#CC79A7", Tropical = "#009E73")
thm = theme(legend.title = element_blank())

# combine dataframes for plotting
d1 = df2.tmin.summ %>% 
  rename(val.median = tair.min.med,
         val.sd = tair.min.sd) %>% 
  mutate(var = "Minimum Temperature (\u00b0C)")

d2 = df2.tmax.summ %>% 
  dplyr::select(f_lat, frelhgt, tair.max.med, tair.max.sd) %>% 
  rename(val.median = tair.max.med,
         val.sd = tair.max.sd) %>% 
  mutate(var = "Maximum Temperature (\u00b0C)")

# rh during warmest month
d3 = df2.tmax.summ %>% 
  dplyr::select(f_lat, frelhgt, rhmin.med, rhmin.sd) %>% 
  rename(val.median = rhmin.med,
         val.sd = rhmin.sd) %>% 
  mutate(var = "Minimum Relative Humidity (%)")

dfplt = bind_rows(d1, d2, d3, pad.summ)

ggplot(dfplt) +
  #geom_boxplot(data = df2.tmax, aes(x = frelhgt, y = rhmin)) +
  geom_pointrange(aes(x = frelhgt, y = val.median,
                                            ymin = val.median-val.sd,
                                            ymax = val.median+val.sd,
                                            color = f_lat),
                  position = position_dodge(width = 0.2)) +
  geom_line(aes(x = frelhgt, y = val.median, group = f_lat, color = f_lat),
            linewidth = 1) +
  #scale_y_continuous("Minimum Relative Humidity (%)") +
  scale_x_discrete("Relative Height in the Canopy") +
  scale_color_manual(values = pal) +
  coord_flip() +
  facet_wrap(~var, nrow = 1, scales = "free_x", strip.position = "bottom") +
  theme_bw() +
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank())




# VARIABILITY -------------------------------------------------------------

# issue because min and max temps are not in the same location

ggplot(df2.tmax) +
  geom_point(aes(tair.max, rhmin, color = f_lat), alpha = 0.5) +
  theme_classic()

df2.tmax %>% 
  filter(reqhgt == 0.15 | relhgt >= 0.8) %>% 
  mutate(f_hgt = ifelse(reqhgt == 0.15, "ground", "canopy"),
         cat = paste(f_lat, f_hgt, sep = "_")) %>% 
  ggplot() +
  geom_point(aes(tair.max, rhmin, color = cat), alpha = 0.5) +
  theme_classic()

df2.tmin %>% 
  filter(reqhgt == 0.15 | relhgt >= 0.8) %>% 
  mutate(f_hgt = ifelse(reqhgt == 0.15, "ground", "canopy"),
         cat = paste(f_lat, f_hgt, sep = "_")) %>% 
  ggplot() +
  geom_point(aes(tair.min, rhmin, color = cat), alpha = 0.5) +
  theme_classic()
  
  
# try summaryizing by a slightly coarser resolution
pad2 = pad2 %>% 
  mutate(lat2 = cut(lat_lm_a0, breaks = seq(-10,60,0.05), include_lowest = T, right = T),
         lon2 = cut(lon_lm_a0, breaks = seq(-92,140,0.05), include_lowest = T, right = T))

# summarize by 0.05degree grid cell
dfvarplt.tmin = inner_join(pad2, df2.tmin, by = c("shot_num", "region", "frelhgt", "h" = "reqhgt")) %>% 
  filter(h == 0.15 | relhgt > 0.8) %>% 
  group_by(frelhgt, minbin, maxbin, region, region2, lat2, lon2, f_lat) %>% 
  summarise(pavd = mean(pavd), tair.min = median(tair.min))

dfvarplot.tmax = inner_join(pad2, df2.tmax, by = c("shot_num", "region", "frelhgt", "h" = "reqhgt")) %>% 
  filter(h == 0.15 | relhgt > 0.8) %>% 
  group_by(frelhgt, minbin, maxbin, region, region2, lat2, lon2, f_lat) %>% 
  summarise(pavd = mean(pavd), tair.max = median(tair.max), rhmin = median(rhmin))

dfvarplt = inner_join(dfvarplt.tmin, dfvarplot.tmax, by = c("lat2", "lon2", "frelhgt", "minbin", "maxbin", "region", "region2", "f_lat"))

dfvarplt %>% 
  mutate(f_hgt = as.character(frelhgt),
         f_hgt = ifelse(f_hgt == "(0,0.1]", "ground", "canopy")) %>% 
  mutate(cat = paste(region, f_hgt, sep = "_")) %>% 
  ggplot() +
  geom_point(aes(tair.min, rhmin, color = cat)) +
  theme_classic()



fig <- dfvarplt %>% 
  mutate(f_hgt = as.character(frelhgt),
         f_hgt = ifelse(f_hgt == "(0,0.1]", "ground", "canopy")) %>% 
  mutate(cat = paste(region, f_hgt, sep = "_")) %>% 
  plot_ly(x = ~tair.min, y = ~tair.max, z = ~rhmin, color = ~cat)
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'Tmin'),
                                   yaxis = list(title = 'Tmax'),
                                   zaxis = list(title = 'RHmin')))
fig




  
