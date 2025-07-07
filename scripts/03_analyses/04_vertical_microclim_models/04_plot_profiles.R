# Plot boxplots of vertical microclimate profiles
# Then consider running an NMDS to look at bioclim at different latitudes and iin canopy and ground

library(ggthemes)
library(dplyr)
library(ggplot2)
library(foreach)
library(data.table)
library(tidyr)

thm = theme(legend.text = element_text(size = 4),
            legend.title = element_text(size = 6),
            axis.text = element_text(size = 4),
            axis.title = element_text(size = 6),
            strip.text = element_text(size = 6))

# Load in modeling data
micromods = list.files("results/microclim_models2/", recursive = T, full.names = T)
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
# calculate min temp and temp range per day for each point
# take mean of min temps across days for each point
df2.tmin = mm.tmin %>% 
  mutate(doy = yday(obs_time), 
         relhgt = reqhgt/canopy_hgt,
         frelhgt = cut(relhgt, breaks = seq(0,1,0.1), include_lowest = T, right = T)) %>% 
  group_by(shot_num, region, canopy_hgt, frelhgt, relhgt, reqhgt, doy) %>% 
  summarise(tair.min = min(tair), tair.diu = max(tair) - min(tair)) %>% 
  group_by(shot_num, region, frelhgt, relhgt, canopy_hgt, reqhgt) %>% 
  summarise(tair.min = mean(tair.min), tair.diu = mean(tair.diu)) %>% 
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
  summarise(tair.min.med = median(tair.min), tair.min.sd = sd(tair.min),
            tair.tdiu.med = median(tair.diu), tair.tdiu.sd = sd(tair.diu))

pad.summ = pad2 %>% 
  mutate(f_lat = case_when(grepl("boreal", region)~"Boreal",
                           grepl("mark_twain|europe", region)~"Temperate",
                           grepl("danum|pacaya", region)~"Tropical")) %>% 
  group_by(frelhgt, f_lat) %>% 
  summarize(val.median = median(pavd), val.sd = sd(pavd)) %>% 
  mutate(var = "PAVD (m2/m3)")


# VERTICAL PROFILES ----------------------------------------------------------------

pal = c(Boreal =  "#56B4E9", Temperate = "#E08A00", Tropical = "#008234")
thm = theme(legend.title = element_blank())

# combine dataframes for plotting
d1 = df2.tmin.summ %>% 
  dplyr::select(f_lat, frelhgt, tair.min.med, tair.min.sd) %>% 
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

d4 = df2.tmin.summ %>% 
  dplyr::select(f_lat, frelhgt, tair.tdiu.med, tair.tdiu.sd) %>% 
  rename(val.median = tair.tdiu.med,
         val.sd = tair.tdiu.sd) %>% 
  mutate(var = "Diurnal temperature range\ncoldest month (\u00b0C)")

dfplt = bind_rows(d1, d2, d3, d4, pad.summ)

dfplt = dfplt %>% 
  mutate(minbin = str_split_i(frelhgt, ",", 1),
         minbin = gsub("\\(", "", minbin),
         minbin = as.numeric(minbin))

profileplot = ggplot(dfplt) +
  #geom_boxplot(data = df2.tmax, aes(x = frelhgt, y = rhmin)) +
  # geom_pointrange(aes(x = frelhgt, y = val.median,
  #                                           ymin = val.median-val.sd,
  #                                           ymax = val.median+val.sd,
  #                                           color = f_lat),
  #                 position = position_dodge(width = 0.2)) +
  geom_ribbon(aes(x = minbin,
                  ymin = val.median-val.sd,
                  ymax = val.median+val.sd,
                  fill = f_lat), alpha = 0.2) + 
  geom_line(aes(x = minbin, y = val.median, group = f_lat, color = f_lat),
           linewidth = 0.5) +
  #scale_y_continuous("Minimum Relative Humidity (%)") +
  scale_x_continuous("Relative Height in the Canopy") +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  coord_flip() +
  facet_wrap(~var, nrow = 1, scales = "free_x", strip.position = "bottom") +
  theme_bw() +
  thm +
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(2, "mm"))


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

# plot tmin coldest month and rhmin warmest month
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




# plot offsets ------------------------------------------------------------

df2.tmin.offset = df2.tmin %>% 
  group_by(shot_num) %>% 
  mutate(min_tair.min = min(tair.min), max_tair.min = max(tair.min),
         offset_tair.min = max_tair.min - min_tair.min) %>% 
  dplyr::select(shot_num, canopy_hgt, f_lat, region, offset_tair.min) %>% 
  distinct()


df2.tmax.offset = df2.tmax %>% 
  group_by(shot_num) %>% 
  mutate(min_tair.max = min(tair.max), max_tair.max = max(tair.max),
         offset_tair.max = max_tair.max - min_tair.max,
         min_rhmin = min(rhmin), max_rhmin = max(rhmin),
         offset_rhmin = max_rhmin - min_rhmin) %>% 
  dplyr::select(shot_num, canopy_hgt, f_lat, region, offset_tair.max, offset_rhmin) %>% 
  distinct()

p1 = ggplot() +
  geom_boxplot(data = df2.tmax.offset, aes(f_lat, offset_rhmin, fill = region)) +
  theme_bw()

p2 = ggplot() +
  geom_boxplot(data = df2.tmax.offset, aes(f_lat, offset_tair.max, fill = region)) +
  theme_bw()


p3 = ggplot() +
  geom_boxplot(data = df2.tmin.offset, aes(f_lat, offset_tair.min, fill = region)) +
  theme_bw()


# PLOT STRATIFICATION -----------------------------------------------------

d = mm.tmax %>% 
  filter(region == "mark_twain_nf") %>% 
  mutate(obs_time2 = case_when(region == "danum_valley" ~ obs_time+8*3600,
                               region == "pacaya_samiria" ~ obs_time-5*3600,
                               region == "europe_temperate" ~ obs_time+2*3600,
                               region == "mark_twain_nf" ~ obs_time-5*3600,
                               region == "canada_boreal" ~ obs_time-4*3600,
                               region == "russia_boreal" ~ obs_time+10*3600),
         hour = hour(obs_time2)) %>% 
  group_by(shot_num, hour, canopy_hgt, reqhgt, region) %>% 
  summarise(tair = mean(tair), tleaf = mean(tleaf), relhum = mean(relhum)) %>% 
  ggplot() +
  geom_jitter(aes(x = tair, y = reqhgt/canopy_hgt), pch = ".")

d = mm.tmax %>% 
  filter(region == "pacaya_samiria") %>% 
  mutate(hour = hour(obs_time)) %>% 
  group_by(shot_num, hour, canopy_hgt, reqhgt, region) %>% 
  summarise(tair = mean(tair), tleaf = mean(tleaf), relhum = mean(relhum)) %>% 
  ggplot() +
  geom_jitter(aes(x = tair, y = reqhgt/canopy_hgt), pch = ".")

d = mm.tmax %>% 
  mutate(obs_time2 = case_when(region == "danum_valley" ~ obs_time+8*3600,
                               region == "pacaya_samiria" ~ obs_time-5*3600,
                               region == "europe_temperate" ~ obs_time+2*3600,
                               region == "mark_twain_nf" ~ obs_time-5*3600,
                               region == "canada_boreal" ~ obs_time-4*3600,
                               region == "russia_boreal" ~ obs_time+10*3600),
         hour = hour(obs_time2),
         relhgt = reqhgt/canopy_hgt,
         frelhgt = cut(relhgt, breaks = seq(0,1,0.1), include_lowest = T, right = T)) %>%  # break points
  filter(relhgt <= 1) %>% 
  group_by(hour, frelhgt, region) %>% 
  summarize(tair = mean(tair), relhum = mean(relhum), tleaf = mean(tleaf))%>% 
  mutate(site = case_when(region == "danum_valley" ~ "Tropical Asia",
                          region == "pacaya_samiria" ~ "Tropical S America",
                          region == "canada_boreal" ~ "Boreal N America",
                          region == "russia_boreal" ~ "Boreal Asia",
                          region == "mark_twain_nf" ~ "Temperate N America",
                          region == "europe_temperate" ~ "Temperate Europe") )


p1 = ggplot(d) +
  geom_line(aes(hour, tair, color = frelhgt)) +
  facet_wrap(~region)

ggplot(d) +
  geom_line(aes(hour, tleaf, color = frelhgt)) +
  facet_wrap(~region)

ggplot(d) +
  geom_line(aes(hour, relhum, color = frelhgt)) +
  facet_wrap(~region)


d2 = mm.tmax %>% 
  mutate(obs_time2 = case_when(region == "danum_valley" ~ obs_time+8*3600,
                               region == "pacaya_samiria" ~ obs_time-5*3600,
                               region == "europe_temperate" ~ obs_time+2*3600,
                               region == "mark_twain_nf" ~ obs_time-5*3600,
                               region == "canada_boreal" ~ obs_time-4*3600,
                               region == "russia_boreal" ~ obs_time+10*3600),
         hour = hour(obs_time2),
         relhgt = reqhgt/canopy_hgt,
         frelhgt = cut(relhgt, breaks = seq(0,1,0.1), include_lowest = T, right = T)) %>%  # break points
  filter(relhgt <= 1) %>% 
  group_by(hour, frelhgt, region) %>% 
  summarize(tair = mean(tair), relhum = mean(relhum), tleaf = mean(tleaf)) %>% 
  group_by(frelhgt, region) %>% 
  summarize(tair.min = min(tair), tair.max = max(tair),
            relhum.min = min(relhum), relhum.max = max(relhum))

d2.1 = d2 %>% 
  mutate(x = 0)
d2.2 = d2 %>% 
  mutate(x = 23)
d2 = rbind(d2.1, d2.2) %>% 
  mutate(site = case_when(region == "danum_valley" ~ "Tropical Asia",
                          region == "pacaya_samiria" ~ "Tropical S America",
                          region == "canada_boreal" ~ "Boreal N America",
                          region == "russia_boreal" ~ "Boreal Asia",
                          region == "mark_twain_nf" ~ "Temperate N America",
                          region == "europe_temperate" ~ "Temperate Europe") )


tmaxmonth.tair.strata = ggplot() +
  geom_ribbon(data = d2 %>% filter(frelhgt == "(0,0.1]" | frelhgt == "(0.9,1]"),
              aes(x, ymin = tair.min, ymax = tair.max, fill = frelhgt), alpha = 0.1,
              show.legend = F) +
  geom_line(data = d, aes(hour, tair, color = frelhgt), linewidth = 0.25) +
  scale_x_continuous("Time", limits = c(0,23), breaks = seq(0,23,6), expand = c(0,0)) +
  scale_y_continuous("Hourly temperature (\u00b0C)", expand = c(0,0)) +
  scale_color_discrete_sequential(palette = "viridis", 
                                  guide = guide_colorsteps(title = "Relative\nheight", show.limits = T)) +
  scale_fill_manual(values = sequential_hcl(n = 10, palette = "viridis")[c(10,1)]) +
  facet_wrap(~site, nrow = 2, scales = "free_y", dir = "v") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.key.size = unit(5, "mm")) +
  thm

ggplot() +
  geom_ribbon(data = d2 %>% filter(frelhgt == "(0,0.1]" | frelhgt == "(0.9,1]"),
              aes(x, ymin = relhum.min, ymax = relhum.max, fill = frelhgt), alpha = 0.1) +
  geom_line(data = d, aes(hour, relhum, color = frelhgt)) +
  scale_x_continuous("Time", limits = c(0,24), breaks = seq(0,23,6)) +
  facet_wrap(~region) +
  theme_classic()



# MAP ---------------------------------------------------------------------

wd = vect("data/original/rnaturalearth_world.shp") %>% 
  crop(ext(-180,180,-60,90)) %>% 
  project("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m")

gedipts = pad2 %>% 
  dplyr::select(lon_lm_a0, lat_lm_a0, region) %>%
  mutate(region2 = case_when(grepl("boreal", region)~"Boreal",
                             grepl("danum|pac", region) ~ "Tropical",
                             grepl("eur|mark", region) ~ "Temperate"),
         cols = case_when(region2=="Boreal" ~ "#56B4E9",
                          region2=="Temperate" ~ "#E08A00",
                          region2=="Tropical" ~ "#008234")) %>% 
  distinct() %>% 
  vect(geom = c("lon_lm_a0", "lat_lm_a0"), crs = "epsg:4326") %>% 
  project("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m")

map = ggplot() +
  geom_spatvector(data = wd) +
  geom_spatvector(data = gedipts, aes(color = cols)) +
  scale_color_identity() +
  theme_void()


# COMBINE PLOTS -----------------------------------------------------------

(map|tmaxmonth.tair.strata) / profileplot

ggsave("figures/ms_figures/fig4_profiles.png", width = 180, height = 110, dpi = 300, unit = "mm")
