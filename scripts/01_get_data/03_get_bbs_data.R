# download breeding bird survey data for select regions and examine changes in arboreality over time

library(bbsAssistant)
library(data.table)

sb_items
View(region_codes)

grab_bbs_data(bbs_dir = "data/original/bbs_data")

dat = fread("data/original/bbs_data/50-StopData/50-StopData/1997ToPresent_SurveyWide/Fifty3/fifty3.csv")
head(dat)

dat = fread("data/original/bbs_data/50-StopData/50-StopData/1997ToPresent_SurveyWide/Fifty6/fifty6.csv")


# FL state num
fl = region_codes %>% filter(State == "florida")
nc = region_codes %>% filter(State == "north carolina")

dat.fl = dat %>% filter(StateNum == nc$StateNum)

# calculate presab on route
dat.fl[, pa := apply(.SD, 1, sum), .SDcols = Stop1:Stop50]
dat.fl[, pa := ifelse(pa>0, 1, 0)]

# reomve stop columns
dat.fl[, 8:57 := NULL]

# add species names
head(species_list)
splist = species_list %>% dplyr::select(AOU, Scientific_Name)

dat.fl = left_join(dat.fl, splist, by = "AOU")
dat.fl = dat.fl %>% rename(sciname = Scientific_Name)

# verticality data for birds
vert =  read.csv("data/derivative_data/species_data/birds_breedresident_spdat_eltonvert.csv")
vert = vert %>% 
  dplyr::select(sciname, Verticality) %>% 
  mutate(sciname = gsub("_", " ", sciname))

# combine vert data with bbs data for FL
dat.fl = left_join(dat.fl, vert, by = "sciname")

dat.fl.sum = dat.fl %>% 
  filter(pa == 1) %>% 
  group_by(RouteDataID, CountryNum, StateNum, Route, Year) %>% 
  summarise(vert = mean(Verticality, na.rm = T))

ggplot(dat.fl.sum, aes(x = Year, y = vert)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic()
















