# Goal: Make response plots to visualize predictions of SES verticality across predictor variables
# These are conditional plots with all variables held at the mean, while the variable of interest is allowed to vary
# from min to max value

library(sdmTMB)
library(ggeffects)
library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)

amph = readRDS("results/sdmTMB_models2.3/amphibians_sesvert_tempdiu.rds")
amphmod = amph$fullmod$modlist[[1]]

birds = readRDS("results/sdmTMB_models2.3/birds_sesvert_tempdiu.rds")
birdsmod = birds$fullmod

mammals = readRDS("results/sdmTMB_models2.3/mammals_sesvert_tempdiu.rds")
mammalsmod = mammals$fullmod

repts = readRDS("results/sdmTMB_models2.3/reptiles_sesvert_tempdiu.rds")
reptsmod = repts$fullmod

range(birdsmod$data$log_precip_dry)
range(birdsmod$data$tmax_warm)
range(birdsmod$data$tmin_cold)
range(birdsmod$data$canopy_height)
range(birdsmod$data$veg_den)
range(birdsmod$data$precip_warm)
range(birdsmod$data$temp_diu)

var = c(
  "log_precip_dry",
  "tmax_warm",
  "tmin_cold",
  "canopy_height2",
  "veg_den",
  "precip_warm",
  "temp_diu",
  "log_precip_dry,canopy_height2",
  "tmax_warm,canopy_height2",
  "tmin_cold,canopy_height2",
  "precip_warm,canopy_height2",
  "temp_diu,canopy_height2",
  "canopy_height2,tmin_cold"
  )


mods = list(amphmod, birdsmod, mammalsmod, reptsmod) # list of model fits
names(mods) = c("Amphibians", "Birds", "Mammals", "Reptiles")
predfit = list() # store predictions for each taxon

# for each taxon,
# predict varying one variable and holding others constant at their mean
for(m in 1:length(mods)) {
  modi = mods[[m]] # select model
  pred.df = data.frame() # dataframe to store results
  for(i in var) {
    
    # if varying 2 variables (canopy height and a climate var)
    if(grepl(",", i)) {
      v = str_split_1(i, ",")
      v1 = c(-2,0,2)
      v2 = seq(min(modi$data[,v[2]]), max(modi$data[,v[2]]), 0.1)
      var_i = expand.grid(v1,v2)
      colnames(var_i)=c(v[1],v[2])
      nd = matrix(ncol=7, nrow = nrow(var_i))
      colnames(nd) = var[1:7]
      nd[,v[1]] = var_i[,v[1]]
      nd[,v[2]] = var_i[,v[2]]
      nd[is.na(nd)] <- 0
      nd = as.data.frame(nd)
      
      # predict while holding all variables except one constant at their mean
      # mean is zero because all variables were standardized
      pred = predict(modi, nd, re_form = NA, se_fit = T)
      pred = pred %>% 
        # calculate 95% CI
        # checked this calculation against ggeffects - same answer, but ggeffects takes much longer
        mutate(conf.high = est-qnorm(0.975)*est_se,
               conf.low = est+qnorm(0.975)*est_se,
               var = v[2],
               group_var = v[1]) %>% 
        dplyr::select(any_of(c(v[1], v[2], "est", "conf.low", "conf.high", "pred", "var", "group_var"))) %>% 
        rename(x = v[2], group = v[1])
    
    # if varying one variable  
    } else {
      var_i = seq(min(modi$data[,i]), max(modi$data[,i]), 0.1) # vary variable from min to max at interval of 0.1
      nd = matrix(ncol = 7, nrow = length(var_i)) # make prediction matrix
      colnames(nd) = var[1:7]
      nd[,i] = var_i
      nd[is.na(nd)] <- 0
      nd = as.data.frame(nd)
      
      # predict while holding all variables except one constant at their mean
      # mean is zero because all variables were standardized
      pred = predict(modi, nd, re_form = NA, se_fit = T)
      pred = pred %>% 
        # calculate 95% CI
        # checked this calculation against ggeffects - same answer, but ggeffects takes much longer
        mutate(conf.high = est-qnorm(0.975)*est_se,
               conf.low = est+qnorm(0.975)*est_se,
               var = i,
               group_var = NA,
               group = NA) %>% 
        dplyr::select(any_of(c(i, "est", "conf.low", "conf.high", "pred", "var"))) %>% 
        rename(x = i)
    }
    pred.df = bind_rows(pred.df, pred)
    
  }
  pred.df$taxa = names(mods)[m]
  predfit[[names(mods)[m]]] = pred.df
}

saveRDS(predfit, file = "results/sdmTMB_models2.3/conditional_effects/conditional_effects.rds")

test = predfit[["Mammals"]]
test %>% filter(is.na(group_var)) %>% 
  ggplot() +
  geom_line(aes(x, est)) +
  geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  facet_wrap(~var) +
  theme_bw()

test %>% filter(group_var == "canopy_height2") %>% 
  ggplot() +
  geom_line(aes(x, est, color = factor(group))) +
  geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high, fill = factor(group)), alpha = 0.2) +
  theme_bw()

