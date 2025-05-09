
# PLOT COEFFICIENTS FOR MODEL COMPARISONS ---------------------------------

plot_compMods_coefs = function(mods, fname) {
  # mods: list of named models to plot coefs
  # fname: filename to save plot to
  coefs = foreach(i = 1:length(mods), .combine = "rbind") %do% {
    df = tidy(mods[[i]], conf.int = T, conf.level = 0.95) %>% 
      mutate(model = names(mods)[i],
             term = factor(term, levels = c("log_clim_velocity",
                                            "precip_wet", "log_precip_dry", 
                                            "I(tmax_warm^2)", "tmax_warm","tmin_cold",
                                            "veg_den", "canopy_height")))
  }
  
  p = coefs %>% 
    filter(term != "(Intercept)") %>% 
    ggplot(aes(x = estimate, y = term, xmin = conf.low, xmax = conf.high, color = model)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_pointrange(position = position_dodge2(width = 0.3)) +
    scale_x_continuous("Estimate") +
    theme_bw() + 
    theme(axis.title.y = element_blank())
  
  png(fname, height = 180, width = 180, units = "mm", res = 300) 
  p
  dev.off()
}



# CHECK RESIDUAL PLOTS ----------------------------------------------------

plot_resids = function(mod, response_var, fpath, integer_response) {
  # mod: sdmTMB model to check residuals on
  # response variable used in the model
  # file path without file extension
  

  pred = fitted(mod)
  y = mod$response[,1]
  s = simulate(mod, nsim = 500, seed = 12345, type = "mle-mvn")
  
  # simulate quantile residuals
  r <- DHARMa::createDHARMa(
    simulatedResponse = s,
    observedResponse = y,
    fittedPredictedResponse = pred,
    integerResponse = integer_response
  )
  
  png(paste0(fpath, "01.png"), height = 105, width = 210, res = 300, units = "mm")
  plot(r)
  dev.off()
  

  # need to subset residuals to test spatial autocor - otherwise too computationally intensive
  sub = sample(1:nrow(s), 10000)
  r.sub = DHARMa::createDHARMa(
    simulatedResponse = s[sub,],
    observedResponse = y[sub],
    fittedPredictedResponse = pred[sub]
  )
  
  px = mod$data$x[sub]
  py = mod$data$y[sub]
  png(paste0(fpath,  "07.png"), height = 150, width = 300, res = 300, units = "mm")
  testSpatialAutocorrelation(r.sub, x = px, y = py)
  dev.off()

  
  png(paste0(fpath, "02.png"), height = 100, width = 150, res = 300, units = "mm")
  testQuantiles(r)
  dev.off()
  
  png(paste0(fpath, "03.png"), height = 100, width = 125, res = 300, units = "mm")
  hist(r$scaledResiduals)
  dev.off()
  
  png(paste0(fpath, "04.png"), height = 150, width = 150, res = 300, units = "mm")
  plotResiduals(r, rank = T, form = mod$data[,response_var])
  dev.off()
  
  png(paste0(fpath, "06.png"), height = 150, width = 150, res = 300, units = "mm")
  plot(mod$data[,response_var], pred, xlab = response_var, ylab = paste0(response_var, " prediction"))
  abline(a=0, b=1, col = "red")
  dev.off()
  
}

