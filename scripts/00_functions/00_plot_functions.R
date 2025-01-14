
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
  
  
  # predtest = predict(mod)
  pred = predict(mod, type = "response")
  sum(dat$p.arb == 0) / length(dat$p.arb)
  sum(s==0) / length(s)
  s = simulate(mod, nsim = 500, seed = 12345, type = "mle-mvn")
  #s = sdmTMB:::simulate.sdmTMB(mod, nsim = 500, seed = 12345, type = "mle-mvn")
  
  dharma_resid = dharma_residuals(s, mod, return_DHARMa = T)
  plot(dharma_resid)
  testResiduals(dharma_resid)
  
  simulate(mod, nsim = 500, type = "mle-mvn") |>
    dharma_residuals(mod)
  
  richmat = matrix(rep(c(pred$rich), 500), ncol = 500)
  stest = s/richmat
  
  # simulate quantile residuals
  r <- DHARMa::createDHARMa(
    simulatedResponse = s,
    observedResponse = pred[,response_var],
    fittedPredictedResponse = pred[,"est"],
    integerResponse = integer_response
  )
  

  png(paste0(fpath, "01.png"), height = 105, width = 210, res = 300, units = "mm")
  plot(r)
  dev.off()
  
  png(paste0(fpath, "02.png"), height = 100, width = 150, res = 300, units = "mm")
  testQuantiles(r)
  dev.off()
  
  png(paste0(fpath, "03.png"), height = 100, width = 300, res = 300, units = "mm")
  testResiduals(r)
  dev.off()
  
  png(paste0(fpath, "04.png"), height = 150, width = 150, res = 300, units = "mm")
  plotResiduals(r, rank = F, form = pred[,response_var])
  dev.off()
  
  png(paste0(fpath, "05.png"), height = 150, width = 150, res = 300, units = "mm")
  plot(pred[,response_var], pred$est_non_rf, xlab = response_var, ylab = paste0(response_var, " prediction without spatial random field"))
  abline(a=0, b=1, col = "red")
  dev.off()
  
  png(paste0(fpath, "06.png"), height = 150, width = 150, res = 300, units = "mm")
  plot(pred[,response_var], pred$est, xlab = response_var, ylab = paste0(response_var, " prediction"))
  abline(a=0, b=1, col = "red")
  dev.off()
  
  # need to subset predictions to calculate spatial autocorrelation - otherwise there are too many points
  sub = sample(1:nrow(s), 10000)
  r.sub = DHARMa::createDHARMa(
    simulatedResponse = s[sub,],
    observedResponse = pred[,response_var][sub],
    fittedPredictedResponse = pred[,"est"][sub]
  )
  
  png(paste0(fpath,  "07.png"), height = 150, width = 300, res = 300, units = "mm")
  testSpatialAutocorrelation(r.sub, x = pred$x[sub], y = pred$y[sub])
  dev.off()
}

