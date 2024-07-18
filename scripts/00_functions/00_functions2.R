mat = matrix(c(0,1,1,0,0,
              0,1,0,1,0,
              1,1,1,1,1,
              0,0,0,1,1), byrow = T, nrow = 4)

vert = c(1,0.5,0,1,1)
p.arb = c(1,0,0,1,1)
names(vert) = names(mat)

colnames(mat) = paste0("sp", 1:ncol(mat))
rownames(mat) = paste0("site", 1:nrow(mat))

mat = as.data.table(mat)



trait_summ = function(m, t, fun) {
  # summarise trait value per site
  # replace zeros with NAs
  # fun: function to summarize - if fun = "prop" calculates proportions of ones for a binary variable
  
  if(sum(!colnames(m) %in% names(t)) > 0) {
    stop("species differ between community matrix and trait vector")
  }
  if(sum(!colnames(t) %in% names(m)) > 0) {
    stop("species differ between community matrix and trait vector")
  }
  
  if(sum(colnames(m) != names(t)) > 0){
    t[order(match(names(t), colnames(m)))]
  }
  
  m[m == 0] = NA
  
  # replace vert scores in matrix and calculate summary trait value for each site
  if(fun == "prop") {
    apply(t(m)*t, 2, function(x) sum(x, na.rm = T)/sum(!is.na(x)))
  } else {
    apply(t(m)*t, 2, fun, na.rm = T)
  }
}

sestrait = function(m, t, fun, reps, quantiles = c(0.025, 0.975)) {
  # m = site x species matrix
  # t = vector of traits (names of trait vector must correspond to column names of matrix)
  # fun = function to summarize trait values per site (e.g. mean will provide mean and SES of mean trait)
  # number of replicates to generate random communities
  # quantiles = vector of quantiles to summarize replicates
  
  if(sum(!colnames(m) %in% names(t)) > 0) {
    stop("species differ between community matrix and trait vector")
  }
  if(sum(!colnames(t) %in% names(m)) > 0) {
    stop("species differ between community matrix and trait vector")
  }
  
  if(sum(colnames(m) != names(t)) > 0){
    t[order(match(names(t), colnames(m)))]
  }
  
  t_obs = trait_summ(m, t, fun) #NA values due to richness = 0
  t_rand = replicate(reps, trait_summ(randomizeMatrix(m, null.model = "richness"), fun = fun, t))
  
  t_rand_mean = apply(X = t_rand, MARGIN = 1, FUN = mean, na.rm = T)
  t_rand_sd = apply(X = t_rand, MARGIN = 1, FUN = sd, na.rm = T)
  t_ses = (t_obs - t_rand_mean)/t_rand_sd
  # calculate quantiles
  t_quant = t(apply(X = t_rand, MARGIN = 1, FUN = quantile, probs = quantiles, na.rm = T))
  
  return(cbind(t_obs, t_rand_mean, t_rand_sd, t_ses, t_quant))

}




