library(foreach)
library(tidyterra)
library(ggplot2)
library(terra)
library(colorspace)
library(patchwork)
library(fmesher)
library(sdmTMB)

# testing data
# occ = matrix(c(0,1,1,0,0,0,
#              0,1,0,1,0,0,
#              1,1,1,1,1,0,
#              0,0,0,1,1,0,
#              1,0,1,1,0,0,
#              0,1,1,1,1,0,
#              0,0,0,1,1,1,
#              0,0,0,1,0,1,
#              0,1,0,0,0,1), byrow = T, nrow = 9)
# 
# tr = data.frame(Verticality = c(1,0.5,0,1,1,0.5),
#            Arb = c(1,0,0,1,1,0),
#            Fos = c(0,0,1,0,0,0),
#            Ter = c(0,1,0,0,0,1),
#            sciname = paste0("sp", 1:ncol(occ)))
# colnames(occ) = paste0("sp", 1:ncol(occ))
# rownames(occ) = paste0("site", 1:nrow(occ))
# 
# env = data.frame(realm = c(rep("a",4), NA, rep("b", 4)), biome = c(rep("b", 4), NA, rep("b2", 4)))

gridcell_dat = function(occ, tr, env, reps = 100) {
  # test if species names in trait_dat are in same order as occ
  if(sum(colnames(occ)[3:ncol(occ)] != tr$sciname) > 0) {
    stop("Species names in occurrence dataframe and trait dataframe do not match")
  }
  
  # test match between occ and env
  if(all.equal(env[,c("x", "y")], occ[,c("x", "y")]) != TRUE){
    stop("env and occ coordinates do not match")
  }
  
  # name occ rows and remove xy columns
  rownames(occ) = paste0(occ$x, "_", occ$y)
  rownames(env) = paste0(env$x, "_", env$y)
  occ = occ[,-c(1,2)]
  
  realms = unique(env$realm)
  realms = realms[which(!is.na(realms))]
  
  cl = makeCluster(7)
  registerDoParallel(cl)
  
  dat = foreach(i = realms,
          .packages = "picante") %do% {
    # filter to realm
    occ.r = occ[which(env$realm==i),]
    colnames(occ.r) = colnames(occ)
    #occ.r = as.matrix(occ.r)
    # remove species not in the realm
    #test = apply(occ.r, 2, sum, na.rm = T)
    #test = which(test > 0)
    occ.r = occ.r[,which(apply(occ.r, 2, sum, na.rm = T)>0)]
    tr.r = tr[which(tr$sciname %in% colnames(occ.r)),]
    env.r = env[which(env$realm==i),]
    
    # calculate richness
    rich = apply(occ.r, 1, sum)
    
    arb = tr.r$Arb
    names(arb) = tr.r$sciname
    p.arb = trait_summ(occ.r, arb, fun = "prop")
    
    fos = tr.r$Fos
    names(fos) = tr.r$sciname
    p.fos = trait_summ(occ.r, fos, fun = "prop")
    
    ter = tr.r$Ter
    names(ter) = tr.r$sciname
    p.ter = trait_summ(occ.r, ter, fun = "prop")
    
    vert = tr.r$Verticality
    names(vert) = tr.r$sciname
    sesvert = sestrait(occ.r, vert, fun = "mean", reps = reps, quantiles = c(0.025, 0.975))
    colnames(sesvert) = c("vert.mean", "vert.rand.mean", "vert.rand.sd", 
                            "vert.mean.ses", "vert.rand.mean.025", "vert.rand.mean.975")
    
    celldat.r = cbind(rich, p.arb, p.fos, p.ter, sesvert, env.r)
    occ.r$x = celldat.r$x
    occ.r$y = celldat.r$y
    
    list(celldat = celldat.r, occ = occ.r)
    
  }
  stopCluster(cl)
  
  return(dat)
}



# FIT MESH ----------------------------------------------------------------

fit_mesh = function(f, dat, range, v, family, wts = NULL) {
  # Determine optimal mesh size and sensitivity of results to changes in the mesh
  # optimal mesh size occurs when the range of the resulting model is at least 5 times the max.edge length
  
  # f = formula to fit the model
  # dat = data frame with response and predictor variables and x,y coords (colnames must be named "x" and "y")
  # range = starting range for the mesh (max.edge is 1/5 of the starting range)
  # vector of world outline
  
    # initialize max.edge ratio
  edge.ratio = 0.2
  
  # initialize iteration
  iter = 1
  
  # initialize list to store models and meshes
  mesh = list()
  mod = list()
  
  while(edge.ratio >= 0.2) {
    print(iter)
    max.edge = range/5
      
    # estimate mesh
    m = fm_mesh_2d(loc = dat[,c('x', 'y')], cutoff = max.edge/5, 
                   max.edge = c(max.edge, max.edge*5),  
                   offset = c(max.edge, max.edge*10),
                   crs = fm_crs("+proj=cea +datum=WGS84"))
    mesh[[iter]] = make_mesh(dat[,c('x', 'y')], xy_cols = c("x", "y"), mesh = m)

    iter2 = 1
    all_ok = FALSE
    while(!all_ok & iter2 <= 2) {
      print(paste0("iter2 = ", iter2))
      if(iter2 == 1) {
            mod[[iter]] = sdmTMB(f,
                 data = dat,
                 weights = wts,
                 mesh = mesh[[iter]],
                 spatial = "on",
                 reml = T, 
                 family = family,
                 spatial_varying = ~ 0 + canopy_height,
                 control = sdmTMBcontrol(eval.max = 8000, iter.max = 4000))

      } else {
        mod[[iter]] = sdmTMB(f,
                             data = dat,
                             weights = wts,
                             mesh = mesh[[iter]],
                             spatial = "on",
                             reml = T, 
                             family = family,
                             spatial_varying = ~ 0 + canopy_height,
                             control = sdmTMBcontrol(eval.max = 8000, iter.max = 4000, start = list(ln_phi = ln_phi,
                                                                                                    ln_tau_O = ln_tau_o,
                                                                                                    ln_tau_Z = ln_tau_z,
                                                                                                    ln_kappa = ln_kappa)))
        
      }
      
      all_ok = sanity(mod[[iter]])$all_ok
      pars = get_pars(mod[[iter]])
      ln_phi = pars$ln_phi
      ln_tau_o = pars$ln_tau_O
      ln_tau_z = pars$ln_tau_Z
      ln_kappa = pars$ln_kappa
      iter2 = iter2+1
      gc()
    }
    
    
    # update mesh with range parameter from the previous model
    ranpars = tidy(mod[[iter]], "ran_pars", conf.int = TRUE)
    range = ranpars$estimate[ranpars$term == "range"]
    
    edge.ratio = max.edge/range
    iter = iter + 1
    gc()

  }
  
  # compare model coefs
  coefs = foreach(i = 1:length(mod), .combine = "rbind") %do% {
    coefs = tidy(mod[[i]], conf.int = T, conf.level = 0.95) %>% 
      mutate(model = paste0("mesh",i))
  }
  
  ranpars = foreach(i = 1:length(mod), .combine = "rbind") %do% {
    ranpars = tidy(mod[[i]], "ran_pars", conf.int = T, conf.level = 0.95) %>% 
      mutate(model = paste0("mesh",i))
  }
  
  coefs.plt = ggplot(coefs, aes(x = estimate, y = term, xmin = conf.low, xmax = conf.high, color = model)) +
    geom_pointrange(position = position_dodge2(width = 0.5)) +
    theme_bw()
  
  
  ranpars.plt = ggplot(ranpars, aes(x = estimate, y = term, xmin = conf.low, xmax = conf.high, color = model)) +
    geom_pointrange(position = position_dodge2(width = 0.5)) +
    theme_bw()
  
  plot_spatial_varying = function(mod, var, v, legend_title) {
    coefs = tidy(mod)
    coef.var = as.numeric(coefs[which(coefs$term == var), "estimate"])
    pred = predict(mod, type = "response")
    datpred = pred %>% 
      mutate(x = x*1e5, y = y*1e5,
             canopy_height_est_slope = zeta_s_canopy_height + coef.var) %>% 
      relocate(x, y, .before = rich) %>% 
      rast(crs = "+proj=cea + datum=WGS84") %>% 
      project("epsg:4326")
    
    
    ggplot() +
      #geom_spatvector(data = v, color = "black", fill = "gray70") +
      geom_spatraster(data = datpred, aes(fill = canopy_height_est_slope)) +
      #geom_spatvector(data = v, color = "black", fill = NA) +
      scale_fill_continuous_divergingx("BrBG", na.value = NA, guide = guide_colorbar(legend_title)) +
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_rect(color = "black", fill = NA),
            panel.grid = element_blank(),
            legend.key.height = unit(3, "mm"))
    
  }
  
  svc.plts = foreach(i = 1:length(mod)) %do% {
    plot_spatial_varying(mod = mod[[i]], var = "canopy_height", v = v, legend_title = "coefficient") +
      ggtitle(paste0("model_", i))
  }
  
  svcplt = wrap_plots(svc.plts, ncol = 2)
  
  # the last model and mesh in the list will be the optimal one
  return(list(meshes = mesh, mods = mod, max.edge = max.edge,
              coefs.plt = coefs.plt, ranpars.plt = ranpars.plt, svc.plot = svcplt))
  
}

# FIT MESH RICHNESS ----------------------------------------------------------------

fit_mesh_richness = function(f, dat, range, v, family, wts = NULL, edge.ratio) {
  # Determine optimal mesh size and sensitivity of results to changes in the mesh
  # optimal mesh size occurs when the range of the resulting model is at least 5 times the max.edge length
  
  # f = formula to fit the model
  # dat = data frame with response and predictor variables and x,y coords (colnames must be named "x" and "y")
  # range = starting range for the mesh (max.edge is 1/5 of the starting range)
  # vector of world outline
  
  # initialize max.edge ratio
  er = edge.ratio # the ratio of max edge to spatial range. Generally max edge should be less than 1/5 of the spatial range (so edge ratio should be 0.2), but sometimes this does not converge
  
  # initialize iteration
  iter = 1
  
  # initialize list to store models and meshes
  mesh = list()
  mod = list()
  
  while(er >= edge.ratio) {
    print(iter)
    max.edge = range/5
    
    # estimate mesh
    m = fm_mesh_2d(loc = dat[,c('x', 'y')], cutoff = max.edge/5, 
                   max.edge = c(max.edge, max.edge*5),  
                   offset = c(max.edge, max.edge*10),
                   crs = fm_crs("+proj=cea +datum=WGS84"))
    mesh[[iter]] = make_mesh(dat[,c('x', 'y')], xy_cols = c("x", "y"), mesh = m)
    
    iter2 = 1
    all_ok = FALSE
    while(!all_ok & iter2 <= 2) {
      print(paste0("iter2 = ", iter2))
      if(iter2 == 1) {
        mod[[iter]] = sdmTMB(f,
                             data = dat,
                             weights = wts,
                             mesh = mesh[[iter]],
                             spatial = "on",
                             reml = T, 
                             family = family,
                             control = sdmTMBcontrol(eval.max = 8000, iter.max = 4000))
        
      } else {
        mod[[iter]] = sdmTMB(f,
                             data = dat,
                             weights = wts,
                             mesh = mesh[[iter]],
                             spatial = "on",
                             reml = T, 
                             family = family,
                             control = sdmTMBcontrol(eval.max = 8000, iter.max = 4000, start = list(ln_phi = ln_phi,
                                                                                                    ln_tau_O = ln_tau_o,
                                                                                                    ln_tau_Z = ln_tau_z,
                                                                                                    ln_kappa = ln_kappa)))
        
      }
      
      all_ok = sanity(mod[[iter]])$all_ok
      pars = get_pars(mod[[iter]])
      ln_phi = pars$ln_phi
      ln_tau_o = pars$ln_tau_O
      ln_tau_z = pars$ln_tau_Z
      ln_kappa = pars$ln_kappa
      iter2 = iter2+1
      gc()
    }
    
    
    # update mesh with range parameter from the previous model
    ranpars = tidy(mod[[iter]], "ran_pars", conf.int = TRUE)
    range = ranpars$estimate[ranpars$term == "range"]
    
    er = max.edge/range
    iter = iter + 1
    gc()
    
  }
  
  # compare model coefs
  coefs = foreach(i = 1:length(mod), .combine = "rbind") %do% {
    coefs = tidy(mod[[i]], conf.int = T, conf.level = 0.95) %>% 
      mutate(model = paste0("mesh",i))
  }
  
  ranpars = foreach(i = 1:length(mod), .combine = "rbind") %do% {
    ranpars = tidy(mod[[i]], "ran_pars", conf.int = T, conf.level = 0.95) %>% 
      mutate(model = paste0("mesh",i))
  }
  
  coefs.plt = ggplot(coefs, aes(x = estimate, y = term, xmin = conf.low, xmax = conf.high, color = model)) +
    geom_pointrange(position = position_dodge2(width = 0.5)) +
    theme_bw()
  
  
  ranpars.plt = ggplot(ranpars, aes(x = estimate, y = term, xmin = conf.low, xmax = conf.high, color = model)) +
    geom_pointrange(position = position_dodge2(width = 0.5)) +
    theme_bw()
  
  
  
  # the last model and mesh in the list will be the optimal one
  return(list(meshes = mesh, mods = mod, max.edge = max.edge,
              coefs.plt = coefs.plt, ranpars.plt = ranpars.plt))
  
}



# COMPARE MODELS USING AIC ------------------------------------------------

compareMods_AIC = function(f, dat, mesh, taxon, response_var, family, wts = NULL) {
  # f: formula for sdmTMB model
  # dat = dataframe for sdmTMB model
  # mesh = mesh used for sdmTMB model
  # taxon: character, taxon being modeled - used for labelling output
  # response_var: character; response variable for the model - used for labelling output
  all_ok = FALSE
  iter = 1
  while(!all_ok & iter <= 2) {
    print(iter)
    if(iter == 1) {
      mod = sdmTMB(f, 
                   data = dat,
                   weights = wts,
                   mesh = mesh,
                   spatial = "on",
                   reml = T, 
                   family = family,
                   control = sdmTMBcontrol(eval.max = 6000, iter.max = 3000))
    } else {
      mod = sdmTMB(f, 
                   data = dat,
                   weights = wts,
                   mesh = mesh,
                   spatial = "on",
                   reml = T,
                   family = family,
                   control = sdmTMBcontrol(eval.max = 6000, iter.max = 3000, 
                                           start = list(ln_phi = ln_phi,
                                                        ln_tau_O = ln_tau_o)))
    }
    all_ok = sanity(mod)$all_ok
    pars = get_pars(mod)
    ln_phi = pars$ln_phi
    ln_tau_o = pars$ln_tau_O
    iter = iter+1
    gc()
  }
  
  
  # model with svc but no random effect
  all_ok = FALSE
  iter = 1
  while(!all_ok & iter <= 2) {
    print(iter)
    if(iter == 1) {
      mod.svc = sdmTMB(f, 
                       data = dat,
                       weights = wts,
                       mesh = mesh,
                       spatial = "on",
                       reml = T, 
                       spatial_varying = ~ 0 + canopy_height,
                       family = family,
                       control = sdmTMBcontrol(eval.max = 6000, iter.max = 3000))
    } else {
      mod.svc = sdmTMB(f, 
                       data = dat,
                       weights = wts,
                       mesh = mesh,
                       spatial = "on",
                       reml = T, 
                       spatial_varying = ~ 0 + canopy_height,
                       family = family,
                       control = sdmTMBcontrol(eval.max = 6000, iter.max = 3000,
                                               start = list(ln_phi = ln_phi)))
    }
    all_ok = sanity(mod.svc)$all_ok
    ln_phi = get_pars(mod.svc)$ln_phi
    iter = iter+1
    gc()
  }
  
  
  # model with realm random intercept but no svc
  all_ok = FALSE
  iter = 1
  while(!all_ok & iter <= 2) {
    print(iter)
    if(iter == 1) {
      mod.realm = sdmTMB(update(f, ~ . + (1|realm)), 
                         data = dat,
                         weights = wts,
                         mesh = mesh,
                         spatial = "on",
                         reml = T,
                         family = family,
                         control = sdmTMBcontrol(eval.max = 6000, iter.max = 3000))
    } else {
      mod.realm = sdmTMB(update(f, ~ . + (1|realm)), 
                         data = dat,
                         weights = wts,
                         mesh = mesh,
                         spatial = "on",
                         reml = T,
                         family = family,
                         control = sdmTMBcontrol(eval.max = 6000, iter.max = 3000,
                                                 start = list(ln_phi = ln_phi,
                                                              ln_kappa = ln_kappa,
                                                              ln_tau_O = ln_tau_o)))
    }
    all_ok = sanity(mod.realm)$all_ok
    pars = get_pars(mod.realm)
    ln_phi = pars$ln_phi
    ln_kappa = pars$ln_kappa
    ln_tau_o = pars$ln_tau_O
    iter = iter+1
    gc()
  }
  
  
  
  # model with realm random intercept and svc
  all_ok = FALSE
  iter = 1
  while(!all_ok & iter <= 2) {
    print(iter)
    if(iter == 1) {
      # run model first time
      mod.realm.svc = sdmTMB(update(f, ~ . + (1|realm)), 
                             data = dat,
                             weights = wts,
                             mesh = mesh,
                             spatial = "on",
                             reml = T, 
                             spatial_varying = ~ 0 + canopy_height,
                             family = family,
                             control = sdmTMBcontrol(eval.max = 6000, iter.max = 3000))
    } else {
      # if model is not ok and has to run a second time
      # set ln_phi to the estimation of ln_phi from the previous model
      mod.realm.svc = sdmTMB(update(f, ~ . + (1|realm)), 
                             data = dat,
                             weights = wts,
                             mesh = mesh,
                             spatial = "on",
                             reml = T, 
                             spatial_varying = ~ 0 + canopy_height,
                             family = family,
                             control = sdmTMBcontrol(eval.max = 6000, iter.max = 3000,
                                                     start = list(ln_phi = ln_phi)))
    }
    
    all_ok = sanity(mod.realm.svc)$all_ok
    ln_phi = get_pars(mod.realm.svc)$ln_phi
    iter = iter+1
    gc()
  }

  
  # compare AIC and log-likelihood of models
  aic = AIC(mod.realm.svc, mod.realm, mod.svc, mod)
  min_aic = min(aic$AIC)
  aic = aic %>% 
    mutate(deltaAIC = aic$AIC - min_aic) %>% 
    rownames_to_column(var = "model")
  
  ll = sort(c(mod.realm.svc = as.numeric(logLik(mod.realm.svc)),
              mod.realm = as.numeric(logLik(mod.realm)),
              mod.svc = as.numeric(logLik(mod.svc)),
              mod = as.numeric(logLik(mod)))) %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "model") %>% 
    rename(logLik = ".")
  
  comp = inner_join(ll, aic, by = "model") %>% 
    mutate(taxon = taxon, response_var = response_var) %>% 
    arrange(deltaAIC)
  
  return(list(
    mods = list(mod = mod, mod.realm = mod.realm, mod.svc = mod.svc, mod.realm.svc = mod.realm.svc),
    compare = comp
  ))
}



# COMPARE MODELS USING CROSS VALIDATION -----------------------------------

compare_cv = function(f, dat, mesh, folds, parallel = F, taxon, response_var) {
  # f: formula for sdmTMB model
  # dat = dataframe for sdmTMB model
  # mesh = mesh used for sdmTMB model
  # folds for cross validation
  # logical: whether or not to run in parallel - default is false
  # taxon: character, taxon being modeled - used for labelling output
  # response_var: character; response variable for the model - used for labelling output
  
  # model with no svc or random intercept
  all_ok = FALSE
  iter = 1
  while(sum(all_ok) < 5  & iter <= 2) {
    if(iter == 1) {
      if(parallel) {plan(multisession)}
      mod.cv = sdmTMB_cv(f, 
                         data = dat,
                         mesh = mesh,
                         spatial = "on",
                         reml = T,
                         fold_ids = folds,
                         k_folds = length(unique(folds)),
                         control = sdmTMBcontrol(eval.max = 6000, iter.max = 3000),
                         parallel = parallel,
                         use_initial_fit = T)
      if(parallel) {plan(sequential)}
    } else {
      if(parallel) {plan(multisession)}
      mod.cv = sdmTMB_cv(f, 
                         data = dat,
                         mesh = mesh,
                         spatial = "on",
                         reml = T, 
                         fold_ids = folds,
                         k_folds = length(unique(folds)),
                         control = sdmTMBcontrol(eval.max = 6000, iter.max = 3000,
                                                 start = list(ln_phi = ln_phi, 
                                                              ln_tau_o = ln_tau_o,
                                                              ln_kappa = ln_kappa)),
                         future_globals = c("ln_phi", "ln_tau_o", "ln_kappa"),
                         parallel = parallel,
                         use_initial_fit = T)
      if(parallel) {plan(sequential)}
    }
    all_ok = lapply(mod.cv$models, sanity)
    all_ok = sapply(all_ok, "[[", "all_ok")
    pars = lapply(mod.cv$models, get_pars)
    ln_phi = mean(sapply(pars, "[[", "ln_phi"))
    ln_tau_o = mean(sapply(pars, "[[", "ln_tau_O"))
    ln_kappa = sapply(pars, "[[", "ln_kappa")
    ln_kappa = matrix(apply(ln_kappa, 1, mean), nrow = 2)
    iter = iter+1
    gc()
  }
  mod.cv_ok = all_ok
  save(mod.cv, mod.cv_ok, file = "tempfiles/rept_mod1.RData")
  
  
  
  all_ok = FALSE
  iter = 1
  while(sum(all_ok) < 5  & iter <= 2) {
    print(iter)
    if(iter == 1) {
      # run model first time
      if(parallel) {plan(multisession)}
      mod.svc.cv = sdmTMB_cv(f, 
                             data = dat,
                             mesh = mesh,
                             spatial = "on",
                             reml = T, 
                             spatial_varying = ~ 0 + canopy_height,
                             fold_ids = folds,
                             k_folds = length(unique(folds)),
                             control = sdmTMBcontrol(eval.max = 6000, iter.max = 3000),
                             parallel = parallel)
      if(parallel) {plan(sequential)}
    } else {
      # if model is not ok and has to run a second time
      # set ln_phi to the estimation of ln_phi from the previous model
      #plan(multisession)
      mod.svc.cv = sdmTMB_cv(f, 
                             data = dat,
                             mesh = mesh,
                             spatial = "on",
                             reml = T, 
                             spatial_varying = ~ 0 + canopy_height,
                             fold_ids = folds,
                             k_folds = length(unique(folds)),
                             control = sdmTMBcontrol(eval.max = 6000, iter.max = 3000,
                                                     start = list(ln_phi = ln_phi, 
                                                                  ln_tau_O = ln_tau_o,
                                                                  ln_kappa = ln_kappa)),
                             future_globals = c("ln_phi", "ln_tau_o", "ln_kappa"),
                             parallel = parallel)
      if(parallel) {plan(sequential)}
    }
    
    all_ok = lapply(mod.svc.cv$models, sanity)
    all_ok = sapply(all_ok, "[[", "all_ok")
    pars = lapply(mod.svc.cv$models, get_pars)
    ln_phi = mean(sapply(pars, "[[", "ln_phi"))
    ln_tau_o = mean(sapply(pars, "[[", "ln_tau_O"))
    ln_kappa = sapply(pars, "[[", "ln_kappa")
    ln_kappa = matrix(apply(ln_kappa, 1, mean), nrow = 2)
    iter = iter+1
    gc()
  }
  mod.svc.cv_ok = all_ok
  save(mod.svc.cv, mod.svc.cv_ok, file = "tempfiles/rept_mod2.RData")
  
  # model with realm but no svc
  all_ok = FALSE
  iter = 1
  while(sum(all_ok) < 5  & iter <= 2) {
    print(iter)
    if(iter == 1) {
      if(parallel) {plan(multisession)}
      mod.realm.cv = sdmTMB_cv(update(f, ~ . + (1|realm)), 
                               data = dat,
                               mesh = mesh,
                               spatial = "on",
                               reml = T, 
                               fold_ids = folds,
                               k_folds = length(unique(folds)),
                               control = sdmTMBcontrol(eval.max = 6000, iter.max = 3000),
                               parallel = parallel)
      if(parallel) {plan(sequential)}
    } else {
      if(parallel) {plan(multisession)}
      mod.realm.cv = sdmTMB_cv(update(f, ~ . + (1|realm)), 
                               data = dat,
                               mesh = mesh,
                               spatial = "on",
                               reml = T, 
                               fold_ids = folds,
                               k_folds = length(unique(folds)),
                               control = sdmTMBcontrol(eval.max = 6000, iter.max = 3000,
                                                       start = list(ln_phi = ln_phi,
                                                                    ln_tau_O = ln_tau_o,
                                                                    ln_kappa = ln_kappa)),
                               parallel = parallel)
      if(parallel) {plan(sequential)}
    }
    all_ok = lapply(mod.realm.cv$models, sanity)
    all_ok = sapply(all_ok, "[[", "all_ok")
    pars = lapply(mod.realm.cv$models, get_pars)
    ln_phi = mean(sapply(pars, "[[", "ln_phi"))
    ln_tau_o = mean(sapply(pars, "[[", "ln_tau_O"))
    ln_kappa = sapply(pars, "[[", "ln_kappa")
    ln_kappa = matrix(apply(ln_kappa, 1, mean), nrow = 2)
    iter = iter+1
    gc()
  }
  mod.realm.cv_ok = all_ok
  save(mod.realm.cv, mod.realm.cv_ok, file = "tempfiles/rept_mod3.RData")
  
  
  # model with realm and svc
  # model with realm and svc
  all_ok = FALSE
  iter = 1
  while(sum(all_ok) < 5  & iter <= 2) {
    if(iter == 1) {
      # run model first time
      if(parallel) {plan(multisession)}
      mod.realm.svc.cv = sdmTMB_cv(update(f, ~ . + (1|realm)), 
                                   data = dat,
                                   mesh = mesh,
                                   spatial = "on",
                                   reml = T, 
                                   spatial_varying = ~ 0 + canopy_height,
                                   fold_ids = folds,
                                   k_folds = length(unique(folds)),
                                   control = sdmTMBcontrol(eval.max = 6000, iter.max = 3000),
                                   parallel = parallel)
      if(parallel) {plan(sequential)}
    } else {
      # if model is not ok and has to run a second time
      # set ln_phi to the estimation of ln_phi from the previous model
      #plan(multisession)
      mod.realm.svc.cv = sdmTMB_cv(update(f, ~ . + (1|realm)), 
                                   data = dat,
                                   mesh = mesh,
                                   spatial = "on",
                                   reml = T, 
                                   spatial_varying = ~ 0 + canopy_height,
                                   fold_ids = folds,
                                   k_folds = length(unique(folds)),
                                   control = sdmTMBcontrol(eval.max = 8000, iter.max = 4000,
                                                           start = list(ln_phi = ln_phi,
                                                                        ln_tau_O = 3,
                                                                        ln_kappa = ln_kappa)),
                                   parallel = parallel)
      
      if(parallel) {plan(sequential)}
    }
    
    all_ok = lapply(mod.realm.svc.cv$models, sanity)
    all_ok = sapply(all_ok, "[[", "all_ok")
    pars = lapply(mod.realm.svc.cv$models, get_pars)
    ln_phi = mean(sapply(pars, "[[", "ln_phi"))
    ln_tau_o = mean(sapply(pars, "[[", "ln_tau_O"))
    ln_kappa = sapply(pars, "[[", "ln_kappa")
    ln_kappa = matrix(apply(ln_kappa, 1, mean), nrow = 2)
    iter = iter+1
    gc()
  }
  mod.realm.svc.cv_ok = all_ok
  save(mod.realm.svc.cv, mod.realm.svc.cv_ok, file = "tempfiles/rept_mod4.RData")
  
  # since some of the folds did not converge, get average fold logLik across the models that did converge
  mod.cv_ll = sum(mod.cv$fold_loglik[mod.cv_ok])/sum(mod.cv_ok)
  mod.realm.cv_ll = sum(mod.realm.cv$fold_loglik[mod.realm.cv_ok])/sum(mod.realm.cv_ok)
  mod.svc.cv_ll = sum(mod.svc.cv$fold_loglik[mod.svc.cv_ok])/sum(mod.svc.cv_ok)
  mod.realm.svc.cv_ll = sum(mod.realm.svc.cv$fold_loglik[mod.realm.svc.cv_ok])/sum(mod.realm.svc.cv_ok)
  
  comp_ll = data.frame(model = c("mod", "mod.realm", "mod.svc", "mod.realm.svc"),
                       avg_logLik = c(mod.cv_ll, mod.realm.cv_ll, mod.svc.cv_ll, mod.realm.svc.cv_ll)) %>% 
    arrange(desc(avg_logLik)) %>% 
    mutate(taxon = taxon, response_var = response_var)
  
  return(list(
    mods = list(mod.cv = mod.cv, mod.realm.cv = mod.realm.cv, mod.svc.cv = mod.svc.cv, mod.realm.svc.cv = mod.realm.svc.cv),
    all_ok = c(mod.cv = mod.cv_ok, mod.realm.cv = mod.realm.cv_ok, mod.svc.cv = mod.svc.cv_ok, mod.realm.svc.cv = mod.realm.svc.cv_ok),
    compMods_cv = comp_ll
  ))
  
}

# beta distribution for mean verticality
compare_cv_beta = function(f, dat, mesh, folds, parallel = F, taxon, response_var) {
  # f: formula for sdmTMB model
  # dat = dataframe for sdmTMB model
  # mesh = mesh used for sdmTMB model
  # folds for cross validation
  # logical: whether or not to run in parallel - default is false
  # taxon: character, taxon being modeled - used for labelling output
  # response_var: character; response variable for the model - used for labelling output
  
  # model with no svc or random intercept
  all_ok = FALSE
  iter = 1
  while(sum(all_ok) < 5  & iter <= 2) {
    print(iter)
    if(iter == 1) {
      if(parallel) {plan(multisession)}
      mod.cv = sdmTMB_cv(f, 
                         data = dat,
                         mesh = mesh,
                         spatial = "on",
                         reml = T, 
                         family = Beta(),
                         fold_ids = folds,
                         k_folds = length(unique(folds)),
                         control = sdmTMBcontrol(eval.max = 6000, iter.max = 3000),
                         parallel = parallel,
                         use_initial_fit = T)
      if(parallel) {plan(sequential)}
    } else {
      if(parallel) {plan(multisession)}
      mod.cv = sdmTMB_cv(f, 
                         data = dat,
                         mesh = mesh,
                         spatial = "on",
                         reml = T, 
                         fold_ids = folds,
                         k_folds = length(unique(folds)),
                         control = sdmTMBcontrol(eval.max = 6000, iter.max = 3000,
                                                 start = list(ln_phi = ln_phi, 
                                                              ln_tau_O = ln_tau_o,
                                                              ln_kappa = ln_kappa)),
                         future_globals = c("ln_phi", "ln_tau_o", "ln_kappa"),
                         parallel = parallel,
                         family = Beta(),
                         use_initial_fit = T)
      if(parallel) {plan(sequential)}
    }
    all_ok = lapply(mod.cv$models, sanity)
    all_ok = sapply(all_ok, "[[", "all_ok")
    pars = lapply(mod.cv$models, get_pars)
    ln_phi = mean(sapply(pars, "[[", "ln_phi"))
    ln_tau_o = mean(sapply(pars, "[[", "ln_tau_O"))
    ln_kappa = sapply(pars, "[[", "ln_kappa")
    ln_kappa = matrix(apply(ln_kappa, 1, mean), nrow = 2)
    iter = iter+1
    gc()
  }
  mod.cv_ok = all_ok
  save(mod.cv, mod.cv_ok, file = "tempfiles/rept_mod1.RData")
  
  
  
  all_ok = FALSE
  iter = 1
  while(sum(all_ok) < 5  & iter <= 2) {
    print(iter)
    if(iter == 1) {
      # run model first time
      if(parallel) {plan(multisession)}
      mod.svc.cv = sdmTMB_cv(f, 
                             data = dat,
                             mesh = mesh,
                             spatial = "on",
                             reml = T, 
                             spatial_varying = ~ 0 + canopy_height,
                             fold_ids = folds,
                             k_folds = length(unique(folds)),
                             family = Beta(),
                             control = sdmTMBcontrol(eval.max = 6000, iter.max = 3000),
                             parallel = parallel)
      if(parallel) {plan(sequential)}
    } else {
      # if model is not ok and has to run a second time
      # set ln_phi to the estimation of ln_phi from the previous model
      #plan(multisession)
      mod.svc.cv = sdmTMB_cv(f, 
                             data = dat,
                             mesh = mesh,
                             spatial = "on",
                             reml = T, 
                             spatial_varying = ~ 0 + canopy_height,
                             fold_ids = folds,
                             k_folds = length(unique(folds)),
                             family = Beta(),
                             control = sdmTMBcontrol(eval.max = 6000, iter.max = 3000,
                                                     start = list(ln_phi = ln_phi, 
                                                                  ln_tau_O = ln_tau_o,
                                                                  ln_kappa = ln_kappa)),
                             future_globals = c("ln_phi", "ln_tau_o", "ln_kappa"),
                             parallel = parallel)
      if(parallel) {plan(sequential)}
    }
    
    all_ok = lapply(mod.svc.cv$models, sanity)
    all_ok = sapply(all_ok, "[[", "all_ok")
    pars = lapply(mod.svc.cv$models, get_pars)
    ln_phi = mean(sapply(pars, "[[", "ln_phi"))
    ln_tau_o = mean(sapply(pars, "[[", "ln_tau_O"))
    ln_kappa = sapply(pars, "[[", "ln_kappa")
    ln_kappa = matrix(apply(ln_kappa, 1, mean), nrow = 2)
    iter = iter+1
    gc()
  }
  mod.svc.cv_ok = all_ok
  save(mod.svc.cv, mod.svc.cv_ok, file = "tempfiles/rept_mod2.RData")
  
  # model with realm but no svc
  all_ok = FALSE
  iter = 1
  while(sum(all_ok) < 5  & iter <= 2) {
    print(iter)
    if(iter == 1) {
      if(parallel) {plan(multisession)}
      mod.realm.cv = sdmTMB_cv(update(f, ~ . + (1|realm)), 
                               data = dat,
                               mesh = mesh,
                               spatial = "on",
                               reml = T, 
                               fold_ids = folds,
                               k_folds = length(unique(folds)),
                               family = Beta(),
                               control = sdmTMBcontrol(eval.max = 6000, iter.max = 3000),
                               parallel = parallel)
      if(parallel) {plan(sequential)}
    } else {
      if(parallel) {plan(multisession)}
      mod.realm.cv = sdmTMB_cv(update(f, ~ . + (1|realm)), 
                               data = dat,
                               mesh = mesh,
                               spatial = "on",
                               reml = T, 
                               fold_ids = folds,
                               k_folds = length(unique(folds)),
                               family = Beta(),
                               control = sdmTMBcontrol(eval.max = 6000, iter.max = 3000,
                                                       start = list(ln_phi = ln_phi,
                                                                    ln_tau_O = ln_tau_o,
                                                                    ln_kappa = ln_kappa)),
                               parallel = parallel)
      if(parallel) {plan(sequential)}
    }
    all_ok = lapply(mod.realm.cv$models, sanity)
    all_ok = sapply(all_ok, "[[", "all_ok")
    pars = lapply(mod.realm.cv$models, get_pars)
    ln_phi = mean(sapply(pars, "[[", "ln_phi"))
    ln_tau_o = mean(sapply(pars, "[[", "ln_tau_O"))
    ln_kappa = sapply(pars, "[[", "ln_kappa")
    ln_kappa = matrix(apply(ln_kappa, 1, mean), nrow = 2)
    iter = iter+1
    gc()
  }
  mod.realm.cv_ok = all_ok
  save(mod.realm.cv, mod.realm.cv_ok, file = "tempfiles/rept_mod3.RData")
  
  
  # model with realm and svc
  # model with realm and svc
  all_ok = FALSE
  iter = 1
  while(sum(all_ok) < 5  & iter <= 2) {
    print(iter)
    if(iter == 1) {
      # run model first time
      if(parallel) {plan(multisession)}
      mod.realm.svc.cv = sdmTMB_cv(update(f, ~ . + (1|realm)), 
                                   data = dat,
                                   mesh = mesh,
                                   spatial = "on",
                                   reml = T, 
                                   spatial_varying = ~ 0 + canopy_height,
                                   fold_ids = folds,
                                   k_folds = length(unique(folds)),
                                   family = Beta(),
                                   control = sdmTMBcontrol(eval.max = 6000, iter.max = 3000),
                                   parallel = parallel)
      if(parallel) {plan(sequential)}
    } else {
      # if model is not ok and has to run a second time
      # set ln_phi to the estimation of ln_phi from the previous model
      #plan(multisession)
      mod.realm.svc.cv = sdmTMB_cv(update(f, ~ . + (1|realm)), 
                                   data = dat,
                                   mesh = mesh,
                                   spatial = "on",
                                   reml = T, 
                                   spatial_varying = ~ 0 + canopy_height,
                                   fold_ids = folds,
                                   k_folds = length(unique(folds)),
                                   family = Beta(),
                                   control = sdmTMBcontrol(eval.max = 8000, iter.max = 4000,
                                                           start = list(ln_phi = ln_phi,
                                                                        ln_tau_O = ln_tau_o,
                                                                        ln_kappa = ln_kappa)),
                                   parallel = parallel)
      if(parallel) {plan(sequential)}
    }
    
    all_ok = lapply(mod.realm.svc.cv$models, sanity)
    all_ok = sapply(all_ok, "[[", "all_ok")
    pars = lapply(mod.realm.svc.cv$models, get_pars)
    ln_phi = mean(sapply(pars, "[[", "ln_phi"))
    ln_tau_o = mean(sapply(pars, "[[", "ln_tau_O"))
    ln_kappa = sapply(pars, "[[", "ln_kappa")
    ln_kappa = matrix(apply(ln_kappa, 1, mean), nrow = 2)
    iter = iter+1
    gc()
  }
  mod.realm.svc.cv_ok = all_ok
  
  save(mod.realm.svc.cv, mod.realm.svc.cv_ok, file = "tempfiles/rept_mod4.RData")
  
  # since some of the folds did not converge, get average fold logLik across the models that did converge
  mod.cv_ll = sum(mod.cv$fold_loglik[mod.cv_ok])/sum(mod.cv_ok)
  mod.realm.cv_ll = sum(mod.realm.cv$fold_loglik[mod.realm.cv_ok])/sum(mod.realm.cv_ok)
  mod.svc.cv_ll = sum(mod.svc.cv$fold_loglik[mod.svc.cv_ok])/sum(mod.svc.cv_ok)
  mod.realm.svc.cv_ll = sum(mod.realm.svc.cv$fold_loglik[mod.realm.svc.cv_ok])/sum(mod.realm.svc.cv_ok)
  
  comp_ll = data.frame(model = c("mod", "mod.realm", "mod.svc", "mod.realm.svc"),
                       avg_logLik = c(mod.cv_ll, mod.realm.cv_ll, mod.svc.cv_ll, mod.realm.svc.cv_ll)) %>% 
    arrange(desc(avg_logLik)) %>% 
    mutate(taxon = taxon, response_var = response_var)
  
  return(list(
    mods = list(mod.cv = mod.cv, mod.realm.cv = mod.realm.cv, mod.svc.cv = mod.svc.cv, mod.realm.svc.cv = mod.realm.svc.cv),
    all_ok = c(mod.cv = mod.cv_ok, mod.realm.cv = mod.realm.cv_ok, mod.svc.cv = mod.svc.cv_ok, mod.realm.svc.cv = mod.realm.svc.cv_ok),
    compMods_cv = comp_ll
  ))
  
}

# binomial distribution for proportion arboreality
compare_cv_binom = function(f, dat, mesh, folds, parallel = F, taxon, response_var) {
  # f: formula for sdmTMB model
  # dat = dataframe for sdmTMB model
  # mesh = mesh used for sdmTMB model
  # folds for cross validation
  # logical: whether or not to run in parallel - default is false
  # taxon: character, taxon being modeled - used for labelling output
  # response_var: character; response variable for the model - used for labelling output
  
  # model with no svc or random intercept
  all_ok = FALSE
  iter = 1
  while(sum(all_ok) < 5  & iter <= 2) {
    print(iter)
    if(iter == 1) {
      if(parallel) {plan(multisession)}
      mod.cv = sdmTMB_cv(f, 
                         data = dat,
                         mesh = mesh,
                         spatial = "on",
                         reml = T, 
                         family = binomial(),
                         fold_ids = folds,
                         k_folds = length(unique(folds)),
                         control = sdmTMBcontrol(eval.max = 6000, iter.max = 3000),
                         parallel = parallel,
                         use_initial_fit = T)
      if(parallel) {plan(sequential)}
    } else {
      if(parallel) {plan(multisession)}
      mod.cv = sdmTMB_cv(f, 
                         data = dat,
                         mesh = mesh,
                         spatial = "on",
                         reml = T, 
                         fold_ids = folds,
                         k_folds = length(unique(folds)),
                         control = sdmTMBcontrol(eval.max = 6000, iter.max = 3000,
                                                 start = list(ln_phi = ln_phi, 
                                                              ln_tau_O = ln_tau_o,
                                                              ln_kappa = ln_kappa)),
                         future_globals = c("ln_phi", "ln_tau_o", "ln_kappa"),
                         parallel = parallel,
                         family = binomial(),
                         use_initial_fit = T)
      if(parallel) {plan(sequential)}
    }
    all_ok = lapply(mod.cv$models, sanity)
    all_ok = sapply(all_ok, "[[", "all_ok")
    pars = lapply(mod.cv$models, get_pars)
    ln_phi = mean(sapply(pars, "[[", "ln_phi"))
    ln_tau_o = mean(sapply(pars, "[[", "ln_tau_O"))
    ln_kappa = sapply(pars, "[[", "ln_kappa")
    ln_kappa = matrix(apply(ln_kappa, 1, mean), nrow = 2)
    iter = iter+1
    gc()
  }
  mod.cv_ok = all_ok

  all_ok = FALSE
  iter = 1
  while(sum(all_ok) < 5  & iter <= 2) {
    print(iter)
    if(iter == 1) {
      # run model first time
      if(parallel) {plan(multisession)}
      mod.svc.cv = sdmTMB_cv(f, 
                             data = dat,
                             mesh = mesh,
                             spatial = "on",
                             reml = T, 
                             spatial_varying = ~ 0 + canopy_height,
                             fold_ids = folds,
                             k_folds = length(unique(folds)),
                             family = binomial(),
                             control = sdmTMBcontrol(eval.max = 6000, iter.max = 3000),
                             parallel = parallel)
      if(parallel) {plan(sequential)}
    } else {
      # if model is not ok and has to run a second time
      # set ln_phi to the estimation of ln_phi from the previous model
      #plan(multisession)
      mod.svc.cv = sdmTMB_cv(f, 
                             data = dat,
                             mesh = mesh,
                             spatial = "on",
                             reml = T, 
                             spatial_varying = ~ 0 + canopy_height,
                             fold_ids = folds,
                             k_folds = length(unique(folds)),
                             family = binomial(),
                             control = sdmTMBcontrol(eval.max = 6000, iter.max = 3000,
                                                     start = list(ln_phi = ln_phi, 
                                                                  ln_tau_O = ln_tau_o,
                                                                  ln_kappa = ln_kappa)),
                             future_globals = c("ln_phi", "ln_tau_o", "ln_kappa"),
                             parallel = parallel)
      if(parallel) {plan(sequential)}
    }
    
    all_ok = lapply(mod.svc.cv$models, sanity)
    all_ok = sapply(all_ok, "[[", "all_ok")
    pars = lapply(mod.svc.cv$models, get_pars)
    ln_phi = mean(sapply(pars, "[[", "ln_phi"))
    ln_tau_o = mean(sapply(pars, "[[", "ln_tau_O"))
    ln_kappa = sapply(pars, "[[", "ln_kappa")
    ln_kappa = matrix(apply(ln_kappa, 1, mean), nrow = 2)
    iter = iter+1
    gc()
  }
  mod.svc.cv_ok = all_ok

  # model with realm but no svc
  all_ok = FALSE
  iter = 1
  while(sum(all_ok) < 5  & iter <= 2) {
    print(iter)
    if(iter == 1) {
      if(parallel) {plan(multisession)}
      mod.realm.cv = sdmTMB_cv(update(f, ~ . + (1|realm)), 
                               data = dat,
                               mesh = mesh,
                               spatial = "on",
                               reml = T, 
                               fold_ids = folds,
                               k_folds = length(unique(folds)),
                               family = binomial(),
                               control = sdmTMBcontrol(eval.max = 6000, iter.max = 3000),
                               parallel = parallel)
      if(parallel) {plan(sequential)}
    } else {
      if(parallel) {plan(multisession)}
      mod.realm.cv = sdmTMB_cv(update(f, ~ . + (1|realm)), 
                               data = dat,
                               mesh = mesh,
                               spatial = "on",
                               reml = T, 
                               fold_ids = folds,
                               k_folds = length(unique(folds)),
                               family = binomial(),
                               control = sdmTMBcontrol(eval.max = 6000, iter.max = 3000,
                                                       start = list(ln_phi = ln_phi,
                                                                    ln_tau_O = ln_tau_o,
                                                                    ln_kappa = ln_kappa)),
                               parallel = parallel)
      if(parallel) {plan(sequential)}
    }
    all_ok = lapply(mod.realm.cv$models, sanity)
    all_ok = sapply(all_ok, "[[", "all_ok")
    pars = lapply(mod.realm.cv$models, get_pars)
    ln_phi = mean(sapply(pars, "[[", "ln_phi"))
    ln_tau_o = mean(sapply(pars, "[[", "ln_tau_O"))
    ln_kappa = sapply(pars, "[[", "ln_kappa")
    ln_kappa = matrix(apply(ln_kappa, 1, mean), nrow = 2)
    iter = iter+1
    gc()
  }
  mod.realm.cv_ok = all_ok

  
  # model with realm and svc
  # model with realm and svc
  all_ok = FALSE
  iter = 1
  while(sum(all_ok) < 5  & iter <= 2) {
    print(iter)
    if(iter == 1) {
      # run model first time
      if(parallel) {plan(multisession)}
      mod.realm.svc.cv = sdmTMB_cv(update(f, ~ . + (1|realm)),
                                   data = dat,
                                   mesh = mesh,
                                   spatial = "on",
                                   reml = T, 
                                   spatial_varying = ~ 0 + canopy_height,
                                   fold_ids = folds,
                                   k_folds = length(unique(folds)),
                                   family = binomial(),
                                   control = sdmTMBcontrol(eval.max = 6000, iter.max = 3000),
                                   parallel = parallel)
      if(parallel) {plan(sequential)}
    } else {
      # if model is not ok and has to run a second time
      # set ln_phi to the estimation of ln_phi from the previous model
      #plan(multisession)
      mod.realm.svc.cv = sdmTMB_cv(update(f, ~ . + (1|realm)),
                                   data = dat,
                                   mesh = mesh,
                                   spatial = "on",
                                   reml = T, 
                                   spatial_varying = ~ 0 + canopy_height,
                                   fold_ids = folds,
                                   k_folds = length(unique(folds)),
                                   family = binomial(),
                                   control = sdmTMBcontrol(eval.max = 8000, iter.max = 4000,
                                                           start = list(ln_phi = ln_phi,
                                                                        ln_tau_O = ln_tau_o,
                                                                        ln_kappa = ln_kappa)),
                                   parallel = parallel)
      if(parallel) {plan(sequential)}
    }
    
    all_ok = lapply(mod.realm.svc.cv$models, sanity)
    all_ok = sapply(all_ok, "[[", "all_ok")
    pars = lapply(mod.realm.svc.cv$models, get_pars)
    ln_phi = mean(sapply(pars, "[[", "ln_phi"))
    ln_tau_o = mean(sapply(pars, "[[", "ln_tau_O"))
    ln_kappa = sapply(pars, "[[", "ln_kappa")
    ln_kappa = matrix(apply(ln_kappa, 1, mean), nrow = 2)
    iter = iter+1
    gc()
  }
  mod.realm.svc.cv_ok = all_ok
  

  # since some of the folds did not converge, get average fold logLik across the models that did converge
  mod.cv_ll = sum(mod.cv$fold_loglik[mod.cv_ok])/sum(mod.cv_ok)
  mod.realm.cv_ll = sum(mod.realm.cv$fold_loglik[mod.realm.cv_ok])/sum(mod.realm.cv_ok)
  mod.svc.cv_ll = sum(mod.svc.cv$fold_loglik[mod.svc.cv_ok])/sum(mod.svc.cv_ok)
  mod.realm.svc.cv_ll = sum(mod.realm.svc.cv$fold_loglik[mod.realm.svc.cv_ok])/sum(mod.realm.svc.cv_ok)
  
  comp_ll = data.frame(model = c("mod", "mod.realm", "mod.svc", "mod.realm.svc"),
                       avg_logLik = c(mod.cv_ll, mod.realm.cv_ll, mod.svc.cv_ll, mod.realm.svc.cv_ll)) %>% 
    arrange(desc(avg_logLik)) %>% 
    mutate(taxon = taxon, response_var = response_var)
  
  return(list(
    mods = list(mod.cv = mod.cv, mod.realm.cv = mod.realm.cv, mod.svc.cv = mod.svc.cv, mod.realm.svc.cv = mod.realm.svc.cv),
    all_ok = c(mod.cv = mod.cv_ok, mod.realm.cv = mod.realm.cv_ok, mod.svc.cv = mod.svc.cv_ok, mod.realm.svc.cv = mod.realm.svc.cv_ok),
    compMods_cv = comp_ll
  ))
  
}

# PREDICT TO FUTURE ENV DATA ----------------------------------------------
# Predicts sdmtmb model to present and future climate data
predict_future = function(mod, newdata, type, fpath) {
  # mod: model to make predictions from
  # newdata: data to predict to
  # type: "response" or "link" for predictions (see documentation in predict.sdmTMB)
  # fpath to save .RData file containing model used for predictions, and predictions to present and future climates
  
  pred.f = predict(mod, newdata = newdata, type = type)
  head(pred.f)
  
  pred = predict(mod, type = "response")
  
  pred.f$est.dif = pred.f$est - pred$est
  
  save(mod, pred, pred.f, file = fpath)
}



# DEVIANCE EXPLAINED ------------------------------------------------------
# see https://github.com/groundfish-climatechange/multispecies-redistribution-climatechange/blob/main/scripts/sdmTMB-model-functions.Rmd

calc_deviance <- function(mod, resp_var){
  # mod = sdmTMB model object
  # resp_var = character; response variable
  fit <- mod
  dat <- fit$data
  f = as.formula(paste(resp_var, "~ 1", sep = " "))
  null <- sdmTMB(vert.mean ~ 1,
                 spatial="off",
                 mesh = fit$spde,
                 data = fit$data)
  null_dev <- -2 * as.numeric(logLik(null))
  
  log_lik <- as.numeric(logLik(mod))
  
  resid_dev <- -2 * log_lik
  
  dev_explained <- 100 * (null_dev - resid_dev)/null_dev
  
  tibble(null_dev=null_dev,resid_dev=resid_dev,dev_explained=dev_explained)
}








