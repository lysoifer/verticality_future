# from DHARMa - fixed plot title

plotResiduals2 <- function (simulationOutput, form = NULL, quantreg = NULL, rank = T, 
          asFactor = NULL, smoothScatter = NULL, quantiles = c(0.25, 
                                                               0.5, 0.75), ...) 
{
  a <- list(...)
  a$ylab = checkDots("ylab", "DHARMa residual", ...)
  a$xlab = checkDots("xlab", ifelse(is.null(form), "Model predictions", 
                                    gsub(".*[$]", "", deparse(substitute(form)))), ...)
  if (rank == T) 
    a$xlab = paste(a$xlab, "(rank transformed)")
  simulationOutput = ensureDHARMa(simulationOutput, convert = T)
  res = simulationOutput$scaledResiduals
  if (inherits(form, "DHARMa")) 
    stop("DHARMa::plotResiduals > argument form cannot be of class DHARMa. Note that the syntax of plotResiduals has changed since DHARMa 0.3.0. See ?plotResiduals.")
  pred = ensurePredictor(simulationOutput, form)
  if (!is.factor(pred)) {
    if (rank == T) {
      pred = rank(pred, ties.method = "average")
      pred = pred/max(pred)
      a$xlim = checkDots("xlim", c(0, 1), ...)
    }
    nuniq = length(unique(pred))
    ndata = length(pred)
    if (is.null(asFactor)) 
      asFactor = (nuniq == 1) | (nuniq < 10 & ndata/nuniq > 
                                   10)
    if (asFactor) 
      pred = factor(pred)
  }
  if (is.null(quantreg)) 
    if (length(res) > 2000) 
      quantreg = FALSE
  else quantreg = TRUE
  switchScatter = 10000
  if (is.null(smoothScatter)) 
    if (length(res) > switchScatter) 
      smoothScatter = TRUE
  else smoothScatter = FALSE
  blackcol = rgb(0, 0, 0, alpha = max(0.1, 1 - 3 * length(res)/switchScatter))
  if (is.factor(pred)) {
    testCategorical(simulationOutput = simulationOutput, 
                    catPred = pred, quantiles = quantiles)
  }
  else if (smoothScatter == TRUE) {
    defaultCol = ifelse(res == 0 | res == 1, 2, blackcol)
    do.call(graphics::smoothScatter, append(list(x = pred, 
                                                 y = res, ylim = c(0, 1), axes = FALSE, colramp = colorRampPalette(c("white", 
                                                                                                                     "darkgrey"))), a))
    points(pred[defaultCol == 2], res[defaultCol == 2], 
           col = "red", cex = 0.5)
    axis(1)
    axis(2, at = c(0, quantiles, 1))
  }
  else {
    defaultCol = ifelse(res == 0 | res == 1, 2, blackcol)
    defaultPch = ifelse(res == 0 | res == 1, 8, 1)
    a$col = checkDots("col", defaultCol, ...)
    a$pch = checkDots("pch", defaultPch, ...)
    do.call(plot, append(list(res ~ pred, ylim = c(0, 1), 
                              axes = FALSE), a))
    axis(1)
    axis(2, at = c(0, quantiles, 1))
  }
  main = checkDots("main", ifelse(is.null(form), "Residual vs. predicted", 
                                  "Residual vs. predictor"), ...)
  out = NULL
  if (is.numeric(pred)) {
    if (quantreg == F) {
      #title(main = main, cex.main = 1)
      abline(h = quantiles, col = "black", lwd = 0.5, 
             lty = 2)
      try({
        lines(smooth.spline(pred, res, df = 10), lty = 2, 
              lwd = 2, col = "red")
        abline(h = 0.5, col = "red", lwd = 2)
      }, silent = T)
    }
    else {
      out = testQuantiles(simulationOutput, pred, quantiles = quantiles, 
                          plot = F)
      if (any(out$pvals < 0.05, na.rm = TRUE)) {
        main = paste(main, "Quantile deviations detected (red curves)", 
                     sep = "\n")
        if (out$p.value <= 0.05) {
          main = paste(main, "Combined adjusted quantile test significant", 
                       sep = "\n")
        }
        else {
          main = paste(main, "Combined adjusted quantile test n.s.", 
                       sep = "\n")
        }
        maincol = "red"
      }
      else {
        main = paste(main, "No significant problems detected", 
                     sep = "\n")
        maincol = "black"
      }
      title(main = main, cex.main = 0.8, col.main = maincol)
      for (i in 1:length(quantiles)) {
        lineCol = ifelse(out$pvals[i] <= 0.05 & !(is.na(out$pvals[i])), 
                         "red", "black")
        filCol = ifelse(out$pvals[i] <= 0.05 & !(is.na(out$pvals[i])), 
                        "#FF000040", "#00000020")
        abline(h = quantiles[i], col = lineCol, lwd = 0.5, 
               lty = 2)
        polygon(c(out$predictions$pred, rev(out$predictions$pred)), 
                c(out$predictions[, 2 * i] - out$predictions[, 
                                                             2 * i + 1], rev(out$predictions[, 2 * i] + 
                                                                               out$predictions[, 2 * i + 1])), col = "#00000020", 
                border = F)
        lines(out$predictions$pred, out$predictions[, 
                                                    2 * i], col = lineCol, lwd = 2)
      }
    }
  }
  invisible(out)
}

checkDots <- function(name, default, ...) {
  args <- list(...)
  if(!name %in% names(args)) {
    ## Default value
    return(default)
  } else {
    ## If the argument was defined in the ... part, return it
    return(args[[name]])
  }
}

ensureDHARMa <- function(simulationOutput,
                         convert = FALSE){
  
  if(inherits(simulationOutput, "DHARMa")){
    return(simulationOutput)
  } else {
    
    if(convert == FALSE) stop("wrong argument to function, simulationOutput must be a DHARMa object!")
    else {
      
      if (class(simulationOutput)[1] %in% getPossibleModels()){
        if (convert == "Model" | convert == TRUE) return(simulateResiduals(simulationOutput))
      } else if(is.vector(simulationOutput, mode = "numeric") & convert == TRUE) {
        out = list()
        out$scaledResiduals = simulationOutput
        out$nObs = length(out$scaledResiduals)
        class(out) = "DHARMa"
        return(out)
      }
    }
  }
  stop("wrong argument to function, simulationOutput must be a DHARMa object or a numeric vector of quantile residuals!")
}


ensurePredictor <- function(simulationOutput,
                            predictor = NULL){
  if(!is.null(predictor)){
    
    if(length(predictor) != length(simulationOutput$scaledResiduals)) stop("DHARMa: residuals and predictor do not have the same length. The issue is possibly that you have NAs in your predictor that were removed during the model fit. Remove the NA values from your predictor.")
    
    if(is.character(predictor)) {
      predictor = factor(predictor)
      warning("DHARMa:::ensurePredictor: character string was provided as predictor. DHARMa has converted to factor automatically. To remove this warning, please convert to factor before attempting to plot with DHARMa.")
    }
    
  } else {
    
    predictor = simulationOutput$fittedPredictedResponse
    if(is.null(predictor)) stop("DHARMa: can't extract predictor from simulationOutput, and no predictor provided")
  }
  return(predictor)
}
