#
#    plot.ppm.S
#
#    $Revision: 2.3 $    $Date: 2006/04/12 11:56:26 $
#
#    plot.ppm()
#         Plot a point process model fitted by ppm().
#        
#
#
plot.ppm <- function(x, ngrid = c(40,40),
		     superimpose = TRUE,
                     trend = TRUE, cif = TRUE, pause = interactive(),
                     how=c("persp","image", "contour"),
                     plot.it=TRUE,
                     locations=NULL, covariates=NULL, ...)
{
  model <- x
#       Plot a point process model fitted by ppm().
#
  verifyclass(model, "ppm")
#
#       find out what kind of model it is
#
  mod <- summary(model)
  stationary <- mod$stationary
  poisson    <- mod$poisson
  marked     <- mod$marked
  multitype  <- mod$multitype
  data       <- mod$entries$data
        
  if(marked) {
    if(!multitype)
      stop("Not implemented for general marked point processes")
    else
      mrkvals <- levels(data$marks)
  } else mrkvals <- 1
  ntypes <- length(mrkvals)
        
#
#        Interpret options
#        -----------------
#        
#        Whether to plot trend, cif
        
  if(!trend && !cif) {
    cat("Nothing plotted - both \'trend\' and \'cif\' are FALSE\n")
    return(invisible(NULL))
  }
#        Suppress uninteresting plots
#        unless explicitly instructed otherwise
  if(missing(trend))
    trend <- !stationary
  if(missing(cif))
    cif <- !poisson
        
  if(!trend && !cif) {
    cat("Nothing plotted -- all plots selected are flat surfaces.\n")
    return(invisible())
  }

#
#        Do the prediction
#        ------------------

  out <- list()
  surftypes <- c("trend","cif")[c(trend,cif)]
  ng <- if(missing(ngrid) && !missing(locations)) NULL else ngrid

  for (ttt in surftypes) {
    p <- predict(model,
                   ngrid=ng, locations=locations, covariates=covariates,
                   type = ttt)
    if(is.im(p))
      p <- list(p)
    out[[ttt]] <- p
  }

#        Make it a plotppm object
#        ------------------------  
  
  class(out) <- "plotppm"
  attr(out, "mrkvals") <- mrkvals

#        Actually plot it if required
#        ----------------------------  
  if(plot.it) {
    if(!superimpose)
      data <- NULL
    plot(out,data=data,trend=trend,cif=cif,how=how,pause=pause, ...)
  }

  
  return(invisible(out)) 
}

