#
# simulation of FITTED model
#
#  $Revision: 1.8 $ $Date: 2004/08/30 05:38:41 $
#
#
rmh.ppm <- function(model, start = NULL, control = NULL, ...,
                    verbose=TRUE, project=TRUE) {
  verifyclass(model, "ppm")

  # convert fitted model object to list of parameters for rmh.default
  X <- rmhmodel.ppm(model, verbose=verbose, project=project)

  # call appropriate simulation routine

  if(X$cif != "poisson") {
    if(is.null(start))
      start <- list(n.start=summary(model, quick=TRUE)$data$n)
    if(is.null(control)) 
      control <- list(nrep=1e6)
    return(rmh.default(X, start=start, control=control, ..., verbose=verbose))
  }
  
  # Poisson process
  intensity <- if(is.null(X$trend)) X$par$beta else X$trend
  if(is.null(X$types))
    return(rpoispp(intensity, win=X$w, ...))
  else
    return(rmpoispp(intensity, win=X$w, types=X$types))
}

