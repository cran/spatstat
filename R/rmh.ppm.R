#
# simulation of FITTED model
#
#  $Revision: 1.7 $ $Date: 2004/01/27 11:24:46 $
#
#
rmh.ppm <- function(model, start = NULL, control = NULL, ...,
                    verbose=TRUE, project=TRUE) {
  verifyclass(model, "ppm")

  # convert fitted model object to list of parameters for rmh.default
  X <- rmhmodel.ppm(model, verbose=verbose, project=project)

  # call appropriate simulation routine

  if(X$cif != "poisson")
    return(rmh.default(X, start=start, control=control, ..., verbose=verbose))

  # Poisson process
  intensity <- if(is.null(X$trend)) X$par$beta else X$trend
  if(is.null(X$types))
    return(rpoispp(intensity, win=X$w, ...))
  else
    return(rmpoispp(intensity, win=X$w, types=X$types))
}

