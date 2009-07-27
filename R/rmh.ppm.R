#
# simulation of FITTED model
#
#  $Revision: 1.18 $ $Date: 2009/07/21 22:12:29 $
#
#
rmh.ppm <- function(model, start = NULL,
                    control = rmhcontrol(),
                    ..., verbose=TRUE, project=TRUE) {
  verifyclass(model, "ppm")

  control <- rmhcontrol(control)
  
  # convert fitted model object to list of parameters for rmh.default
  X <- rmhmodel.ppm(model, verbose=verbose, project=project, control=control)

  # set initial state

  if(is.null(start)) {
    datapattern <- data.ppm(model)
    start <- rmhstart(n.start=datapattern$n)
  }

  return(rmh.default(X, start=start, control=control, ..., verbose=verbose))
}

simulate.ppm <- function(object, nsim=1, ...,
                         start = NULL,
                         control = rmhcontrol(),
                         project=TRUE,
                         verbose=FALSE, progress=TRUE) {
  verifyclass(object, "ppm")

  # Set up parameters for rmh
  rmodel <- rmhmodel(object, verbose=FALSE)
  if(is.null(start)) {
    datapattern <- data.ppm(object)
    start <- rmhstart(n.start=datapattern$n)
  }
  rstart <- rmhstart(start)
  rcontr <- rmhcontrol(control)
  # pre-digest arguments
  rmhinfolist <- rmh(rmodel, rstart, rcontr, preponly=TRUE, verbose=verbose)
  # go
  out <- list()
  if(nsim > 0) {
    for(i in 1:nsim) {
      if(progress) progressreport(i, nsim)
      out[[i]] <- rmhEngine(rmhinfolist, verbose=verbose, ...)
    }
  }
  out <- as.listof(out)
  names(out) <- paste("Simulation", 1:nsim)
  return(out)
}  
