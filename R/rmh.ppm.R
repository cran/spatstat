#
# simulation of FITTED model
#
#  $Revision: 1.22 $ $Date: 2011/11/02 07:33:15 $
#
#
rmh.ppm <- function(model, start = NULL,
                    control = default.rmhcontrol(model, expand=expand),
                    ..., project=TRUE, expand=NULL, verbose=TRUE) {
  verifyclass(model, "ppm")

  control <- rmhcontrol(control)
  
  # convert fitted model object to list of parameters for rmh.default
  X <- rmhmodel(model, verbose=verbose, project=project, control=control)

  # set initial state

  if(is.null(start)) {
    datapattern <- data.ppm(model)
    start <- rmhstart(n.start=datapattern$n)
  }

  return(rmh.default(X, start=start, control=control, ..., verbose=verbose))
}

simulate.ppm <- function(object, nsim=1, ...,
                         start = NULL,
                         control = default.rmhcontrol(object, expand=expand),
                         project=TRUE,
                         expand=NULL,
                         verbose=FALSE,
                         progress=(nsim > 1)) {
  verifyclass(object, "ppm")

  # Set up parameters for rmh
  rcontr <- rmhcontrol(control)
  rmodel <- rmhmodel(object, verbose=FALSE, project=TRUE, control=rcontr)
  if(is.null(start)) {
    datapattern <- data.ppm(object)
    start <- rmhstart(n.start=datapattern$n)
  }
  rstart <- rmhstart(start)
  # pre-digest arguments
  rmhinfolist <- rmh(rmodel, rstart, rcontr, preponly=TRUE, verbose=verbose)
  # go
  out <- list()
  if(nsim > 0) {
    if(progress) {
      cat(paste("Generating", nsim, "simulated", 
                ngettext(nsim, "pattern", "patterns"),
                "..."))
      flush.console()
    }
    for(i in 1:nsim) {
      out[[i]] <- rmhEngine(rmhinfolist, verbose=verbose, ...)
      if(progress) progressreport(i, nsim)
    }
  }
  out <- as.listof(out)
  if(nsim > 0)
    names(out) <- paste("Simulation", 1:nsim)
  return(out)
}  
