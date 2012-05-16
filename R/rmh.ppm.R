#
# simulation of FITTED model
#
#  $Revision: 1.26 $ $Date: 2012/05/12 06:46:42 $
#
#
rmh.ppm <- function(model, start = NULL,
                    control = default.rmhcontrol(model),
                    ..., project=TRUE, verbose=TRUE) {
  verifyclass(model, "ppm")
  argh <- list(...)

  if(is.null(control)) {
    control <- default.rmhcontrol(model)
  } else {
    control <- rmhcontrol(control)
  }

  # override 
  if(length(list(...)) > 0)
    control <- update(control, ...)
  
  # convert fitted model object to list of parameters for rmh.default
  X <- rmhmodel(model, verbose=verbose, project=project, control=control)

  # set initial state

  if(is.null(start)) {
    datapattern <- data.ppm(model)
    start <- rmhstart(n.start=datapattern$n)
  }
  
  # call rmh.default 
  # passing only arguments unrecognised by rmhcontrol
  known <- names(argh) %in% names(formals(rmhcontrol.default))
  fargs <- argh[!known]

  Y <- do.call("rmh.default",
               append(list(model=X, start=start, control=control,
                           verbose=verbose),
                      fargs))
  return(Y)
}

simulate.ppm <- function(object, nsim=1, ...,
                         start = NULL,
                         control = default.rmhcontrol(object),
                         project=TRUE,
                         verbose=FALSE,
                         progress=(nsim > 1)) {
  verifyclass(object, "ppm")
  argh <- list(...)

  # set up control parameters
  if(missing(control) || is.null(control)) {
    rcontr <- default.rmhcontrol(object)
  } else {
    rcontr <- rmhcontrol(control)
  }
  # override 
  if(length(list(...)) > 0)
    rcontr <- update(rcontr, ...)
  
  # Set up model parameters for rmh
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
    # call rmh
    # passing only arguments unrecognised by rmhcontrol
    known <- names(argh) %in% names(formals(rmhcontrol.default))
    fargs <- argh[!known]
    rmhargs <- append(list(InfoList=rmhinfolist, verbose=verbose), fargs)
    for(i in 1:nsim) {
      out[[i]] <- do.call("rmhEngine", rmhargs)
      if(progress) progressreport(i, nsim)
    }
  }
  out <- as.listof(out)
  if(nsim > 0)
    names(out) <- paste("Simulation", 1:nsim)
  return(out)
}  
