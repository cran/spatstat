#
# simulation of FITTED model
#
#  $Revision: 1.15 $ $Date: 2005/03/10 21:30:55 $
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

