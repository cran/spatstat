#
# simulation of FITTED model
#
#  $Revision: 1.5 $ $Date: 2003/03/14 07:17:08 $
#
rmh.ppm <- function(model, start, control, ...) {
  verifyclass(model, "ppm")

  if(!is.stationary.ppm(model))
    stop("Simulation of non-stationary models is not yet implemented")
  
  #### POISSON CASE must be treated separately
  # as the model does not contain an 'interaction' object
  # and rmh() does not simulate Poisson processes anyway.
  
  if(is.poisson.ppm(model)) {
    lambda <- exp(model$coef[[1]])
    X <- model$Q$data
    if(!is.marked(X)) {
      # uniform Poisson, unmarked
      return(rpoispp(lambda, win=X$window))
    } else {
      if(!is.factor(X$marks))
        stop("Internal error: marks are not a factor")
      # uniform Poisson, multitype
      mu <- markspace.integral(X)
      Xsim <- rpoispp(lambda * mu, win=X$window)
      if(Xsim$n > 0) {
        lev <- levels(X$marks)
        marques <- sample(lev, Xsim$n, replace=TRUE)
        Xsim$marks <- factor(marques, levels=lev)
      }
      return(Xsim)
    }
  }

  ######  GENERAL non-Poisson CASE #####################################
  # extract the  translator function for the interaction in the model
  inte <- model$interaction
  if(is.null(inte$rmhmodel))
    stop(paste("Simulation of a fitted \'", inte$name,
               "\' has not yet been implemented in rmh"))
  # apply the translation
  dmodel <- (inte$rmhmodel)(model, inte)
  if(is.null(dmodel))
    stop("Internal problem: rmhmodel returns NULL")
  # invoke Metropolis-Hastings engine
  rmh.default(dmodel, start, control, ...)
}



