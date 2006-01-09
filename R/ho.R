#
#  ho.R
#
#  Huang-Ogata method 
#
#  $Revision: 1.2 $ $Date: 2006/01/09 05:14:55 $
#

ho.engine <- function(model, ..., nsim=100, nrmh=1e5,
                        start=NULL,
                        control=list(nrep=nrmh), verb=TRUE) {
  verifyclass(model, "ppm")

  if(is.null(start)) 
    start <- list(n.start=data.ppm(model)$n)
  
  # check that the model can be simulated
  inte <- model$interaction
  valid <- is.null(inte) || inte$valid(coef(model), inte)
  if(!valid) {
    warning("Fitted model cannot be simulated")
    return(NULL)
  }
  
  # compute the observed value of the sufficient statistic
  X <- data.ppm(model)
  sobs <- suffstat(model, X)
  
  # generate 'nsim' realisations of the fitted model
  # and compute the sufficient statistics of the model
  rmodel <- rmhmodel(model)
  rstart <- rmhstart(start)
  rcontrol <- rmhcontrol(control)
  if(verb) cat("Simulating... ")
  for(i in 1:nsim) {
    if(verb) cat(paste(i, " ", if(i %% 10 == 0) "\n", sep=""))
    Xi <- rmh(rmodel, rstart, rcontrol, verbose=FALSE)
    v <- suffstat(model,Xi)
    if(i == 1) 
      svalues <- matrix(, nrow=nsim, ncol=length(v))
    svalues[i, ] <- v
  }
  if(verb) cat("Done.\n\n")
  # calculate the sample mean and variance of the
  # sufficient statistic for the simulations
  smean <- apply(svalues, 2, mean, na.rm=TRUE)
  svar <- var(svalues, na.rm=TRUE)
  # value of canonical parameter from MPL fit
  theta0 <- coef(model)
  # Newton-Raphson update
  theta <- theta0 + as.vector(solve(svar) %*% (sobs - smean))
  # update model  & return
  newmodel <- model
  newmodel$coef <- theta
  newmodel$coef.mpl <- theta0
  newmodel$method <- "ho"
  return(newmodel)
}


