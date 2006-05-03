#
# observed Fisher information matrix
# and asymptotic covariance & correlation matrices
# for (inhom) Poisson models
#
#

vcov.ppm <- function(object, ..., what="vcov", verbose=TRUE) {
  verifyclass(object, "ppm")
  if(!is.poisson.ppm(object))
    stop("Sorry, vcov.ppm is only implemented for Poisson processes")

  what.options <- c("vcov", "corr", "fisher", "Fisher")
  what.map     <- c("vcov", "corr", "fisher", "fisher")
  if(is.na(m <- pmatch(what, what.options)))
    stop(paste("Unrecognised option: what=", sQuote(what)))
  what <- what.map[m]
  
  fi <- fitted(object)
  wt <- quad.ppm(object)$w
  gf <- getglmfit(object)
  if(is.null(gf)) {
    if(verbose) 
      warning("Refitting the model using GLM/GAM")
    object <- update(object, forcefit=TRUE)
    gf <- getglmfit(object)
    if(is.null(gf))
      stop("Internal error - refitting did not yield a glm object")
  }
  gd <- object$internal$glmdata
  fo <- object$trend
  if(is.null(fo)) fo <- (~1)
  mof <- model.frame(fo, gd)
  mom <- model.matrix(fo, mof)
  fisher <- 0
  for(i in 1:nrow(mom)) {
    ro <- mom[i, ]
    v <- outer(ro, ro, "*") * fi[i] * wt[i]
    fisher <- fisher + v
  }
  momnames <- dimnames(mom)[[2]]
  dimnames(fisher) <- list(momnames, momnames)

  switch(what,
         fisher = { return(fisher) },
         vcov   = { return(solve(fisher)) },
         corr={
           co <- solve(fisher)
           sd <- sqrt(diag(co))
           return(co / outer(sd, sd, "*"))
         })
}
