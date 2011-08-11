#
# observed Fisher information matrix
# and asymptotic covariance & correlation matrices
# for (inhom) Poisson models
#
#  $Revision: 1.29 $  $Date: 2011/08/05 01:53:48 $
#

vcov.ppm <- function(object, ..., what="vcov", verbose=TRUE,
                     gam.action="warn", matrix.action="warn",
                     hessian=FALSE) {
  verifyclass(object, "ppm")

  stopifnot(length(what) == 1 && is.character(what))
  what.options <- c("vcov", "corr", "fisher", "Fisher", "internals")
  what.map     <- c("vcov", "corr", "fisher", "fisher", "internals")
  if(is.na(m <- pmatch(what, what.options)))
    stop(paste("Unrecognised option: what=", sQuote(what)))
  what <- what.map[m]

  method <- resolve.defaults(list(...), list(method="C"))$method
  
  # Fisher information *may* be contained in object
  fisher <- object$fisher
  varcov <- object$varcov
  
  # Do we need to go into the guts?
  needguts <-
    (is.null(fisher) && what=="fisher") ||
    (is.null(varcov) && what %in% c("vcov", "corr")) ||
    (what == "internals")

  # In general it is not true that varcov = solve(fisher)
  # because we might use different estimators,
  # or the parameters might be a subset of the canonical parameter

  ############## guts ##############################

  if(needguts) {
    if(!is.poisson.ppm(object) && !hessian)
      stop(paste("Sorry, asymptotic variances are not available",
                 "for non-Poisson models",
                 "fitted using maximum pseudolikelihood"))

    gf <- getglmfit(object)
    # we need a glm or gam
    if(is.null(gf)) {
      if(verbose) 
        warning("Refitting the model using GLM/GAM")
      object <- update(object, forcefit=TRUE)
      gf <- getglmfit(object)
      if(is.null(gf))
        stop("Internal error - refitting did not yield a glm object")
    }
    # if it's a gam, issue a warning
    if(inherits(gf, "gam")) {
      switch(gam.action,
             fatal={
               stop(paste("model was fitted by gam();",
                          "execution halted because fatal=TRUE"),
                    call.=FALSE)
             },
             warn={
               warning(paste("model was fitted by gam();",
                             "asymptotic variance calculation ignores this"),
                       call.=FALSE)
             },
             silent={})
    }
    # compute related stuff
    fi <- fitted(object, type="trend", check=FALSE, drop=TRUE)
    wt <- quad.ppm(object, drop=TRUE)$w
    # extract model matrix
    mom <- model.matrix(object, keepNA=FALSE) 
    # compute Fisher information if not known
    if(is.null(fisher)) {
      switch(method,
             C = {
               fisher <- sumouter(mom, fi * wt)
             },
             interpreted = {
               for(i in 1:nrow(mom)) {
                 ro <- mom[i, ]
                 v <- outer(ro, ro, "*") * fi[i] * wt[i]
                 if(!any(is.na(v)))
                   fisher <- fisher + v
               }
               momnames <- dimnames(mom)[[2]]
               dimnames(fisher) <- list(momnames, momnames)
             })
    }
  }

  ############## end of guts ####################################
  
  if(what %in% c("vcov", "corr") && is.null(varcov)) {
    # Need variance-covariance matrix.
    # Derive from Fisher information
    msilent <- (matrix.action == "silent")
    mfatal  <- (matrix.action == "fatal")
    varcov <- try(solve(fisher), silent=msilent)
    if(inherits(varcov, "try-error")) {
      if(mfatal)
        stop("Inverse of matrix not available")
      return(NULL)
    }
  }
         
  switch(what,
         fisher = { return(fisher) },
         vcov   = {
           return(varcov)
         },
         corr={
           sd <- sqrt(diag(varcov))
           return(varcov / outer(sd, sd, "*"))
         },
         internals = {
           return(list(fisher=fisher, suff=mom))
         })
}


suffloc <- function(object) {
  return(vcov.ppm(object, what="internals")$suff)
}

