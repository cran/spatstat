#
# observed Fisher information matrix
# and asymptotic covariance & correlation matrices
# for (inhom) Poisson models
#
#  $Revision: 1.18 $  $Date: 2007/04/02 06:04:51 $
#

vcov.ppm <- function(object, ..., what="vcov", verbose=TRUE, gamaction="warn") {
  verifyclass(object, "ppm")

  stopifnot(length(what) == 1 && is.character(what))
  what.options <- c("vcov", "corr", "fisher", "Fisher", "internals")
  what.map     <- c("vcov", "corr", "fisher", "fisher", "internals")
  if(is.na(m <- pmatch(what, what.options)))
    stop(paste("Unrecognised option: what=", sQuote(what)))
  what <- what.map[m]

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
    if(!is.poisson.ppm(object))
      stop(paste("Sorry, vcov.ppm is not implemented for non-Poisson models",
                 "fitted using maximum pseudolikelihood"))

    fi <- fitted(object, type="trend", check=FALSE)
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
    if(!inherits(gf, "gam")) {
      # model fitted by glm
      gd <- object$internal$glmdata
      fo <- object$trend
      if(is.null(fo)) fo <- (~1)
    } else {
      # model fitted by gam
      switch(gamaction,
             fatal={
               stop("Sorry, vcov.ppm is not implemented for models fitted by gam")
             },
             warn={
               warning(paste("model was fitted by gam();",
                             "asymptotic variance calculation ignores this"))
               # proceed to Rolf's code
             },
             silent={})
      # Rolf's code
      gd <- as.data.frame(predict(gf, type="terms"))
      names(gd) <- gsub("\\)", "", gsub("\\(", "", names(gd)))
      rhs <- paste(colnames(gd), collapse="+")
      fo <- as.formula(paste("~", rhs))
      the.rest <- object$internal$glmdata[, c(".mpl.W", ".mpl.Y", ".mpl.SUBSET")]
      gd <- cbind(the.rest, gd)
    }
    mof <- model.frame(fo, gd)
    mom <- model.matrix(fo, mof)
    
    # compute Fisher information if not known
    if(is.null(fisher)) {
      fisher <- 0
      for(i in 1:nrow(mom)) {
        ro <- mom[i, ]
        v <- outer(ro, ro, "*") * fi[i] * wt[i]
        if(!any(is.na(v)))
          fisher <- fisher + v
      }
      momnames <- dimnames(mom)[[2]]
      dimnames(fisher) <- list(momnames, momnames)
    }
  }

  ############## end of guts ####################################
  

  switch(what,
         fisher = { return(fisher) },
         vcov   = {
           # Derive variance-covariance from Fisher information unless given
           if(is.null(varcov)) {
             varcov <- try(solve(fisher))
             if(inherits(varcov, "try-error"))
               return(NULL)
           }
           return(varcov)
         },
         corr={
           # Derive variance-covariance from Fisher information unless given
           if(is.null(varcov)) {
             varcov <- try(solve(fisher))
             if(inherits(varcov, "try-error"))
               return(NULL)
           }
           sd <- sqrt(diag(varcov))
           return(varcov / outer(sd, sd, "*"))
         },
         internals = {
           nbg <- matrowany(is.na(gd))
           if(any(nbg)) {
             suff <- matrix( , nrow=nrow(gd), ncol=ncol(mom))
             suff[nbg,] <- NA
             suff[!nbg, ] <- mom
           } else suff <- mom
           return(list(fisher=fisher, suff=suff))
         })
}


suffloc <- function(object) {
  return(vcov.ppm(object, what="internals")$suff)
}

