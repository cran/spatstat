#
#  density.ppp.R
#
#  Method for 'density' for point patterns
#
#  $Revision: 1.35 $    $Date: 2011/06/17 01:43:21 $
#

ksmooth.ppp <- function(x, sigma, ..., edge=TRUE) {
  .Deprecated("density.ppp", package="spatstat")
  density.ppp(x, sigma, ..., edge=edge)
}

density.ppp <- function(x, sigma=NULL, ...,
                        weights=NULL, edge=TRUE, varcov=NULL,
                        at="pixels", leaveoneout=TRUE,
                        adjust=1, diggle=FALSE) {
  verifyclass(x, "ppp")

  output <- pickoption("output location type", at,
                       c(pixels="pixels",
                         points="points"))
  
  ker <- resolve.2D.kernel(sigma, ..., varcov=varcov, x=x, adjust=adjust)
  sigma <- ker$sigma
  varcov <- ker$varcov

  if(output == "points") {
    # VALUES AT DATA POINTS ONLY
    result <- densitypointsEngine(x, sigma, varcov=varcov,
                                  weights=weights, edge=edge,
                                  leaveoneout=leaveoneout,
                                  diggle=diggle)
    if(!is.null(uhoh <- attr(result, "warnings"))) {
      switch(uhoh,
             underflow=warning("underflow due to very small bandwidth"),
             warning(uhoh))
    }
    return(result)
  }
  
  # VALUES AT PIXELS
  if(!edge) {
    # no edge correction
    edg <- NULL
    raw <- second.moment.calc(x, sigma, what="smooth", ...,
                              weights=weights, varcov=varcov)
    raw$v <- raw$v/(raw$xstep * raw$ystep)
    smo <- raw
  } else if(!diggle) {
    # edge correction e(u)
    both <- second.moment.calc(x, sigma, what="smoothedge", ...,
                              weights=weights, varcov=varcov)
    raw <- both$smooth
    edg <- both$edge
    raw$v <- raw$v/(raw$xstep * raw$ystep)
    smo <- eval.im(raw/edg)
  } else {
    # edge correction e(x_i)
    edg <- second.moment.calc(x, sigma, what="edge", ..., varcov=varcov)
    wi <- 1/safelookup(edg, x, warn=FALSE)
    wi[!is.finite(wi)] <- 0
    # edge correction becomes weight attached to points
    if(is.null(weights)) {
      newweights <- wi
    } else {
      stopifnot(length(weights) == npoints(x))
      newweights <- weights * wi
    }
    raw <- second.moment.calc(x, sigma, what="smooth", ...,
                              weights=newweights, varcov=varcov)
    raw$v <- raw$v/(raw$xstep * raw$ystep)
    smo <- raw
  }
  
  result <- smo[x$window, drop=FALSE]

  # internal use only
  spill <- list(...)$spill
  if(!is.null(spill)) 
    return(list(sigma=sigma, varcov=varcov, raw = raw, edg=edg))

  # normal return
  attr(result, "sigma") <- sigma
  attr(result, "varcov") <- varcov
  return(result)
}

densitypointsEngine <- function(x, sigma, ...,
                                weights=NULL, edge=TRUE, varcov=NULL,
                                leaveoneout=TRUE, diggle=FALSE) {
  if(is.null(varcov)) {
    const <- 1/(2 * pi * sigma^2)
  } else {
    detSigma <- det(varcov)
    Sinv <- solve(varcov)
    const <- 1/(2 * pi * sqrt(detSigma))
  }
  # Leave-one-out computation
  # contributions from pairs of distinct points
  # closer than 8 standard deviations
  sd <- if(is.null(varcov)) sigma else sqrt(sum(diag(varcov)))
  cutoff <- 8 * sd
  # detect very small bandwidth
  if(min(nndist(x)) > cutoff) {
    npts <- npoints(x)
    result <- if(leaveoneout) rep(0, npts) else rep(const, npts)
    attr(result, "sigma") <- sigma
    attr(result, "varcov") <- varcov
    attr(result, "warnings") <- "underflow"
    return(result)
  }
  # evaluate edge correction weights at points 
  if(edge) {
    win <- x$window
    if(is.null(varcov) && win$type == "rectangle") {
      # evaluate Gaussian probabilities directly
      xr <- win$xrange
      yr <- win$yrange
      xx <- x$x
      yy <- x$y
      xprob <-
        pnorm(xr[2], mean=xx, sd=sigma) - pnorm(xr[1], mean=xx, sd=sigma)
      yprob <-
        pnorm(yr[2], mean=yy, sd=sigma) - pnorm(yr[1], mean=yy, sd=sigma)
      edgeweight <- xprob * yprob
    } else {
      edg <- second.moment.calc(x, sigma=sigma, what="edge", varcov=varcov)
      edgeweight <- safelookup(edg, x, warn=FALSE)
    }
    if(diggle) {
      # Diggle edge correction
      # edgeweight is attached to each point
      if(is.null(weights)) {
        weights <- 1/edgeweight
      } else {
        stopifnot(length(weights) == npoints(x) || length(weights) == 1)
        weights <- weights/edgeweight
      }
    }
  }
  
  if(spatstat.options("densityC")) {
    # .................. new C code ...........................
    npts <- npoints(x)
    result <- numeric(npts)
    # sort into increasing order of x coordinate (required by C code)
    oo <- order(x$x)
    xx <- x$x[oo]
    yy <- x$y[oo]
    if(is.null(varcov)) {
      # isotropic kernel
      if(is.null(weights)) {
        zz <- .C("denspt",
                 nxy     = as.integer(npts),
                 x       = as.double(xx),
                 y       = as.double(yy),
                 rmaxi   = as.double(cutoff),
                 sig     = as.double(sd),
                 result   = as.double(double(npts)),
                 PACKAGE = "spatstat")
        result[oo] <- zz$result
      } else {
        wtsort <- weights[oo]
        zz <- .C("wtdenspt",
                 nxy     = as.integer(npts),
                 x       = as.double(xx),
                 y       = as.double(yy),
                 rmaxi   = as.double(cutoff),
                 sig     = as.double(sd),
                 weight  = as.double(wtsort),
                 result   = as.double(double(npts)),
                 PACKAGE = "spatstat")
        result[oo] <- zz$result
      }
    } else {
      # anisotropic kernel
      flatSinv <- as.vector(t(Sinv))
      if(is.null(weights)) {
        zz <- .C("adenspt",
                 nxy     = as.integer(npts),
                 x       = as.double(xx),
                 y       = as.double(yy),
                 rmaxi   = as.double(cutoff),
                 detsigma = as.double(detSigma),
                 sinv    = as.double(flatSinv),
                 result   = as.double(double(npts)),
                 PACKAGE = "spatstat")
        result[oo] <- zz$result
      } else {
        wtsort <- weights[oo]
        zz <- .C("awtdenspt",
                 nxy     = as.integer(npts),
                 x       = as.double(xx),
                 y       = as.double(yy),
                 rmaxi   = as.double(cutoff),
                 detsigma = as.double(detSigma),
                 sinv    = as.double(flatSinv),
                 weight  = as.double(wtsort),
                 result   = as.double(double(npts)),
                 PACKAGE = "spatstat")
        result[oo] <- zz$result
      }
    }
  } else {
      # ..... interpreted code .........................................
    close <- closepairs(x, cutoff)
    i <- close$i
    j <- close$j
    d <- close$d
    # evaluate contribution from each close pair (i,j)
    if(is.null(varcov)) {
      contrib <- const * exp(-d^2/(2 * sigma^2))
    } else {
      # anisotropic kernel
      dx <- close$dx
      dy <- close$dy
      contrib <- const * exp(-(dx * (dx * Sinv[1,1] + dy * Sinv[1,2])
                               + dy * (dx * Sinv[2,1] + dy * Sinv[2,2]))/2)
    }
    # multiply by weights
    if(!is.null(weights))
      contrib <- contrib * weights[j]
    # sum
    result <- tapply(contrib, factor(i, levels=1:(x$n)), sum)
    result[is.na(result)] <- 0
    #
  }
  # ----- contribution from point itself ----------------
  if(!leaveoneout) {
    # add contribution from point itself
    self <- const
    if(!is.null(weights))
      self <- self * weights
    result <- result + self
  }
  # ........  Edge correction ........................................
  if(edge && !diggle) 
    result <- result/edgeweight

  # ............. validate .................................
  result <- as.numeric(result)
  npts <- npoints(x)
  if(length(result) != npts) 
    stop(paste("Internal error: incorrect number of lambda values",
               "in leave-one-out method:",
               "length(lambda) = ", length(result),
               "!=", npts, "= npoints"))
  if(any(is.na(result))) {
    nwrong <- sum(is.na(result))
    stop(paste("Internal error:", nwrong, "NA or NaN",
               ngettext(nwrong, "value", "values"),
               "generated in leave-one-out method"))
  }
  # tack on bandwidth
  attr(result, "sigma") <- sigma
  attr(result, "varcov") <- varcov
  # 
  return(result)
}

smooth.ppp <- function(X, ..., weights=rep(1, npoints(X)), at="pixels") {
  verifyclass(X, "ppp")
  if(!is.marked(X, dfok=TRUE))
    stop("X should be a marked point pattern")
  at <- pickoption("output location type", at,
                   c(pixels="pixels",
                     points="points"))
  # determine smoothing parameters
  ker <- resolve.2D.kernel(..., x=X)
  sigma <- ker$sigma
  varcov <- ker$varcov
  #
  if(weightsgiven <- !missing(weights)) {
    check.nvector(weights, npoints(X)) 
    # rescale weights to avoid numerical gremlins
    weights <- weights/median(abs(weights))
  }
  # get marks
  marx <- marks(X, dfok=TRUE)
  #
  if(!is.data.frame(marx)) {
    # ........ vector of marks ...................
    values <- marx
    if(is.factor(values)) {
      warning("Factor valued marks were converted to integers")
      values <- as.numeric(values)
    }
    switch(at,
           points={
             if(!weightsgiven)
               weights <- NULL
             result <-
               do.call("smoothpointsEngine",
                       resolve.defaults(list(x=X),
                                        list(values=values, weights=weights),
                                        list(sigma=sigma, varcov=varcov),
                                        list(...)))
           },
           pixels={
             numerator <-
               do.call("density.ppp",
                       resolve.defaults(list(x=X, at="pixels"),
                                        list(weights = values * weights),
                                        list(sigma=sigma, varcov=varcov),
                                        list(...),
                                        list(edge=FALSE)))
             denominator <-
               do.call("density.ppp",
                       resolve.defaults(list(x=X, at="pixels"),
                                        list(weights = weights),
                                        list(sigma=sigma, varcov=varcov),
                                        list(...),
                                        list(edge=FALSE)))
             ratio <- eval.im(numerator/denominator)
             result <- eval.im(ifelse(is.finite(ratio), ratio, NA))
             attr(result, "warnings") <- attr(numerator, "warnings")
           })
  } else {
    # ......... data frame of marks ..................
    nmarx <- ncol(marx)
    result <- list()
    # compute denominator
    denominator <-
      do.call("density.ppp",
              resolve.defaults(list(x=X, at=at),
                               list(weights = weights),
                               list(sigma=sigma, varcov=varcov),
                               list(...),
                               list(edge=FALSE)))
    # smooth each column of marks in turn
    for(j in 1:nmarx) {
      values <- marx[,j] 
      if(is.factor(values)) {
        warning(paste("Factor valued marks in column", j,
                      "were converted to integers"))
        values <- as.numeric(values)
      }
      # compute j-th numerator
      numerator <-
        do.call("density.ppp",
                resolve.defaults(list(x=X, at=at),
                                 list(weights = values * weights),
                                 list(sigma=sigma, varcov=varcov),
                                 list(...),
                                 list(edge=FALSE)))
      switch(at,
             points={
               if(is.null(uhoh <- attr(numerator, "warnings"))) {
                 ratio <- numerator/denominator
                 ratio <- ifelse(is.finite(ratio), ratio, NA)
               } else {
                 warning("returning original values")
                 ratio <- values
                 attr(ratio, "warnings") <- uhoh
               }
             },
             pixels={
               ratio <- eval.im(numerator/denominator)
               ratio <- eval.im(ifelse(is.finite(ratio), ratio, NA))
               attr(ratio, "warnings") <- attr(numerator, "warnings")
           })
    # store results
      result[[j]] <- ratio
    }
    # wrap up
    names(result) <- colnames(marx)
    result <- switch(at,
                     pixels=as.listof(result),
                     points=as.data.frame(result))
    attr(result, "warnings") <-
      unlist(lapply(result, function(x){ attr(x, "warnings") }))
  }
  attr(result, "sigma") <- sigma
  attr(result, "varcov") <- varcov
  return(result)
}

smoothpointsEngine <- function(x, values, sigma, ...,
                               weights=NULL, varcov=NULL,
                               leaveoneout=TRUE) {
  if(is.null(varcov)) {
    const <- 1/(2 * pi * sigma^2)
  } else {
    detSigma <- det(varcov)
    Sinv <- solve(varcov)
    const <- 1/(2 * pi * sqrt(detSigma))
  }
  # Contributions from pairs of distinct points
  # closer than 8 standard deviations
  sd <- if(is.null(varcov)) sigma else sqrt(sum(diag(varcov)))
  cutoff <- 8 * sd
  # detect very small bandwidth
  if(min(nndist(x)) > cutoff) {
    npts <- npoints(x)
    warning("Very small bandwidth; original values returned")
    result <- values
    attr(result, "sigma") <- sigma
    attr(result, "varcov") <- varcov
    attr(result, "warnings") <- "underflow"
    return(result)
  }
  if(spatstat.options("densityC")) {
    # .................. new C code ...........................
    npts <- npoints(x)
    result <- numeric(npts)
    stopifnot(is.logical(leaveoneout))
    # sort into increasing order of x coordinate (required by C code)
    oo <- order(x$x)
    xx <- x$x[oo]
    yy <- x$y[oo]
    vv <- values[oo]
    if(is.null(varcov)) {
      # isotropic kernel
      if(is.null(weights)) {
        zz <- .C("smoopt",
                 nxy     = as.integer(npts),
                 x       = as.double(xx),
                 y       = as.double(yy),
                 v       = as.double(vv),
                 self    = as.integer(!leaveoneout),
                 rmaxi   = as.double(cutoff),
                 sig     = as.double(sd),
                 result   = as.double(double(npts)),
                 PACKAGE = "spatstat")
        result[oo] <- zz$result
      } else {
        wtsort <- weights[oo]
        zz <- .C("wtsmoopt",
                 nxy     = as.integer(npts),
                 x       = as.double(xx),
                 y       = as.double(yy),
                 v       = as.double(vv),
                 self    = as.integer(!leaveoneout),
                 rmaxi   = as.double(cutoff),
                 sig     = as.double(sd),
                 weight  = as.double(wtsort),
                 result   = as.double(double(npts)),
                 PACKAGE = "spatstat")
        result[oo] <- zz$result
      }
    } else {
      # anisotropic kernel
      flatSinv <- as.vector(t(Sinv))
      if(is.null(weights)) {
        zz <- .C("asmoopt",
                 nxy     = as.integer(npts),
                 x       = as.double(xx),
                 y       = as.double(yy),
                 v       = as.double(vv),
                 self    = as.integer(!leaveoneout),
                 rmaxi   = as.double(cutoff),
                 detsigma = as.double(detSigma),
                 sinv    = as.double(flatSinv),
                 result   = as.double(double(npts)),
                 PACKAGE = "spatstat")
        result[oo] <- zz$result
      } else {
        wtsort <- weights[oo]
        zz <- .C("awtsmoopt",
                 nxy     = as.integer(npts),
                 x       = as.double(xx),
                 y       = as.double(yy),
                 v       = as.double(vv),
                 self    = as.integer(!leaveoneout),
                 rmaxi   = as.double(cutoff),
                 detsigma = as.double(detSigma),
                 sinv    = as.double(flatSinv),
                 weight  = as.double(wtsort),
                 result   = as.double(double(npts)),
                 PACKAGE = "spatstat")
        result[oo] <- zz$result
      }
    }
  } else {
    # previous, partly interpreted code
    # compute weighted densities
    if(is.null(weights)) {
      # weights are implicitly equal to 1
      numerator <- do.call("density.ppp",
                         resolve.defaults(list(x=x, at="points"),
                                          list(weights = values),
                                          list(sigma=sigma, varcov=varcov),
                                          list(leaveoneout=leaveoneout),
                                          list(...),
                                          list(edge=FALSE)))
      denominator <- do.call("density.ppp",
                             resolve.defaults(list(x=x, at="points"),
                                              list(sigma=sigma, varcov=varcov),
                                              list(leaveoneout=leaveoneout),
                                              list(...),
                                              list(edge=FALSE)))
    } else {
      numerator <- do.call("density.ppp",
                           resolve.defaults(list(x=x, at="points"),
                                            list(weights = values * weights),
                                            list(sigma=sigma, varcov=varcov),
                                            list(leaveoneout=leaveoneout),
                                            list(...),
                                            list(edge=FALSE)))
      denominator <- do.call("density.ppp",
                             resolve.defaults(list(x=x, at="points"),
                                              list(weights = weights),
                                              list(sigma=sigma, varcov=varcov),
                                              list(leaveoneout=leaveoneout),
                                              list(...),
                                              list(edge=FALSE)))
    }
    if(is.null(uhoh <- attr(numerator, "warnings"))) {
      result <- numerator/denominator
      result <- ifelse(is.finite(result), result, NA)
    } else {
      warning("returning original values")
      result <- values
      attr(result, "warnings") <- uhoh
    }
  }
  # pack up and return
  attr(result, "sigma") <- sigma
  attr(result, "varcov") <- varcov
  return(result)
}


markmean <- function(X, ...) { smooth.ppp(X, ...) }

markvar  <- function(X, ...) {
  if(!is.marked(X, dfok=FALSE))
    stop("X should have (one column of) marks")
  E1 <- smooth.ppp(X, ...)
  X2 <- X %mark% marks(X)^2
  E2 <- smooth.ppp(X2, ...)
  V <- eval.im(E2 - E1^2)
  return(V)
}

resolve.2D.kernel <- function(sigma=NULL, ..., varcov=NULL, x, mindist=NULL,
                              adjust=1) {
  sigma.given <- !is.null(sigma)
  varcov.given <- !is.null(varcov)
  if(sigma.given) {
    stopifnot(is.numeric(sigma))
    stopifnot(length(sigma) %in% c(1,2))
    stopifnot(all(sigma > 0))
  }
  if(varcov.given)
    stopifnot(is.matrix(varcov) && nrow(varcov) == 2 && ncol(varcov)==2 )
  # reconcile
  ngiven <- varcov.given + sigma.given
  switch(ngiven+1,
         {
           # default
           w <- x$window
           sigma <- (1/8) * min(diff(w$xrange), diff(w$yrange))
         },
         {
           if(sigma.given && length(sigma) == 2) 
             varcov <- diag(sigma^2)
           if(!is.null(varcov))
             sigma <- NULL
         },
         {
           stop(paste("Give only one of the arguments",
                      sQuote("sigma"), "and", sQuote("varcov")))
         })
  # apply adjustments 
  if(!is.null(sigma))  sigma <- adjust * sigma
  if(!is.null(varcov)) varcov <- sqrt(adjust) * sigma
  #
  sd <- if(is.null(varcov)) sigma else sqrt(sum(diag(varcov)))
  cutoff <- 8 * sd
  uhoh <- if(!is.null(mindist) && mindist < cutoff) "underflow" else NULL
  result <- list(sigma=sigma, varcov=varcov, cutoff=cutoff, warnings=uhoh)
  return(result)
}

