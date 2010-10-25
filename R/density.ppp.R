#
#  density.ppp.R
#
#  Method for 'density' for point patterns
#
#  $Revision: 1.23 $    $Date: 2010/10/24 12:41:56 $
#

ksmooth.ppp <- function(x, sigma, ..., edge=TRUE) {
  .Deprecated("density.ppp", package="spatstat")
  density.ppp(x, sigma, ..., edge=edge)
}

density.ppp <- function(x, sigma, ..., weights=NULL, edge=TRUE, varcov=NULL,
                        at="pixels", leaveoneout=TRUE) {
  verifyclass(x, "ppp")
  sigma.given <- !missing(sigma) && !is.null(sigma)
  varcov.given <- !is.null(varcov)
  if(sigma.given) {
    stopifnot(is.numeric(sigma))
    stopifnot(length(sigma) %in% c(1,2))
    stopifnot(all(sigma > 0))
  }
  if(varcov.given)
    stopifnot(is.matrix(varcov) && nrow(varcov) == 2 && ncol(varcov)==2 )    
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
      
  output <- pickoption("output location type", at,
                       c(pixels="pixels",
                         points="points"))

  if(output == "points") {
    # VALUES AT DATA POINTS ONLY
    result <- densitypointsEngine(x, sigma, varcov=varcov,
                                  weights=weights, edge=edge,
                                  leaveoneout=leaveoneout)
    return(result)
  }
  
  # VALUES AT PIXELS
  if(!edge) {
    raw <- second.moment.calc(x, sigma, what="smooth", ...,
                              weights=weights, varcov=varcov)
    raw$v <- raw$v/(raw$xstep * raw$ystep)
    smo <- raw
  } else {
    both <- second.moment.calc(x, sigma, what="smoothedge", ...,
                              weights=weights, varcov=varcov)
    raw <- both$smooth
    edg <- both$edge
    raw$v <- raw$v/(raw$xstep * raw$ystep)
    smo <- eval.im(raw/edg)
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
                                leaveoneout=TRUE) {
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
                 value   = as.double(double(npts)),
                 PACKAGE = "spatstat")
        result[oo] <- zz$value
      } else {
        wtsort <- weights[oo]
        zz <- .C("wtdenspt",
                 nxy     = as.integer(npts),
                 x       = as.double(xx),
                 y       = as.double(yy),
                 rmaxi   = as.double(cutoff),
                 sig     = as.double(sd),
                 weight  = as.double(wtsort),
                 value   = as.double(double(npts)),
                 PACKAGE = "spatstat")
        result[oo] <- zz$value
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
                 value   = as.double(double(npts)),
                 PACKAGE = "spatstat")
        result[oo] <- zz$value
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
                 value   = as.double(double(npts)),
                 PACKAGE = "spatstat")
        result[oo] <- zz$value
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
    result <- result/edgeweight
  }
  result <- as.numeric(result)
  # ............. validate .................................
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

smooth.ppp <- function(X, ..., weights=rep(1,X$n), at="pixels") {
  verifyclass(X, "ppp")
  if(!is.marked(X, dfok=FALSE))
    stop("X should be a marked point pattern")
  values <- marks(X, dfok=FALSE)
  if(is.factor(values)) {
    warning("Factor values were converted to integers")
    values <- as.numeric(values)
  }
  if(!missing(weights)) {
    # rescale weights to avoid numerical gremlins
    weights <- weights/median(abs(weights))
  }
  numerator <- density(X, ..., at=at, weights= values * weights)
  denominator <- density(X, ..., at=at, weights= weights)
  at <- pickoption("output location type", at,
                   c(pixels="pixels",
                     points="points"))
  switch(at,
         pixels={
           ratio <- eval.im(numerator/denominator)
           ratio <- eval.im(ifelse(is.infinite(ratio), NA, ratio))
         },
         points= {
           ratio <- numerator/denominator
           ratio <- ifelse(is.infinite(ratio), NA, ratio)
         })
  attr(ratio, "sigma") <- attr(numerator, "sigma")
  attr(ratio, "varcov") <- attr(numerator, "varcov")
  return(ratio)
}

markmean <- function(X, ...) { smooth.ppp(X, ...) }

markvar  <- function(X, ...) {
  E1 <- smooth.ppp(X, ...)
  X2 <- X %mark% marks(X)^2
  E2 <- smooth.ppp(X2, ...)
  V <- eval.im(E2 - E1^2)
  return(V)
}

