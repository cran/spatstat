#
#           Kmeasure.R
#
#           $Revision: 1.31 $    $Date: 2008/12/19 19:43:30 $
#
#     pixellate()        convert a point pattern to a pixel image
#
#     Kmeasure()         compute an estimate of the second order moment measure
#
#     Kest.fft()        use Kmeasure() to form an estimate of the K-function
#
#     second.moment.calc()    underlying algorithm
#
#     This file uses the temporary 'image' class defined in images.R

pixellate <- function(x, ..., W=NULL, weights=NULL, padzero=FALSE)
{
    verifyclass(x, "ppp")

    if(!is.null(W))
      W <- as.mask(W)
    else {
      # determine W using as.mask
      dotargs <- list(...)
      namesargs <- names(dotargs)
      matched <- namesargs %in% names(formals(as.mask))
      W <- do.call("as.mask", append(list(x$window), dotargs[matched]))
    } 

    if(x$n == 0) {
      zeroimage <- as.im(as.double(0), W)
      if(padzero) # map NA to 0
        zeroimage <- na.handle.im(zeroimage, 0)
      return(zeroimage)
    }
    
    pixels <- nearest.raster.point(x$x, x$y, W)
    nr <- W$dim[1]
    nc <- W$dim[2]
    if(is.null(weights)) {
    ta <- table(row = factor(pixels$row, levels = 1:nr), col = factor(pixels$col,
        levels = 1:nc))
    } else {
        ta <- tapply(weights, list(row = factor(pixels$row, levels = 1:nr),
                    col = factor(pixels$col, levels=1:nc)), sum)
        ta[is.na(ta)] <- 0
    }
    out <- im(ta, xcol = W$xcol, yrow = W$yrow, unitname=unitname(W))
    # clip to window of data
    if(!padzero)
      out <- out[W, drop=FALSE]
    return(out)
}

Kmeasure <- function(X, sigma, edge=TRUE, ..., varcov=NULL) {
  stopifnot(is.ppp(X))
  
  sigma.given <- !missing(sigma) && !is.null(sigma)
  varcov.given <- !is.null(varcov)
  ngiven <- sigma.given + varcov.given
  if(ngiven == 2)
    stop(paste("Give only one of the arguments",
               sQuote("sigma"), "and", sQuote("varcov")))
  if(ngiven == 0)
    stop(paste("Please specify smoothing bandwidth", sQuote("sigma"),
               "or", sQuote("varcov")))
  if(varcov.given) {
    stopifnot(is.matrix(varcov) && nrow(varcov) == 2 && ncol(varcov)==2 )
    sigma <- NULL
  } else {
    stopifnot(is.numeric(sigma))
    stopifnot(length(sigma) %in% c(1,2))
    stopifnot(all(sigma > 0))
    if(length(sigma) == 2) {
      varcov <- diag(sigma^2)
      sigma <- NULL
    }
  }  

  second.moment.calc(X, sigma=sigma, edge, "Kmeasure", varcov=varcov)
}

second.moment.calc <- function(x, sigma=NULL, edge=TRUE,
                               what="Kmeasure", debug=FALSE, ..., varcov=NULL)
{
  choices <- c("kernel", "smooth", "Kmeasure", "Bartlett", "edge")
  if(!(what %in% choices))
    stop(paste("Unknown choice: what = ", sQuote(what),
               "; available options are:",
               paste(sQuote(choices), collapse=", ")))
  if(is.ppp(x)) 
  # convert list of points to mass distribution 
    X <- pixellate(x, ..., padzero=TRUE)
  else if(is.im(x))
    X <- x
  else stop("internal error: unrecognised format for x")
  Y <- X$v
  xw <- X$xrange
  yw <- X$yrange
  # pad with zeroes
  nr <- nrow(Y)
  nc <- ncol(Y)
  Ypad <- matrix(0, ncol=2*nc, nrow=2*nr)
  Ypad[1:nr, 1:nc] <- Y
  lengthYpad <- 4 * nc * nr
  # corresponding coordinates
  xw.pad <- xw[1] + 2 * c(0, diff(xw))
  yw.pad <- yw[1] + 2 * c(0, diff(yw))
  xcol.pad <- xw[1] + X$xstep * (1/2 + 0:(2*nc-1))
  yrow.pad <- yw[1] + X$ystep * (1/2 + 0:(2*nr-1))
  # set up Gauss kernel
  xcol.G <- X$xstep * c(0:(nc-1),-(nc:1))
  yrow.G <- X$ystep * c(0:(nr-1),-(nr:1))
  xx <- matrix(xcol.G[col(Ypad)], ncol=2*nc, nrow=2*nr)
  yy <- matrix(yrow.G[row(Ypad)], ncol=2*nc, nrow=2*nr)
  if(!is.null(sigma)) {
    if(max(abs(diff(xw)),abs(diff(yw))) < 6 * sigma)
      warning("sigma is too large for this window")
    Kern <- exp(-(xx^2 + yy^2)/(2 * sigma^2))/(2 * pi * sigma^2) * X$xstep * X$ystep
  } else if(!is.null(varcov)) {
    # anisotropic kernel
    detSigma <- det(varcov)
    Sinv <- solve(varcov)
    const <- X$xstep * X$ystep/(2 * pi * sqrt(detSigma))
    Kern <- const * exp(-(xx * (xx * Sinv[1,1] + yy * Sinv[1,2])
                          + yy * (xx * Sinv[2,1] + yy * Sinv[2,2]))/2)
  } else 
    stop("Must specify either sigma or varcov")

  if(what=="kernel") {
    # return the kernel
    # first rearrange it into spatially sensible order (monotone x and y)
    rtwist <- ((-nr):(nr-1)) %% (2 * nr) + 1
    ctwist <- (-nc):(nc-1) %% (2*nc) + 1
    if(debug) {
      if(any(order(xcol.G) != rtwist))
        cat("something round the twist\n")
    }
    Kermit <- Kern[ rtwist, ctwist]
    ker <- im(Kermit, xcol.G[ctwist], yrow.G[ rtwist], unitname=unitname(x))
    return(ker)
  }
  # convolve using fft
  fY <- fft(Ypad)
  fK <- fft(Kern)
  sm <- fft(fY * fK, inverse=TRUE)/lengthYpad
  if(debug) {
    cat(paste("smooth: maximum imaginary part=", signif(max(Im(sm)),3), "\n"))
    cat(paste("smooth: mass error=", signif(sum(Mod(sm))-x$n,3), "\n"))
  }
  if(what=="smooth") {
    # return the smoothed point pattern
    smo <- im(Re(sm)[1:nr, 1:nc], xcol.pad[1:nc], yrow.pad[1:nr],
              unitname=unitname(x))
    return(smo)
  }

  bart <- Mod(fY)^2 * fK
  if(what=="Bartlett") {
     # rearrange into spatially sensible order (monotone x and y)
    rtwist <- ((-nr):(nr-1)) %% (2 * nr) + 1
    ctwist <- (-nc):(nc-1) %% (2*nc) + 1
    bart <- bart[ rtwist, ctwist]
    return(im(Mod(bart),(-nc):(nc-1), (-nr):(nr-1)))
  }
  
  mom <- fft(bart, inverse=TRUE)/lengthYpad
  if(debug) {
    cat(paste("2nd moment measure: maximum imaginary part=",
              signif(max(Im(mom)),3), "\n"))
    cat(paste("2nd moment measure: mass error=",
              signif(sum(Mod(mom))-x$n^2, 3), "\n"))
  }
  mom <- Mod(mom)
  # subtract (delta_0 * kernel) * npoints
#  browser()
  mom <- mom - x$n * Kern
  # edge correction
  if(edge) {
    # compute kernel-smoothed set covariance
    M <- as.mask(x$window, dimyx=c(nr, nc))$m
    # previous line ensures M has same dimensions and scale as Y 
    Mpad <- matrix(0, ncol=2*nc, nrow=2*nr)
    Mpad[1:nr, 1:nc] <- M
    lengthMpad <- 4 * nc * nr
    fM <- fft(Mpad)
    co <- fft(Mod(fM)^2 * fK, inverse=TRUE)/lengthMpad
    co <- Mod(co) 
    a <- sum(M)
    wt <- a/co
    me <- spatstat.options("maxedgewt")
    weight <- matrix(pmin(me, wt), ncol=2*nc, nrow=2*nr)
    if(debug) browser()
    mom <- mom * weight
  # set to NA outside 'reasonable' region
    mom[wt > 10] <- NA
  }
 # rearrange into spatially sensible order (monotone x and y)
  rtwist <- ((-nr):(nr-1)) %% (2 * nr) + 1
  ctwist <- (-nc):(nc-1) %% (2*nc) + 1
  mom <- mom[ rtwist, ctwist]
  if(debug) {
    if(any(order(xcol.G) != rtwist))
      cat("something round the twist\n")
  }
  if(what=="edge") {
    # return convolution of window with kernel
    # (evaluated inside window only)
    con <- fft(fM * fK, inverse=TRUE)/lengthMpad
    return(Mod(con[1:nr, 1:nc]))
  }
  # divide by number of points * lambda
  mom <- mom * area.owin(x$window) / x$n^2
  # return it
  mm <- im(mom, xcol.G[ctwist], yrow.G[rtwist], unitname=unitname(x))
  return(mm)
}

Kest.fft <- function(X, sigma, r=NULL, breaks=NULL) {
  verifyclass(X, "ppp")
  W <- X$window
  lambda <- X$n/area.owin(W)
  rmaxdefault <- rmax.rule("K", W, lambda)        
  bk <- handle.r.b.args(r, breaks, W, rmaxdefault=rmaxdefault)
  breaks <- bk$val
  rvalues <- bk$r
  u <- Kmeasure(X, sigma)
  xx <- rasterx.im(u)
  yy <- rastery.im(u)
  rr <- sqrt(xx^2 + yy^2)
  tr <- whist(rr, breaks, u$v)
  K  <- cumsum(tr)
  rmax <- min(rr[is.na(u$v)])
  K[rvalues >= rmax] <- NA
  result <- data.frame(r=rvalues, theo=pi * rvalues^2, border=K)
  w <- X$window
  alim <- c(0, min(diff(w$xrange), diff(w$yrange))/4)
  out <- fv(result,
            "r", substitute(Kfft(r), NULL),
            "border",
              cbind(border, theo) ~ r, alim,
              c("r", "Kpois(r)", "Kbord(r)"),
              c("distance argument r",
                "theoretical Poisson K(r)",
                "border-corrected estimate of K(r)"),
              unitname=unitname(X)
            )
  return(out)
}

