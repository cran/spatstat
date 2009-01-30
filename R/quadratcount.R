#
#  quadratcount.R
#
#  $Revision: 1.22 $  $Date: 2009/01/30 00:39:28 $
#

quadratcount <- function(X, ...) {
  UseMethod("quadratcount")
}

quadratcount.splitppp <- function(X, ...) {
  as.listof(lapply(X, quadratcount, ...))
}

quadratcount.ppp <- function(X, nx=5, ny=nx, ...,
                             xbreaks=NULL, ybreaks=NULL,
                             tess=NULL)  {
  verifyclass(X, "ppp")
  W <- X$window

  # determine tessellation
  if(is.null(tess)) {
    # rectangular boundaries 
    if(!is.numeric(nx))
      stop("nx should be numeric")
    tess <- quadrats(X, nx=nx, ny=ny, xbreaks=xbreaks, ybreaks=ybreaks)
  } else if(!inherits(tess, "tess"))
    stop("The argument tess should be a tessellation", call.=FALSE)

  if(tess$type == "rect") {
    # rectangular quadrats - fast code
    Xcount <- rectquadrat.countEngine(X$x, X$y, tess$xgrid, tess$ygrid)
  } else {
    # quadrats are another type of tessellation
    Y <- cut(X, tess)
    if(any(is.na(marks(Y))))
      warning("Tessellation does not contain all the points of X")
    Xcount <- table(tile=marks(Y))
  } 
  attr(Xcount, "tess") <- tess
  class(Xcount) <- c("quadratcount", class(Xcount))
  return(Xcount)
}

plot.quadratcount <- function(x, ...,
                              add=FALSE, entries=as.vector(t(as.table(x))),
                              dx=0, dy=0) {
  xname <- deparse(substitute(x))
  tess <- attr(x, "tess")
  do.call("plot.tess",
          resolve.defaults(list(tess),
                           list(...),
                           list(main=xname, add=add),
                           .StripNull=TRUE))
  if(!is.null(entries)) {
    labels <- paste(as.vector(entries))
    til <- tiles(tess)
    incircles <- lapply(til, incircle)
    x0 <- unlist(lapply(incircles, function(z) { z$x }))
    y0 <- unlist(lapply(incircles, function(z) { z$y }))
    ra <- unlist(lapply(incircles, function(z) { z$r }))
    do.call.matched("text.default",
                    resolve.defaults(list(x=x0 + dx * ra, y = y0 + dy * ra),
                                     list(labels=labels),
                                     list(...)))
  }
  return(invisible(NULL))
}

rectquadrat.breaks <- function(xr, yr, nx=5, ny=nx, xbreaks=NULL, ybreaks=NULL) {
  if(is.null(xbreaks))
    xbreaks <- seq(xr[1], xr[2], length=nx+1)
  else if(min(xbreaks) > xr[1] || max(xbreaks) < xr[2])
    stop("xbreaks do not span the range of x coordinates in the window")
  if(is.null(ybreaks))
    ybreaks <- seq(yr[1], yr[2], length=ny+1)
  else if(min(ybreaks) > yr[1] || max(ybreaks) < yr[2])
    stop("ybreaks do not span the range of y coordinates in the window")
  return(list(xbreaks=xbreaks, ybreaks=ybreaks))
}

rectquadrat.countEngine <- function(x, y, xbreaks, ybreaks, weights) {
  if(min(x) < min(xbreaks) || max(x) > max(xbreaks))
    stop("xbreaks do not span the actual range of x coordinates in data")
  if(min(y) < min(ybreaks) || max(y) > max(ybreaks))
    stop("ybreaks do not span the actual range of y coordinates in data")
  xg <- cut(x, breaks=xbreaks, include.lowest=TRUE)
  yg <- cut(y, breaks=ybreaks, include.lowest=TRUE)
  if(missing(weights)) 
    sumz <- table(list(y=yg, x=xg))
  else {
    sumz <- tapply(weights, list(y=yg, x=xg), sum)
    if(any(nbg <- is.na(sumz)))
      sumz[nbg] <- 0
  }
  # reverse order of y 
  sumz <- sumz[rev(seq(nrow(sumz))), ]
  sumz <- as.table(sumz)
  #
  attr(sumz, "xbreaks") <- xbreaks
  attr(sumz, "ybreaks") <- ybreaks
  return(sumz)
}

quadrats <- function(X, nx=5, ny=nx, xbreaks = NULL, ybreaks = NULL) {
  W <- as.owin(X)
  xr <- W$xrange
  yr <- W$yrange
  b <- rectquadrat.breaks(xr, yr, nx, ny, xbreaks, ybreaks)
  Z <- tess(xgrid=b$xbreaks, ygrid=b$ybreaks)
  if(W$type != "rectangle")
    Z <- intersect.tess(Z, W)
  return(Z)
}

