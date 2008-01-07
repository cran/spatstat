#
#  quadratcount.R
#
#  $Revision: 1.11 $  $Date: 2007/12/19 17:53:02 $
#

quadratcount <- function(X, nx=5, ny=nx, xbreaks=NULL, ybreaks=NULL)  {
  verifyclass(X, "ppp")

  W <- X$window
  xr <- W$xrange
  yr <- W$yrange
  b <- quadrat.breaks(xr, yr, nx, ny, xbreaks, ybreaks)
  Xcount <- quadrat.countEngine(X$x, X$y, b$xbreaks, b$ybreaks)
  attr(Xcount, "window") <- W
  class(Xcount) <- c("quadratcount", class(Xcount))
  return(Xcount)
}

plot.quadratcount <- function(x, ..., add=FALSE, entries=as.table(x), dx=0, dy=0) {
  xname <- deparse(substitute(x))
  W <- attr(x, "window")
  if(!add)
    do.call("plot.owin",
            resolve.defaults(list(W),
                             list(col=NULL),
                             list(...),
                             list(main=xname),
                             .StripNull=TRUE))
  xbk <- attr(x, "xbreaks")
  ybk <- attr(x, "ybreaks")
  xr <- W$xrange
  yr <- W$yrange
  do.call.matched("segments",
                  resolve.defaults(list(x0=xbk, y0=yr[1], x1=xbk, y1=yr[2]),
                                   list(...)))
  do.call.matched("segments",
                  resolve.defaults(list(x0=xr[1], y0=ybk, x1=xr[2], y1=ybk),
                                   list(...)))
  if(!is.null(entries)) {
    labels <- paste(as.vector(entries))
    xmid <- xbk[-1] - diff(xbk) * (1/2 - dx)
    ymid <- ybk[-1] - diff(ybk) * (1/2 - dy)
    xy <- expand.grid(x=xmid, y=ymid)
    do.call.matched("text.default",
                    resolve.defaults(list(x=xy$x, y = xy$y),
                                     list(labels=labels),
                                     list(...)))
  }
  invisible(NULL)
}



quadrat.breaks <- function(xr, yr, nx=5, ny=nx, xbreaks=NULL, ybreaks=NULL) {
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

quadrat.countEngine <- function(x, y, xbreaks, ybreaks, weights) {
  if(min(x) < min(xbreaks) || max(x) > max(xbreaks))
    stop("xbreaks do not span the actual range of x coordinates in data")
  if(min(y) < min(ybreaks) || max(y) > max(ybreaks))
    stop("ybreaks do not span the actual range of y coordinates in data")
  xg <- cut(x, breaks=xbreaks, include.lowest=TRUE)
  yg <- cut(y, breaks=ybreaks, include.lowest=TRUE)
  if(missing(weights)) 
    sumz <- table(list(x=xg, y=yg))
  else {
    sumz <- tapply(weights, list(x=xg, y=yg), sum)
    if(any(nbg <- is.na(sumz)))
      sumz[nbg] <- 0
  }
  attr(sumz, "xbreaks") <- xbreaks
  attr(sumz, "ybreaks") <- ybreaks
  return(sumz)
}
