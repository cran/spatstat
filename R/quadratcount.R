#
#  quadratcount.R
#
#  $Revision: 1.2 $  $Date: 2005/04/19 01:18:53 $
#

quadratcount <- function(X, nx=5, ny=nx, xbreaks, ybreaks)  {
  verifyclass(X, "ppp")

  xr <- X$window$xrange
  yr <- X$window$yrange

  if(missing(xbreaks))
    xbreaks <- seq(xr[1], xr[2], length=nx+1)
  else if(min(xbreaks) > xr[1] || max(xbreaks) < xr[2])
    stop("xbreaks do not span the range of x coordinates in the window")

  if(missing(ybreaks))
    ybreaks <- seq(xr[1], xr[2], length=nx+1)
  else if(min(ybreaks) > xr[1] || max(ybreaks) < yr[2])
    stop("ybreaks do not span the range of y coordinates in the window")

  xg <- cut(X$x, breaks=xbreaks, include.lowest=TRUE)
  yg <- cut(X$y, breaks=ybreaks, include.lowest=TRUE)
  sumz <- table(list(x=xg, y=yg))
#  attr(sumz, "xbreaks") <- xbreaks
#  attr(sumz, "ybreaks") <- ybreaks
  return(sumz)
}
