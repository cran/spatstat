#
#	ripras.S	Ripley-Rasson estimator of domain
#
#
#	$Revision: 1.1 $	$Date: 2002/04/07 09:15:39 $
#
#
#
#
#-------------------------------------
ripras <- function(x, y=NULL) {
  # extract x, y coordinates
  z <- xy.coords(x, y)
  x <- z$x
  y <- z$y
  # convex hull
  h <- rev(chull(x, y))  # must be anticlockwise
  # temporary window
  w <- owin(poly=list(x=x[h], y=y[h]))
  # centroid
  ce <- centroid.owin(w)
  # expansion factor
  f <- 1/sqrt(1 - length(h)/length(x))
  # new polygon
  xp <- f * (x[h] - ce$x) + ce$x
  yp <- f * (y[h] - ce$y) + ce$y
  W <- owin(poly=list(x=xp, y=yp))
  return(W)
}
