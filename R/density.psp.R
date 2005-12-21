#
#
#  density.psp.R
#
#  $Revision: 1.2 $    $Date: 2005/12/19 13:51:49 $
#
#

density.psp <- function(x, sigma, ..., edge=TRUE) {
  verifyclass(x, "psp")
  if(missing(sigma))
    sigma <- 0.1 * diameter(x$window)
  len <- lengths.psp(x)
  w <- as.mask(x$window, ...)
  xx <- as.vector(raster.x(w))
  yy <- as.vector(raster.y(w))
  for(i in seq(x$n)) {
    en <- x$ends[i,]
    coz <- (en$x1 - en$x0)/len[i]
    zin <- (en$y1 - en$y0)/len[i]
    dx <- xx - en$x0
    dy <- yy - en$y0
    u1 <- dx * coz + dy * zin
    u2 <- - dx * zin + dy * coz
    value <- dnorm(u2, sd=sigma) *
             (pnorm(u1, sd=sigma) - pnorm(u1-len[i], sd=sigma))
    totvalue <- if(i == 1) value else (value + totvalue)
  }
  dens <- im(totvalue, w$xcol, w$yrow)
  if(edge) {
    edg <- second.moment.calc(midpoints.psp(x), sigma, what="edge", ...)
    dens <- eval.im(dens/edg)
  }
  return(dens)
}
