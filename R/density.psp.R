#
#
#  density.psp.R
#
#  $Revision: 1.4 $    $Date: 2010/06/07 06:36:57 $
#
#

density.psp <- function(x, sigma, ..., edge=TRUE) {
  verifyclass(x, "psp")
  w <- x$window
  if(missing(sigma))
    sigma <- 0.1 * diameter(w)
  w <- as.mask(w, ...)
  len <- lengths.psp(x)
  if(x$n == 0 || all(len == 0))
    return(as.im(0, w))
  #
  ang <- angles.psp(x)
  coz <- cos(ang)
  zin <- sin(ang)
  xx <- as.vector(raster.x(w))
  yy <- as.vector(raster.y(w))
  # compute matrix contribution from each segment 
  for(i in seq(x$n)) {
    en <- x$ends[i,]
    dx <- xx - en$x0
    dy <- yy - en$y0
    u1 <- dx * coz[i] + dy * zin[i]
    u2 <- - dx * zin[i] + dy * coz[i]
    value <- dnorm(u2, sd=sigma) *
      (pnorm(u1, sd=sigma) - pnorm(u1-len[i], sd=sigma))
    totvalue <- if(i == 1) value else (value + totvalue)
  }
  dens <- im(totvalue, w$xcol, w$yrow)
  if(edge) {
    edg <- second.moment.calc(midpoints.psp(x), sigma, what="edge", ...)
    dens <- eval.im(dens/edg)
  }
  dens <- dens[x$window, drop=FALSE]
  return(dens)
}
