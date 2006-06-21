#
# randomseg.R
#
# $Revision: 1.2 $ $Date: 2006/06/16 00:43:23 $
#

rpoisline <- function(lambda, win=owin()) {
  win <- as.owin(win)
  if(win$type == "mask")
    stop("Not implemented for masks")
  # determine circumcircle
  width <- diff(win$xrange)
  height <- diff(win$yrange)
  rmax <- sqrt(width^2 + height^2)/2
  xmid <- mean(win$xrange)
  ymid <- mean(win$yrange)
  # generate poisson lines through circumcircle
  n <- rpois(1, lambda * 2 * pi * rmax)
  if(n == 0)
    return(psp(numeric(0), numeric(0), numeric(0), numeric(0),
               window=win))
  theta <- runif(n, max= 2 * pi)
  p <- runif(n, max=rmax)
  # compute intersection points with circle
  q <- sqrt(rmax^2 - p^2)
  co <- cos(theta)
  si <- sin(theta)
  X <- psp(x0= xmid + p * co + q * si,
           y0= ymid + p * si - q * co,
           x1= xmid + p * co - q * si,
           y1= ymid + p * si + q * co,
           window=win)
  # clip to window
  X <- X[, win]
  return(X)
}
