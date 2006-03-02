#
# randomseg.R
#
# $Revision$ $Date$
#

rpoisline <- function(lambda, win=owin()) {
  win <- as.owin(win)
  if(win$type != "rectangle")
    stop("Only implemented for rectangular windows")
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
  q <- sqrt(1 - p^2)
  co <- cos(theta)
  si <- sin(theta)
  X <- psp(x0= xmid + p * co + q * si,
           y0= xmid + p * si - q * co,
           x1= xmid + p * co - q * si,
           y1= xmid + p * si + q * co,
           window=win)
  # clip to window
  X <- X[, win]
  return(X)
}
