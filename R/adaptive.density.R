#
#  adaptive.density.R
#
#  $Revision: 1.2 $   $Date: 2008/10/01 00:50:28 $
#
#

adaptive.density <- function(X, f=0.1, ..., nrep=1) {
  stopifnot(is.ppp(X))
  npoints <- X$n
  stopifnot(is.numeric(f) && length(f) == 1 && f > 0 & f < 1)
  ntess <- floor(f * npoints)
  if(ntess == 0) {
    # naive estimate of intensity
    W <- X$window
    lam <- npoints/area.owin(W)
    return(as.im(lam, W, ...))
  }
  if(nrep > 1) {
    # estimate is the average of nrep randomised estimates
    total <- 0
    for(i in seq(nrep)) {
      estimate <- adaptive.density(X, f, ..., nrep=1)
      total <- eval.im(total + estimate)
    }
    average <- eval.im(total/nrep)
    return(average)
  }
  ncount <- npoints - ntess
  fcount <- ncount/npoints
  itess <- sample(seq(npoints), ntess, replace=FALSE)
  Xtess <- X[itess]
  Xcount <- X[-itess]
  tes <- dirichlet(Xtess)
  meanintensity <- function(x) { x$n/area.owin(x$window) }
  lam <- unlist(lapply(split(Xcount, tes), meanintensity))
  tesim <- as.im(tes, ...)
  out <- eval.im(lam[as.integer(tesim)]/fcount)
  return(out)
}
