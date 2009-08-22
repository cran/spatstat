#
# areadiff.R
#
#  $Revision: 1.10 $  $Date: 2009/08/22 07:26:39 $
#
# Computes sufficient statistic for area-interaction process
#
# Invokes areadiff.c
#

areaGain <- function(u, X, r, W=NULL, ngrid=spatstat.options("ngrid.disc")) {
  verifyclass(X, "ppp")
  u <- as.ppp(u, W=as.owin(X))
  stopifnot(is.numeric(r) && all(is.finite(r)) && all(r >= 0))
  #
  constrain <- !is.null(W)
  if(constrain && (W$type != "rectangle"))
    stop("only implemented for W=rectangle or W=NULL")
  #
  nu <- u$n
  nr <- length(r)
  if(nr == 0)
    return(numeric(0))
  rmax <- max(r)
  #
  xx <- X$x
  yy <- X$y
  result <- matrix(, nrow=nu, ncol=nr)
  DUP <- spatstat.options("dupC")
  #
  for(i in 1:nu) {
    # shift u[i] to origin
    xu <- u$x[i]
    yu <- u$y[i]
    xshift <- xx - xu
    yshift <- yy - yu
    # find points within distance 2 rmax of origin
    close <- (xshift^2 + yshift^2 < 4 * rmax^2)
    nclose <- sum(close)
    # invoke C routine
    if(!constrain) {
      z <- .C("areadifs",
              rad = as.double(r),
              nrads = as.integer(nr),
              x   = as.double(xshift[close]),
              y   = as.double(yshift[close]),
              nn  = as.integer(nclose),
              ngrid = as.integer(ngrid),
              answer = as.double(numeric(nr)),
              DUP=DUP,
              PACKAGE="spatstat")
      result[i,] <- z$answer
    } else {
      z <- .C("areaBdif",
              rad = as.double(r),
              nrads = as.integer(nr),
              x   = as.double(xshift[close]),
              y   = as.double(yshift[close]),
              nn  = as.integer(nclose),
              ngrid = as.integer(ngrid),
              x0 = as.double(W$xrange[1] - xu),
              y0 = as.double(W$yrange[1] - yu),
              x1 = as.double(W$xrange[2] - xu),
              y1 = as.double(W$yrange[2] - yu),
              answer = as.double(numeric(nr)),
              DUP=DUP,
              PACKAGE="spatstat")
      result[i,] <- z$answer
    }
  }
  return(result)
}

areaLoss <- function(X, r, ngrid=spatstat.options("ngrid.disc"), subset=NULL) {
  verifyclass(X, "ppp")
  n <- X$n
  W <- as.owin(X)
  indices <- if(is.null(subset)) 1:n else (1:n)[subset]
  answer <- matrix(, nrow=length(indices), ncol=length(r))
  for(k in seq(indices)) {
    i <- indices[k]
    answer[k,] <- areaGain(X[i], X[-i], r, W=W, ngrid=ngrid)
  }
  return(answer)
}

