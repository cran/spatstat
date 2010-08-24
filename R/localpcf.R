#
#   localpcf.R
#
#  $Revision: 1.4 $  $Date: 2010/08/24 06:36:37 $
#
#

localpcfmatrix <- function(X, i=seq(npoints(X)), delta, rmax,
                           ..., nr=512, adjust=1) {
  missi <- missing(i)
  nX <- npoints(X)
  W <- as.owin(X)
  lambda <- nX/area.owin(W)
  if(missing(delta)) 
    delta <- adjust * 0.15/sqrt(lambda)
  if(missing(rmax)) 
    rmax <- rmax.rule("K", W, lambda)
  # sort points in increasing order of x coordinate
  oX <- order(X$x)
  Xsort <- X[oX]
  idXsort <- (1:nX)[oX]
  if(missi) {
    Y <- X
    oY <- oX
    Ysort   <- Xsort
    idYsort <- idXsort
  } else {
    # i is some kind of index
    Y <- X[i]
    oY <- order(Y$x)
    Ysort <- Y[oY]
    idYsort <- ((1:nX)[i])[oY]
  }
  nY <- npoints(Y)
  force(nr)
  # call C
  zz <- .C("locpcfx",
           nn1 = as.integer(nY),
           x1  = as.double(Ysort$x),
           y1  = as.double(Ysort$y),
           id1 = as.integer(idYsort),
           nn2 = as.integer(nX),
           x2  = as.double(Xsort$x),
           y2  = as.double(Xsort$y),
           id2 = as.integer(idXsort),
           nnr = as.integer(nr),
           rmaxi=as.double(rmax),
           del=as.double(delta),
           pcf=as.double(double(nr * nY)),
           PACKAGE="spatstat")
  out <- matrix(zz$pcf, nr, nY)
  if(!missi) {
    # reorder columns to match original
    out[, oY] <- out
  }
  out <- out/(2 * pi * lambda)
  attr(out, "r") <- seq(0, rmax, length=nr)
  attr(out, "delta") <- delta
  class(out) <- c("localpcfmatrix", class(out))
  return(out)
}

print.localpcfmatrix <- function(x, ...) {
  cat("Matrix of local pair correlation estimates\n")
  nc <- ncol(x)
  nr <- nrow(x)
  cat(paste("pcf estimates for", nc, ngettext(nc, "point", "points"), "\n"))
  rval <- attr(x, "r")
  cat(paste("r values from 0 to", max(rval), "in", nrow(x), "steps\n"))
  return(invisible(NULL))
}

plot.localpcfmatrix <- function(x, ...) {
  xname <- deparse(substitute(x))
  rval <- attr(x, "r")
  do.call("matplot",
          resolve.defaults(list(rval, x),
                           list(...),
                           list(type="l", main=xname,
                                xlab="r", ylab="pair correlation")))
}

"[.localpcfmatrix" <-
  function(x, i, ...) {
    r     <- attr(x, "r")
    delta <- attr(x, "delta")
    class(x) <- "matrix"
    if(missing(i)) {
      x <- x[ , ...]
    } else {
      x <- x[i, ...]
      if(is.matrix(i))
        return(x)
      r <- r[i]
    }
    if(!is.matrix(x))
      x <- matrix(x, nrow=length(r))
    attr(x, "r") <- r
    attr(x, "delta") <- delta
    class(x) <- c("localpcfmatrix", class(x))
    return(x)
}

