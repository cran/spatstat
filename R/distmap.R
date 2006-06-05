#
#
#      distmap.R
#
#      $Revision: 1.4 $     $Date: 2006/06/01 02:24:45 $
#
#
#     Distance transform
#
#
distmap <- function(X, ...) {
  UseMethod("distmap")
}

distmap.ppp <- function(X, ...) {
  verifyclass(X, "ppp")
  e <- exactdt(X, ...)
  W <- e$w
  dmat <- e$d
  imat <- e$i
  V <- im(dmat, W$xcol, W$yrow)
  I <- im(imat, W$xcol, W$yrow)
  if(X$window$type == "rectangle") {
    # distance to frame boundary
    bmat <- e$b
    B <- im(bmat, W$xcol, W$yrow)
  } else {
    # distance to window boundary, not frame boundary
    bmat <- bdist.pixels(W, coords=FALSE)
    B <- im(bmat, W$xcol, W$yrow)
    # clip all to window
    V <- V[W, drop=FALSE]
    I <- I[W, drop=FALSE]
    B <- B[W, drop=FALSE]
  }
  attr(V, "index") <- I
  attr(V, "bdry")  <- B
  return(V)
}

distmap.owin <- function(X, ...) {
  verifyclass(X, "owin")
  X <- as.mask(X, ...)
  xc <- X$xcol
  yr <- X$yrow
  U <- exactPdt(X)
  V <- im(U$d, xc, yr)
  B <- im(U$b, xc, yr)
  attr(V, "bdry") <- B
  return(V)
}

