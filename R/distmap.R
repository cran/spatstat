#
#
#      distmap.R
#
#      $Revision: 1.6 $     $Date: 2006/11/20 04:41:33 $
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
  uni <- units(W)
  dmat <- e$d
  imat <- e$i
  V <- im(dmat, W$xcol, W$yrow, units=uni)
  I <- im(imat, W$xcol, W$yrow, units=uni)
  if(X$window$type == "rectangle") {
    # distance to frame boundary
    bmat <- e$b
    B <- im(bmat, W$xcol, W$yrow, units=uni)
  } else {
    # distance to window boundary, not frame boundary
    bmat <- bdist.pixels(W, coords=FALSE)
    B <- im(bmat, W$xcol, W$yrow, units=uni)
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
  uni <- units(X)
  xc <- X$xcol
  yr <- X$yrow
  U <- exactPdt(X)
  V <- im(U$d, xc, yr, units=uni)
  B <- im(U$b, xc, yr, units=uni)
  attr(V, "bdry") <- B
  return(V)
}

