#
#
#      distmap.R
#
#      $Revision: 1.2 $     $Date: 2004/10/01 03:08:02 $
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
  U <- exactdt(X, ...)
  xc <- U$xcol
  yr <- U$yrow
  V <- im(U$d, xc, yr)
  attr(V, "index") <- im(U$i, xc, yr)
  attr(V, "bdry")  <- im(U$b, xc, yr)
  return(V)
}

distmap.owin <- function(X, ...) {
  verifyclass(X, "owin")
  X <- as.mask(X, ...)
  xc <- X$xcol
  yr <- X$yrow
  U <- exactPdt(X)
  V <- im(U$d, xc, yr)
  attr(V, "bdry") <- im(U$b, xc, yr)
  return(V)
}

