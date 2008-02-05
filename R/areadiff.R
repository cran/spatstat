#
# areadiff.R
#
#  $Revision: 1.3 $  $Date: 2008/02/01 19:56:36 $
#
# Computes sufficient statistic for area-interaction process
#
# Invokes areadiff.c
#

areadiff <- function(u, X, r, ngrid=256) {
  verifyclass(X, "ppp")
  stopifnot(is.numeric(u) && length(u) == 2)
  stopifnot(is.numeric(r) && length(r) == 1 && r > 0)
  # find points of X within distance 2r of location u
  close <- ((X$x - u[1])^2 + (X$y - u[2])^2 < 4 * r^2)
  X <- X[close]
  # invoke C routine
  z <- .C("areadiff",
          ux = as.double(u[1]),
          uy = as.double(u[2]),
          rad = as.double(r),
          x   = as.double(X$x),
          y   = as.double(X$y),
          nn  = as.integer(X$n),
          ngrid = as.integer(ngrid),
          answer = as.double(numeric(1)),
          PACKAGE="spatstat")
  z$answer
}
