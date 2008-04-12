#
# areadiff.R
#
#  $Revision: 1.5 $  $Date: 2008/04/02 13:37:04 $
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
  DUP <- spatstat.options("dupC")
  z <- .C("areadiff",
          ux = as.double(u[1]),
          uy = as.double(u[2]),
          rad = as.double(r),
          x   = as.double(X$x),
          y   = as.double(X$y),
          nn  = as.integer(X$n),
          ngrid = as.integer(ngrid),
          answer = as.double(numeric(1)),
          DUP=DUP,
          PACKAGE="spatstat")
  z$answer
}
