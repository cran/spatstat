#
#  random.R
#
#  Random point pattern generators for a linear network
#
#  $Revision: 1.1 $   $Date: 2010/06/09 05:11:52 $
#

rpoislpp <- function(lambda, L, ...) {
  verifyclass(L, "linnet")
  LL <- L$lines
  X <- rpoisppOnLines(lambda, LL, ...)
  Y <- lpp(X, L)
  return(Y)
}

runiflpp <- function(n, L) {
  verifyclass(L, "linnet")
  LL <- L$lines
  X <- runifpointOnLines(n, LL)
  Y <- lpp(X, L)
  return(Y)
}
