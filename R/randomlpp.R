#
#  random.R
#
#  Random point pattern generators for a linear network
#
#  $Revision: 1.2 $   $Date: 2012/06/06 09:55:15 $
#

rpoislpp <- function(lambda, L, ...) {
  verifyclass(L, "linnet")
  X <- datagen.rpoisppOnLines(lambda, as.psp(L), ...)
  Y <- lpp(X, L)
  return(Y)
}

runiflpp <- function(n, L) {
  if(!is.numeric(n) || (length(n) != 1) || (n < 0) || (n %% 1 != 0))
    stop("n should be a single nonnegative integer")
  verifyclass(L, "linnet")
  X <- datagen.runifpointOnLines(n, as.psp(L))
  Y <- lpp(X, L)
  return(Y)
}
