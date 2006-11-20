#
#
#   rescale.R
#
#   $Revision: 1.2 $ $Date: 2006/11/20 05:00:49 $
#
#

rescale <- function(X, s) {
  UseMethod("rescale")
}

rescale.ppp <- function(X, s) {
  Y <- affine.ppp(X, mat=diag(c(1/s,1/s)))
  units(Y) <- rescale(units(X), s)
  return(Y)
}

rescale.owin <- function(X, s) {
  Y <- affine.owin(X, mat=diag(c(1/s,1/s)))
  units(Y) <- rescale(units(X), s)
  return(Y)
}

rescale.psp <- function(X, s) {
  Y <- affine.psp(X, mat=diag(c(1/s,1/s)))
  units(Y) <- rescale(units(X), s)
  return(Y)
}
  
rescale.units <- function(X, s) {
  if(summary(X)$vanilla)
    return(X)
  if(!is.numeric(s) || length(s) != 1 || s <= 0)
    stop("s should be a positive number")
  X$multiplier <- s * X$multiplier
  return(X)
}


