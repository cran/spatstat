#
#
#   rescale.R
#
#   $Revision: 1.3 $ $Date: 2007/10/24 09:41:15 $
#
#

rescale <- function(X, s) {
  UseMethod("rescale")
}

rescale.ppp <- function(X, s) {
  Y <- affine.ppp(X, mat=diag(c(1/s,1/s)))
  unitname(Y) <- rescale(unitname(X), s)
  return(Y)
}

rescale.owin <- function(X, s) {
  Y <- affine.owin(X, mat=diag(c(1/s,1/s)))
  unitname(Y) <- rescale(unitname(X), s)
  return(Y)
}

rescale.psp <- function(X, s) {
  Y <- affine.psp(X, mat=diag(c(1/s,1/s)))
  unitname(Y) <- rescale(unitname(X), s)
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


