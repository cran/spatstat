#
#   reach.R
#
#  $Revision: 1.1 $   $Date: 2005/03/12 01:32:31 $
#

reach <- function(x, ...) {
  UseMethod("reach")
}

reach.interact <- function(x, ...) {
  verifyclass(x, "interact")
  irange <- x$irange
  if(is.null(irange))
    return(Inf)
  if(!is.function(irange))
    stop("Internal error - x$irange is not a function")
  ir <- irange(x)
  if(is.na(ir))
    ir <- Inf
  return(ir)
}

reach.ppm <- function(x, ..., epsilon=0) {
  verifyclass(x, "ppm")
  coeffs <- coef(x)
  inte <- x$interaction
  irange <- inte$irange
  if(is.null(irange))
    return(NA)
  if(!is.function(irange))
    stop("Internal error - x$interaction$irange is not a function")
  ir <- irange(inte, coeffs=coeffs, epsilon=epsilon)
  if(is.na(ir))
    ir <- Inf
  return(ir)
}

  
  
