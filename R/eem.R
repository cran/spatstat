# eem.R
#
# Computes the Stoyan-Grabarnik "exponential energy weights" 
#
# $Revision: 1.2 $ $Date: 2006/06/02 08:39:24 $
#

eem <- function(fit, check=TRUE) {
  verifyclass(fit, "ppm")
  lambda <- fitted.ppm(fit, check=check)
  Q <- quad.ppm(fit)
  Z <- is.data(Q)
  eemarks <- 1/lambda[Z]
  attr(eemarks, "type") <- "eem"
  attr(eemarks, "typename") <- "exponential energy marks"
  return(eemarks)
}
