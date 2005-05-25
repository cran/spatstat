# eem.R
#
# Computes the Stoyan-Grabarnik "exponential energy weights" 
#
# $Revision: 1.1 $ $Date: 2005/05/11 19:49:02 $
#

eem <- function(fit) {
  verifyclass(fit, "ppm")
  lambda <- fitted.ppm(fit)
  Q <- quad.ppm(fit)
  Z <- is.data(Q)
  eemarks <- 1/lambda[Z]
  attr(eemarks, "type") <- "eem"
  attr(eemarks, "typename") <- "exponential energy marks"
  return(eemarks)
}
