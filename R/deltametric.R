#
#   deltametric.R
#
#   Delta metric
#
#   $Revision: 1.1 $  $Date: 2009/08/28 05:44:58 $
#

deltametric <- function(A, B, p=2, c=Inf, ...) {
  stopifnot(is.numeric(p) && length(p) == 1 && p > 0)
  dA <- distmap(A, ...)
  dB <- distmap(B, ...)
  if(!is.infinite(c)) {
    dA <- eval.im(pmin(dA, c))
    dB <- eval.im(pmin(dB, c))
  }
  if(is.infinite(p)) {
    # L^infinity
    Z <- eval.im(abs(dA-dB))
    delta <- summary(Z)$max
  } else {
    # L^p
    Z <- eval.im(abs(dA-dB)^p)
    iZ <- summary(Z)$mean
    delta <- iZ^(1/p)
  }
  return(delta)
}





