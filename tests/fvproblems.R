# fvproblems.R

require(spatstat)

# This appears in the workshop notes
# Problem detected by Martin Bratschi

Jdif <- function(X, ..., i) {
  Jidot <- Jdot(X, ..., i=i)
  J <- Jest(X, ...)
  dif <- eval.fv(Jidot - J)
  return(dif)
}
data(amacrine)

Z <- Jdif(amacrine, i="on")
