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

#
#  Test mathlegend code
#

data(cells)
K <- Kest(cells)
plot(K)
plot(K, . ~ r)
plot(K, . - theo ~ r)
plot(K, sqrt(./pi)  ~ r)
plot(K, cbind(iso, theo) ~ r)
plot(K, cbind(iso, theo) - theo ~ r)
plot(K, sqrt(cbind(iso, theo)/pi)  ~ r)
plot(K, cbind(iso/2, -theo) ~ r)
plot(K, cbind(iso/2, trans/2) - theo ~ r)

# test expansion of .x and .y
plot(K, . ~ .x)
plot(K, . - theo ~ .x)
plot(K, .y - theo ~ .x)
plot(K, sqrt(.y) - sqrt(theo) ~ .x)

# problems with parsing weird strings in levels(marks(X))
# noted by Ulf Mehlig

levels(marks(amacrine)) <- c("Nastricreechia krorluppia", "Homo habilis")
plot(Kcross(amacrine))
plot(alltypes(amacrine, "K"))
plot(alltypes(amacrine, "J"))
plot(alltypes(amacrine, pcfcross))


     
