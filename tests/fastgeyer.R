# checks validity of fast C implementation of Geyer interaction
require(spatstat)
data(redwood)
X <- redwood
Q <- quadscheme(X)
U <- union.quad(Q)
EP <- equalpairs.quad(Q)
gg <- Geyer(0.11, 2)
# The value r=0.11 is chosen to avoid hardware numerical effects (gcc bug 323).
# It avoids being close any value of pairdist(redwood).
# The nearest such values are 0.1077.. and 0.1131..
# By contrast if r = 0.1 there are values differing from 0.1 by 3e-17
a <- pairsat.family$eval(X,U,EP,gg$pot,gg$par,"border")
b <-         gg$fasteval(X,U,EP,gg$pot,gg$par,"border")
if(!all(a==b))
  stop("Results of Geyer()$fasteval and pairsat.family$eval do not match")
