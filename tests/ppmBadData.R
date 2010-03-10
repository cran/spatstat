# ppmBadData.R
# $Revision: 1.2 $ $Date: 2008/11/16 08:16:33 $

require(spatstat)

SEED <- 42

K <- 101
A <- 500
X <- seq(0, A, length=K)
G <- expand.grid(x=X, y=X)
FOO <- function(x,y) { sin(x)^2 + cos(y)^2 }
M1 <- im(matrix(FOO(G$x, G$y), K, K), xcol=X, yrow=X)
M <- im(matrix(FOO(G$x, G$y), K, K))
BAR <- function(x) { exp(-6.618913 + 5.855337 * x - 8.432483 * x^2) }
V <- im(BAR(M$v), xcol=X, yrow=X)
# V <- eval.im(exp(-6.618913 + 5.855337 * M - 8.432483 * M^2))
set.seed(SEED)
Y <- rpoispp(V)
fY <- ppm(Y, ~cv + I(cv^2), covariates=list(cv=M), correction="translate")
diagnose.ppm(fY)

