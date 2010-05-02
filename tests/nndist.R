# nndist.R
# Check that nndist and nnwhich give
# results consistent with direct calculation from pairdist

require(spatstat)
eps <- sqrt(.Machine$double.eps)
f <- function(mat,k) { apply(mat, 1, function(z,n) { sort(z)[n]  }, n=k+1) }
g <- function(mat,k) { apply(mat, 1, function(z,n) { order(z)[n] }, n=k+1) }

# Two dimensions

X <- runifpoint(42)

nn <- nndist(X)
nnP <- f(pairdist(X), 1)
if(any(abs(nn - nnP) > eps))
  stop("nndist.ppp does not agree with pairdist")

nn5 <- nndist(X, k=5)
nn5P <- f(pairdist(X), 5)
if(any(abs(nn5 - nn5P) > eps))
  stop("nndist.ppp(k=5) does not agree with pairdist")

nw <- nnwhich(X)
nwP <- g(pairdist(X), 1)
if(any(nw != nwP))
  stop("nnwhich.ppp does not agree with pairdist")

nw5 <- nnwhich(X, k=5)
nw5P <- g(pairdist(X), 5)
if(any(nw5 != nw5P))
  stop("nnwhich.ppp(k=5) does not agree with pairdist")

# Three dimensions

X <- runifpoint3(42)

nn <- nndist(X)
nnP <- f(pairdist(X), 1)
if(any(abs(nn - nnP) > eps))
  stop("nndist.pp3 does not agree with pairdist")

nn5 <- nndist(X, k=5)
nn5P <- f(pairdist(X), 5)
if(any(abs(nn5 - nn5P) > eps))
  stop("nndist.pp3(k=5) does not agree with pairdist")

nw <- nnwhich(X)
nwP <- g(pairdist(X), 1)
if(any(nw != nwP))
  stop("nnwhich.pp3 does not agree with pairdist")

nw5 <- nnwhich(X, k=5)
nw5P <- g(pairdist(X), 5)
if(any(nw5 != nw5P))
  stop("nnwhich.pp3(k=5) does not agree with pairdist")

# m dimensions

X <- runifpointx(42, boxx(c(0,1),c(0,1),c(0,1),c(0,1)))

nn <- nndist(X)
nnP <- f(pairdist(X), 1)
if(any(abs(nn - nnP) > eps))
  stop("nndist.ppx does not agree with pairdist")

nn5 <- nndist(X, k=5)
nn5P <- f(pairdist(X), 5)
if(any(abs(nn5 - nn5P) > eps))
  stop("nndist.ppx(k=5) does not agree with pairdist")

nw <- nnwhich(X)
nwP <- g(pairdist(X), 1)
if(any(nw != nwP))
  stop("nnwhich.ppx does not agree with pairdist")

nw5 <- nnwhich(X, k=5)
nw5P <- g(pairdist(X), 5)
if(any(nw5 != nw5P))
  stop("nnwhich.ppx(k=5) does not agree with pairdist")
