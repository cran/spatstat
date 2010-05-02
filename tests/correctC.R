# tests for agreement between C and interpreted code
# for interpoint distances

require(spatstat)
eps <- .Machine$double.eps * 4

# pairdist.ppp
X <- rpoispp(42)
dC <- pairdist(X, method="C")
dR <- pairdist(X, method="interpreted")
if(any(abs(dC - dR) > eps))
  stop("Algorithms for pairdist() do not agree")

dC <- pairdist(X, periodic=TRUE, method="C")
dR <- pairdist(X, periodic=TRUE, method="interpreted")
if(any(abs(dC - dR) > eps))
  stop("Algorithms for pairdist(periodic=TRUE) do not agree")

# crossdist.ppp
Y <- rpoispp(42)
dC <- crossdist(X, Y, method="C")
dR <- crossdist(X, Y, method="interpreted")
if(any(abs(dC - dR) > eps))
  stop("Algorithms for crossdist() do not agree")

dC <- crossdist(X, Y, periodic=TRUE, method="C")
dR <- crossdist(X, Y, periodic=TRUE, method="interpreted")
if(any(abs(dC - dR) > eps))
  stop("Algorithms for crossdist(periodic=TRUE) do not agree")

# nndist.ppp
nnC <- nndist(X, method="C")
nnI <- nndist(X, method="interpreted")
if(any(abs(nnC - nnI) > eps))
  stop("Algorithms for nndist() do not agree")

nn3C <- nndist(X, k=3, method="C")
nn3I <- nndist(X, k=3, method="interpreted")
if(any(abs(nn3C - nn3I) > eps))
  stop("Algorithms for nndist(k=3) do not agree")

# nnwhich.ppp
nwC <- nnwhich(X, method="C")
nwI <- nnwhich(X, method="interpreted")
if(any(nwC != nwI))
  stop("Algorithms for nnwhich() do not agree")

nw3C <- nnwhich(X, k=3, method="C")
nw3I <- nnwhich(X, k=3, method="interpreted")
if(any(nw3C != nw3I))
  stop("Algorithms for nnwhich(k=3) do not agree")







