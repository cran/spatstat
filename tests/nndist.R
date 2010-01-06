
require(spatstat)
eps <- .Machine$double.eps * 4

X <- runifpoint(42)

nnC <- nndist(X, method="C")
nnI <- nndist(X, method="interpreted")
if(any(abs(nnC - nnI) > eps))
  stop("Algorithms for nndist() do not agree")

nn3C <- nndist(X, k=3, method="C")
nn3I <- nndist(X, k=3, method="interpreted")
if(any(abs(nn3C - nn3I) > eps))
  stop("Algorithms for nndist(k=3) do not agree")

nwC <- nnwhich(X, method="C")
nwI <- nnwhich(X, method="interpreted")
if(any(nwC != nwI))
  stop("Algorithms for nnwhich() do not agree")

nw3C <- nnwhich(X, k=3, method="C")
nw3I <- nnwhich(X, k=3, method="interpreted")
if(any(nw3C != nw3I))
  stop("Algorithms for nnwhich(k=3) do not agree")







