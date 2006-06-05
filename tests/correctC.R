# tests for agreement between C and interpreted code
require(spatstat)

# pairdist.ppp
X <- rpoispp(42)
dC <- pairdist(X, method="C")
dR <- pairdist(X, method="interpreted")
range(dC - dR)

dC <- pairdist(X, periodic=TRUE, method="C")
dR <- pairdist(X, periodic=TRUE, method="interpreted")
range(dC - dR)

# crossdist.ppp
Y <- rpoispp(42)
dC <- crossdist(X, Y, method="C")
dR <- crossdist(X, Y, method="interpreted")
range(dC - dR)

dC <- crossdist(X, Y, periodic=TRUE, method="C")
dR <- crossdist(X, Y, periodic=TRUE, method="interpreted")
range(dC - dR)



