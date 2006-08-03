nncross <- function(X, Y) {
  verifyclass(X, "ppp")
  verifyclass(Y, "ppp")
  if(X$n == 0)
    return(data.frame(dist=numeric(0), which=integer(0)))
  if(Y$n == 0)
    return(data.frame(dist=rep(Inf, X$n),
                      which=rep(NA, X$n)))
  oX <- order(X$y)
  X <- X[oX]
  oY <- order(Y$y)
  Y <- Y[oY]
  nndv <- numeric(X$n)
  nnwh <- integer(X$n)
  z <- .C("nnXwhich",
     n1=as.integer(X$n),
     x1=as.double(X$x),
     y1=as.double(X$y),
     n2=as.integer(Y$n),
     x2=as.double(Y$x),
     y2=as.double(Y$y),
     nnd=as.double(nndv),
     nnwhich=as.integer(nnwh),
     huge=as.double(diameter(X$window)),
     PACKAGE="spatstat")
  nndv[oX] <- z$nnd
  nnwcode <- z$nnwhich + 1
  nnwcode[nnwcode < 1] <- NA
  nnwh[oX] <- oY[nnwcode]
  return(data.frame(dist=nndv, which=nnwh))
}

