#
#  distances.psp.R
#
#  Hausdorff distances for psp objects
#
#  $Revision: 1.5 $ $Date: 2007/10/24 09:41:15 $
#
#

pairdist.psp <- function(X, ..., method="Fortran") {
  verifyclass(X, "psp")
  if(X$n == 0)
    return(matrix(, 0, 0))
  D12 <- AsymmHausdorff.psp(X, X, method=method)
  D <- array(pmax(D12, t(D12)), dim=dim(D12))
  return(D)
}

crossdist.psp <- function(X, Y, ..., method="Fortran") {
  verifyclass(X, "psp")
  Y <- as.psp(Y)
  if(X$n * Y$n == 0)
    return(matrix(, X$n, Y$n))

  DXY <- AsymmHausdorff.psp(X, Y, method=method)
  DYX <- AsymmHausdorff.psp(Y, X, method=method)
  D <- array(pmax(DXY, t(DYX)), dim=dim(DXY))
  return(D)
}

nndist.psp <- function(X, ..., k=1, method="Fortran") {
  verifyclass(X, "psp")
  n <- X$n
  if(n == 0)
    return(numeric(0))
  else if(n <= k)
    return(rep(Inf, n))
  D <- pairdist.psp(X, ..., method=method)
  diag(D) <- Inf
  if(k == 1) 
    NND <- apply(D, 1, min)
  else
    NND <- apply(D, 1, function(z,k) { sort(z)[k] }, k=k)
  return(NND)
}

AsymmHausdorff.psp <- function(X, Y, method="Fortran") {
  if(method != "Fortran" && method != "interpreted")
    stop(paste("Unrecognised method", sQuote(method)))
  # Extract endpoints of X
  EX <- endpoints.psp(X, "both")
  idX <- attr(EX, "id")
  # compute min dist from each endpoint of X to each segment of Y
  DPL <- distppll(cbind(EX$x,EX$y), Y$ends, mintype=0, method=method)
  # maximise over each pair of endpoints
  Dist <- as.vector(DPL)
  Point <- as.vector(idX[row(DPL)])
  Segment <- as.vector(col(DPL))
  DXY <- tapply(Dist, list(factor(Point), factor(Segment)), max)
  return(DXY)
}

  

