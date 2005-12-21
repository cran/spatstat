#
#  distances.psp.R
#
#  Hausdorff distances for psp objects
#
#  $Revision: 1.3 $ $Date: 2005/12/21 01:14:49 $
#
#

pairdist.psp <- function(X, ..., method="Fortran") {
  verifyclass(X, "psp")

  D12 <- AsymmHausdorff.psp(X, X, method=method)
  D <- array(pmax(D12, t(D12)), dim=dim(D12))
  return(D)
}

crossdist.psp <- function(X, Y, ..., method="Fortran") {
  verifyclass(X, "psp")
  Y <- as.psp(Y)

  DXY <- AsymmHausdorff.psp(X, Y, method=method)
  DYX <- AsymmHausdorff.psp(Y, X, method=method)
  D <- array(pmax(DXY, t(DYX)), dim=dim(DXY))
  return(D)
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

  

