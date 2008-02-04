#
#  distances.psp.R
#
#  Hausdorff distance and Euclidean separation for psp objects
#
#  $Revision: 1.5 $ $Date: 2007/10/24 09:41:15 $
#
#

pairdist.psp <- function(X, ..., method="Fortran", type="Hausdorff") {
  verifyclass(X, "psp")
  if(X$n == 0)
    return(matrix(, 0, 0))
  type <- pickoption("type", type,
                     c(Hausdorff="Hausdorff",
                       hausdorff="Hausdorff",
                       separation="separation"))

  # Compute min dist from each endpoint of X to each segment of Y
  D12 <- AsymmHausdorff.psp(X, X, method=method)
  
  switch(type,
         Hausdorff={
           # maximum is Hausdorff metric
           D <- array(pmax(D12, t(D12)), dim=dim(D12))
         },
         separation={
           # first take minimum of endpoint-to-segment distances
           D <- array(pmin(D12, t(D12)), dim=dim(D12))
           # Identify pairs of segments which cross
           cross <- test.selfcrossing.psp(X)
           # Assign separation = 0 to such pairs
           D[cross] <- 0
         })
  return(D)
}

crossdist.psp <- function(X, Y, ..., method="Fortran", type="Hausdorff") {
  verifyclass(X, "psp")
  Y <- as.psp(Y)
  if(X$n * Y$n == 0)
    return(matrix(, X$n, Y$n))

  type <- pickoption("type", type,
                     c(Hausdorff="Hausdorff",
                       hausdorff="Hausdorff",
                       separation="separation"))

  
  # Compute min dist from each endpoint of X to each segment of Y
  DXY <- AsymmHausdorff.psp(X, Y, method=method)
  # Compute min dist from each endpoint of Y to each segment of X
  DYX <- AsymmHausdorff.psp(Y, X, method=method)
  
  switch(type,
         Hausdorff={
           # maximum is Hausdorff metric
           D <- array(pmax(DXY, t(DYX)), dim=dim(DXY))
         },
         separation={
           # first take minimum of endpoint-to-segment distances
           D <- array(pmin(DXY, t(DYX)), dim=dim(DXY))
           # Identify pairs of segments which cross
           cross <- test.crossing.psp(X, Y)
           # Assign separation = 0 to such pairs
           D[cross] <- 0
         })
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

  

