#
#
#      distances.R
#
#      $Revision: 1.11 $     $Date: 2005/07/28 04:45:33 $
#
#
#      Interpoint distances
#
#

pairdist <- function(X, ..., method="C") {
  UseMethod("pairdist")
}

pairdist.ppp <- function(X, ..., method="C") {
  verifyclass(X, "ppp")
  return(pairdist.default(X$x, X$y, method="C"))
}

pairdist.default <-
  function(X, Y=NULL, ..., method="C")
{
  xy <- xy.coords(X,Y)[c("x","y")]
  x <- xy$x
  y <- xy$y

  n <- length(x)
  if(length(y) != n)
    stop("lengths of x and y do not match")

  # special cases
  if(n == 0)
    return(matrix(numeric(0), nrow=0, ncol=0))
  else if(n == 1)
    return(matrix(0,nrow=1,ncol=1))
  
  switch(method,
         interpreted={
           xx <- matrix(rep(x, n), nrow = n)
           yy <- matrix(rep(y, n), nrow = n)
           d <- sqrt((xx - t(xx))^2 + (yy - t(yy))^2)
         },
         C={
           d <- numeric( n * n)
           z<- .C("pairdist", n = as.integer(n),
                  x= as.double(x), y= as.double(y), d= as.double(d),
                  PACKAGE="spatstat")
           d <- matrix(z$d, nrow=n, ncol=n)
         },
         stop(paste("Unrecognised method \"", method, "\"", sep=""))
       )
  invisible(d)
}

nndist <- function(X, ..., method="C") {
  UseMethod("nndist")
}

nndist.ppp <- function(X, ..., method="C") {
  verifyclass(X, "ppp")
  return(nndist.default(X$x, X$y, method="C"))
}

nndist.default <-
  function(X, Y=NULL, ..., method="C")
{
	#  computes the vector of nearest-neighbour distances 
	#  for the pattern of points (x[i],y[i])
	#
  xy <- xy.coords(X,Y)[c("x","y")]
  x <- xy$x
  y <- xy$y

  # validate
  n <- length(x)
  if(length(y) != n)
    stop("lengths of x and y do not match")

  # special cases
  if(n == 0)
    return(numeric(0))
  else if(n == 1)
    return(Inf)
  
  switch(method,
         interpreted={
           #  matrix of squared distances between all pairs of points
           sq <- function(a, b) { (a-b)^2 }
           squd <-  outer(x, x, sq) + outer(y, y, sq)
           #  reset diagonal to a large value so it is excluded from minimum
           diag(squd) <- Inf
           #  nearest neighbour distances
           nnd <- sqrt(apply(squd,1,min))
         },
         C={
           n <- length(x)
           nnd<-numeric(n)
           o <- order(y)
           big <- sqrt(.Machine$double.xmax)
           z<- .C("nndistsort",
                  n= as.integer(n),
                  x= as.double(x[o]), y= as.double(y[o]), nnd= as.double(nnd),
                  as.double(big),
                  PACKAGE="spatstat")
           nnd[o] <- z$nnd
         },
         stop(paste("Unrecognised method \"", method, "\"", sep=""))
         )
  invisible(nnd)
}

crossdist <- function(X, Y, ..., method="C") {
  UseMethod("crossdist")
}

crossdist.ppp <- function(X, Y, ..., method="C") {
  verifyclass(X, "ppp")
  return(crossdist.default(X$x, X$y, Y$x, Y$y, method=method))
}

crossdist.default <-
  function(X, Y, x2, y2, ..., method="C")
{
  x1 <- X
  y1 <- Y
  # returns matrix[i,j] = distance from (x1[i],y1[i]) to (x2[j],y2[j])
  if(length(x1) != length(y1))
    stop("lengths of x and y do not match")
  if(length(x2) != length(y2))
    stop("lengths of x2 and y2 do not match")
  n1 <- length(x1)
  n2 <- length(x2)
  if(n1 == 0 || n2 == 0)
    return(matrix(numeric(0), nrow=n1, ncol=n2))
  switch(method,
         interpreted = {
                 X1 <- matrix(rep(x1, n2), ncol = n2)
                 Y1 <- matrix(rep(y1, n2), ncol = n2)
                 X2 <- matrix(rep(x2, n1), ncol = n1)
                 Y2 <- matrix(rep(y2, n1), ncol = n1)
                 d <- sqrt((X1 - t(X2))^2 + (Y1 - t(Y2))^2)
                 return(d)
               },
               C = {
                 z<- .C("crossdist",
                        nfrom = as.integer(n1),
                        xfrom = as.double(x1),
                        yfrom = as.double(y1),
                        nto = as.integer(n2),
                        xto = as.double(x2),
                        yto = as.double(y2),
                        d = as.double(matrix(0, nrow=n1, ncol=n2)),
                        PACKAGE="spatstat")
                 return(matrix(z$d, nrow=n1, ncol=n2))
               },
               stop(paste("Unrecognised method", method))
               )
      }

