#
#
#      distances.R
#
#      $Revision: 1.1 $     $Date: 2002/05/27 11:24:43 $
#
#
#      Interpoint distances
#
#

"pairdist"<-
function(x, y=NULL, method="C")
{
  # extract x and y coordinate vectors
  if(verifyclass(x, "ppp", fatal=FALSE)) 
    xy <- list(x=x$x, y=x$y)
  else 
    xy <- xy.coords(x,y)[c("x","y")]
  x <- xy$x
  y <- xy$y

  n <- length(x)
  if(length(y) != n)
    stop("lengths of x and y do not match")

  # special cases
  if(n == 0)
    return(numeric(0))
  else if(n == 1)
    return(matrix(1,nrow=1,ncol=1))
  
  switch(method,
         interpreted={
           xx <- matrix(rep(x, n), nrow = n)
           yy <- matrix(rep(y, n), nrow = n)
           d <- sqrt((xx - t(xx))^2 + (yy - t(yy))^2)
         },
         C={
           d <- numeric( n * n)
           z<- .C("pairdist", n = as.integer(n),
                  x= as.double(x), y= as.double(y), d= as.double(d))
           d <- matrix(z$d, nrow=n, ncol=n)
         },
         stop(paste("Unrecognised method \"", method, "\"", sep=""))
       )
  invisible(d)
}

"nndist"<-
function(x, y=NULL, method="C")
{
	#  computes the vector of nearest-neighbour distances 
	#  for the pattern of points (x[i],y[i])
	#
  # extract x and y coordinate vectors
  if(verifyclass(x, "ppp", fatal=FALSE)) 
    xy <- list(x=x$x, y=x$y)
  else 
    xy <- xy.coords(x,y)[c("x","y")]
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
    return(NA)
  
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
                  as.double(big))
           nnd[o] <- z$nnd
         },
         stop(paste("Unrecognised method \"", method, "\"", sep=""))
         )
  invisible(nnd)
}


