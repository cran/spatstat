# 	applynbd.R
#
#     $Revision: 1.1 $     $Date: 2002/07/18 10:53:11 $
#
#
# For each point, identify either
#	 - all points within distance R
#        - the closest N points  
#        - those points satisfying some constraint
# and apply the function FUN to them
#
#################################################################


applynbd <- function(X, FUN, N, R, criterion, exclude=FALSE, ...) {

     nopt <- (!missing(N)) + (!missing(R)) + (!missing(criterion))
     if(nopt > 1)
       stop("exactly one of the arguments \"N\", \"R\", \"criterion\" must be given")
     else if(nopt == 0)
       stop("must specify one of the arguments \"N\", \"R\" or \"criterion\"")
     
     X <- as.ppp(X)
     npts <- X$n

     # compute matrix of pairwise distances
     dist <- pairdist(X$x,X$y)	

     # compute row ranks (avoid ties)
     rankit <- function(x) {  u <- numeric(length(x)); u[order(x)] <- seq(x); return(u) }
     drank <- t(apply(dist, 1, rankit)) - 1

     if(!missing(R)) {
	     # select points closer than R
	     included <- (dist <= R)
     } else if(!missing(N)) {
	     # select N closest points
	     if(N < 1)
		stop("Value of N must be at least 1")
	     if(exclude)
		included <- (drank <= N) 
	     else
		included <- (drank <= N-1)
     } else {
            # some funny criterion
	    included <- matrix(, nrow=npts, ncol=npts)
	    for(i in 1:npts) 
		included[i,] <- criterion(dist[i,], drank[i,])
     }
     
    if(exclude) 
	diag(included) <- FALSE

     # bind into an array
     a <- array(c(included, dist, drank, row(included)), dim=c(npts,npts,4))

     # what to do with a[i, , ]
     go <- function(ai, Z, fun, ...) { 
	which <- as.logical(ai[,1])
        distances <- ai[,2]
	dranks <- ai[,3]
        here <- ai[1,4]	
	fun(Z[which], c(x=Z$x[here], y=Z$y[here]), distances[which], dranks[which], ...) 
     }

     result <- apply(a, 1, go, Z=X, fun=FUN, ...)
  
     return(result)
}

