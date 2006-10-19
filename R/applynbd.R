# 	applynbd.R
#
#     $Revision: 1.10 $     $Date: 2006/10/10 04:22:48 $
#
#  applynbd()
# For each point, identify either
#	 - all points within distance R
#        - the closest N points  
#        - those points satisfying some constraint
# and apply the function FUN to them
#
#  markstat()
#      simple application of applynbd
#################################################################


applynbd <- function(X, FUN, N, R, criterion, exclude=FALSE, ...) {

     nopt <- (!missing(N)) + (!missing(R)) + (!missing(criterion))
     if(nopt > 1)
       stop(paste("exactly one of the arguments",
                  paste(sQuote(c("N", "R", "criterion")), collapse=", "),
                  "must be given"))
     else if(nopt == 0)
       stop(paste("must specify one of the arguments",
                  sQuote("N"), ",", sQuote("R"), "or", sQuote("criterion")))
     
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
	fun(Y=Z[which], current=c(x=Z$x[here], y=Z$y[here]),
            dists=distances[which], dranks=dranks[which], ...) 
     }

     result <- apply(a, 1, go, Z=X, fun=FUN, ...)
  
     return(result)
}


markstat <- function(X, fun, N, R, ...) {
  verifyclass(X, "ppp")
  stopifnot(is.function(fun))
  if(!missing(R) && !missing(N))
    stop("Do not specify both R and N")
  if(missing(R) && missing(N))
    stop("either R or N should be given")
  statfun <- function(Y, current, dists, dranks, func, ...)
    { func(marks(Y), ...) }
  if(!missing(R))
    applynbd(X, statfun, R=R, func=fun, ...)
  else
    applynbd(X, statfun, N=N, func=fun, ...)
}
