#
#	Gest.S
#
#	Compute estimates of nearest neighbour distance distribution function G
#
#	$Revision: 4.3 $	$Date: 2001/11/06 07:00:19 $
#
################################################################################
#
"Gest" <-
"nearest.neighbour" <-
function(X, r=NULL, breaks=NULL, ...) {
#	X		point pattern (of class ppp)
#				(unless specified by X$window)
#       r:              (optional) values of argument r  
#	breaks:		(optional) breakpoints for argument r
#
	verifyclass(X, "ppp")
        
#  compute nearest neighbour distances
	nnd <- nndist(X$x, X$y)
		
#  determine breakpoints for r values
        breaks <- handle.r.b.args(r, breaks, X$window)

#  UNCORRECTED e.d.f. of nearest neighbour distances: use with care
        hh <- hist(nnd,breaks=breaks$val,plot=F)$counts
        edf <- cumsum(hh)/sum(hh)
	
#  distance to boundary
        bdry <- bdist.points(X)

#  observations
	o <- pmin(nnd,bdry)
#  censoring indicators
	d <- (nnd <= bdry)
#
# calculate Kaplan-Meier and border correction (Reduced Sample) estimators
	result <- km.rs(o, bdry, d, breaks)
#        
# append uncorrected e.d.f.        
        result$raw <- edf
# append theoretical value for Poisson
        lambda <- X$n/area.owin(X$window)
        result$theo <- 1 - exp( - lambda * pi * result$r^2)

# neaten up and return        
        result$breaks <- NULL

	return(data.frame(result))
}	

nndist <- function(x, y) {
	#  computes the vector of nearest-neighbour distances 
	#  for the pattern of points (x[i],y[i])
	#	
#  matrix of squared distances between all pairs of points
	sq <- function(a, b) { (a-b)^2 }
	squd <-  outer(x, x, sq) + outer(y, y, sq)
#  reset diagonal to a large value so it is excluded from minimum
	diag(squd) <- Inf
#  nearest neighbour distances
	nnd <- sqrt(apply(squd,1,min))
	nnd
}
