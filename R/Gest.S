#
#	Gest.S
#
#	Compute estimates of nearest neighbour distance distribution function G
#
#	$Revision: 4.5 $	$Date: 2002/05/27 11:26:54 $
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
        hh <- hist(nnd,breaks=breaks$val,plot=FALSE)$counts
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

