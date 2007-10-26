#	Gmulti.S
#
#	Compute estimates of nearest neighbour distance distribution functions
#	for multitype point patterns
#
#	S functions:	
#		Gcross                G_{ij}
#		Gdot		      G_{i\bullet}
#		Gmulti	              (generic)
#
#	$Revision: 4.20 $	$Date: 2007/10/24 09:41:15 $
#
################################################################################

"Gcross" <-		
function(X, i=1, j=2, r=NULL, breaks=NULL, ...)
{
#	computes G_{ij} estimates
#
#	X		marked point pattern (of class 'ppp')
#	i,j		the two mark values to be compared
#  
#       r:              (optional) values of argument r  
#	breaks:		(optional) breakpoints for argument r
#
	X <- as.ppp(X)
	if(!is.marked(X))
		stop(paste("point pattern has no", sQuote("marks")))
	window <- X$window
#
        marx <- marks(X, dfok=FALSE)
	I <- (marx == i)
	if(sum(I) == 0) stop("No points are of type i")
	J <- (marx == j)
	if(sum(J) == 0) stop("No points are of type j")
#
	result <- Gmulti(X, I, J, r, breaks, disjoint=(i!=j))
        attr(result, "ylab") <- substitute(Gcross(r), NULL)
        result
}	

"Gdot" <- 	
function(X, i=1, r=NULL, breaks=NULL, ...) {
#  Computes estimate of 
#      G_{i\bullet}(t) = 
#  P( a further point of pattern in B(0,t)| a type i point at 0 )
#	
#	X		marked point pattern (of class ppp)
#  
#       r:              (optional) values of argument r  
#	breaks:		(optional) breakpoints for argument r
#
	X <- as.ppp(X)

	if(!is.marked(X))
		stop(paste("point pattern has no", sQuote("marks")))
#
	I <- (marks(X) == i)
	if(sum(I) == 0) stop("No points are of type i")
	J <- rep(TRUE, X$n)	# i.e. all points
# 
	result <- Gmulti(X, I, J, r, breaks, disjoint=FALSE)
        attr(result, "ylab") <- substitute(Gdot(r), NULL)
        result
}	

	
##########

"Gmulti" <- 	
function(X, I, J, r=NULL, breaks=NULL, ..., disjoint=NULL) {
#
#  engine for computing the estimate of G_{ij} or G_{i\bullet}
#  depending on selection of I, J
#  
#	X		marked point pattern (of class ppp)
#	
#	I,J		logical vectors of length equal to the number of points
#			and identifying the two subsets of points to be
#			compared.
#  
#       r:              (optional) values of argument r  
#	breaks:		(optional) breakpoints for argument r
#
	verifyclass(X, "ppp")
	W <- X$window
        npoints <- X$n
        area <- area.owin(W)
# check I and J
	if(!is.logical(I) || !is.logical(J))
		stop("I and J must be logical vectors")
	if(length(I) != X$n || length(J) != X$n)
		stop("length of I or J does not equal the number of points")
        nI <- sum(I)
        nJ <- sum(J)
	if(nI == 0) stop("No points satisfy condition I")
	if(nJ == 0) stop("No points satisfy condition J")
        if(is.null(disjoint))
          disjoint <- !any(I & J)
#  determine breakpoints for r values
        lamJ <- nJ/area
        rmaxdefault <- rmax.rule("G", W, lamJ)
        breaks <- handle.r.b.args(r, breaks, W, rmaxdefault=rmaxdefault)
        brks <- breaks$val
        rmax <- breaks$max
#  "type I to type J" nearest neighbour distances
        XI <- X[I]
        XJ <- X[J]
        if(disjoint) 
          nnd <- nncross(XI, XJ)$dist
        else {
          seqnp <- seq(npoints)
          iX <- seqnp[I]
          iY <- seqnp[J]
          nnd <- nncross(XI, XJ, iX, iY)$dist
        }
#  distance to boundary from each type i point
        bdry <- bdist.points(XI)
#  observations
	o <- pmin(nnd,bdry)
#  censoring indicators
	d <- (nnd <= bdry)
#
# calculate
	result <- km.rs(o, bdry, d, brks)
        result$breaks <- NULL

#  UNCORRECTED e.d.f. of I-to-J nearest neighbour distances: use with care
        rightmost <- breaks$max
        hh <- hist(nnd[nnd <= rightmost],breaks=brks,plot=FALSE)$counts
        result$raw <- cumsum(hh)/length(nnd)

# theoretical value for marked Poisson processes
        result$theo <- 1 - exp( - lamJ * pi * result$r^2)
        
# convert to class "fv"
        result <- as.data.frame(result)
        Z <- result[, c("r", "theo", "rs", "km", "hazard", "raw")]
        alim <- range(result$r[result$km <= 0.9])
        labl <- c("r", "Gpois(r)", "Gbord(r)", "Gkm(r)",
                  "lambda(r)", "Graw(r)")
        desc <- c("distance argument r",
                  "theoretical Poisson G(r)",
                  "border corrected estimate of G(r)",
                  "Kaplan-Meier estimate of G(r)",
                  "Kaplan-Meier estimate of hazard function lambda(r)",
                  "uncorrected estimate of G(r)")
        Z <- fv(Z, "r", substitute(Gmulti(r), NULL), "km", . ~ r, alim, labl, desc)
        attr(Z, "dotnames") <- c("km", "rs", "theo")
        unitname(Z) <- unitname(X)
	return(Z)
}	


