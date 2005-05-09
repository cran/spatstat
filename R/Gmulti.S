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
#	$Revision: 4.7 $	$Date: 2005/05/05 01:18:21 $
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
		stop("point pattern has no \'marks\'")
	window <- X$window
# 
	I <- (X$marks == i)
	if(sum(I) == 0) stop("No points are of type i")
	J <- (X$marks == j)
	if(sum(J) == 0) stop("No points are of type j")
#
	Gmulti(X, I, J, r, breaks)
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
		stop("point pattern has no \'marks\'")
#
	I <- (X$marks == i)
	if(sum(I) == 0) stop("No points are of type i")
	J <- rep(TRUE, X$n)	# i.e. all points
# 
	Gmulti(X, I, J, r, breaks)
}	

	
##########

"Gmulti" <- 	
function(X, I, J, r=NULL, breaks=NULL, ...) {
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
	if(!is.marked(X))
		stop("point pattern has no \'marks\'")
	window <- X$window
        brks <- handle.r.b.args(r, breaks, window)$val        
# Extract points of group I and J
	if(!is.logical(I) || !is.logical(J))
		stop("I and J must be logical vectors")
	if(length(I) != X$n || length(J) != X$n)
		stop("length of I or J does not equal the number of points")
	if(sum(I) == 0) stop("No points satisfy condition I")
	xI <- X$x[I]
	yI <- X$y[I]
	if(sum(J) == 0) stop("No points satisfy condition J")
	xJ <- X$x[J]
	yJ <- X$y[J]
#  compute squared distance from each type i point to
#  the nearest other point of any type
	xIm <- matrix(xI, nrow=length(xI), ncol=length(xJ))
	yIm <- matrix(yI, nrow=length(yI), ncol=length(yJ))
	xJm <- matrix(xJ, nrow=length(xJ), ncol=length(xI))
	yJm <- matrix(yJ, nrow=length(yJ), ncol=length(yI))
	squd <-  (xIm - t(xJm))^2 + (yIm - t(yJm))^2
# identify which pairs in 'squd' correspond to identical points
	dag <- (diag(1, X$n) != 0)
	dag <- dag[I,J]
#  reset this 'diagonal' to a large value
	oo <- (diff(window$xrange) + diff(window$yrange) )^2
	squd[dag] <- oo
#  "type I to type J" nearest neighbour distances
	nnd <- sqrt(apply(squd,1,min))
#  distance to boundary from each type i point
        bdry <- bdist.points(X[I, ])
#  observations
	o <- pmin(nnd,bdry)
#  censoring indicators
	d <- (nnd <= bdry)
#
# calculate
	result <- km.rs(o, bdry, d, brks)
        result$breaks <- NULL

#  UNCORRECTED e.d.f. of I-to-J nearest neighbour distances: use with care
        hh <- hist(nnd,breaks=brks,plot=FALSE)$counts
        result$raw <- cumsum(hh)/sum(hh)

# theoretical value for marked Poisson processes
        lamJ <- sum(J)/area.owin(window)
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
        Z <- fv(Z, "r", "G(r)", "km", cbind(km, rs, theo) ~ r, alim, labl, desc)
	return(Z)
}	


