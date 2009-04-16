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
#	$Revision: 4.26 $	$Date: 2009/04/14 23:27:43 $
#
################################################################################

"Gcross" <-		
function(X, i, j, r=NULL, breaks=NULL, ...)
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
  stopifnot(is.multitype(X))
#
  marx <- marks(X, dfok=FALSE)
  if(missing(i)) i <- levels(marx)[1]
  if(missing(j)) j <- levels(marx)[2]
#  
  I <- (marx == i)
  if(sum(I) == 0) stop("No points are of type i")
        
  if(i == j)
    result <- Gest(X[I], r=r, breaks=breaks, ...)
  else {
    J <- (marx == j)
    if(sum(J) == 0) stop("No points are of type j")
    result <- Gmulti(X, I, J, r=r, breaks=breaks, disjoint=FALSE, ...)
  }
  result <-
    rebadge.fv(result,
               substitute(Gcross[i,j](r), list(i=paste(i),j=paste(j))),
               "G*",
               new.yexp=substitute(Gcross[list(i,j)](r),
                                   list(i=paste(i),j=paste(j))))
  return(result)
}	

"Gdot" <- 	
function(X, i, r=NULL, breaks=NULL, ...) {
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
  stopifnot(is.multitype(X))
#
  marx <- marks(X, dfok=FALSE)
  if(missing(i)) i <- levels(marx)[1]
  I <- (marx == i)
  if(sum(I) == 0) stop("No points are of type i")
  J <- rep(TRUE, X$n)	# i.e. all points
# 
  result <- Gmulti(X, I, J, r, breaks, disjoint=FALSE, ...)
  result <- rebadge.fv(result,
                       substitute(Gdot[i](r), list(i=paste(i))),
                       "Gm")
  return(result)
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
        labl <- c("r", "%spois(r)", "%sbord(r)", "%skm(r)",
                  "lambda(r)", "%sraw(r)")
        desc <- c("distance argument r",
                  "theoretical Poisson %s",
                  "border corrected estimate of %s",
                  "Kaplan-Meier estimate of %s",
                  "Kaplan-Meier estimate of hazard function lambda(r)",
                  "uncorrected estimate of %s")
        Z <- fv(Z, "r", substitute(Gmulti(r), NULL),
                "km", . ~ r, alim, labl, desc, fname="Gm")
        attr(Z, "dotnames") <- c("km", "rs", "theo")
        unitname(Z) <- unitname(X)
	return(Z)
}	


