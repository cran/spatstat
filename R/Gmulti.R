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
#	$Revision: 4.32 $	$Date: 2010/03/08 08:23:04 $
#
################################################################################

"Gcross" <-		
function(X, i, j, r=NULL, breaks=NULL, ..., correction=c("rs", "km", "han"))
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
  if(!is.marked(X, dfok=FALSE))
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
    result <- Gmulti(X, I, J, r=r, breaks=breaks, disjoint=FALSE, ...,
                     correction=correction)
  }
  result <-
    rebadge.fv(result,
               substitute(Gcross[i,j](r), list(i=paste(i),j=paste(j))),
               "Gcross",
               new.yexp=substitute(Gcross[list(i,j)](r),
                                   list(i=paste(i),j=paste(j))))
  return(result)
}	

"Gdot" <- 	
function(X, i, r=NULL, breaks=NULL, ..., correction=c("km","rs","han")) {
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
  result <- Gmulti(X, I, J, r, breaks, disjoint=FALSE, ...,
                   correction=correction)
  result <- rebadge.fv(result,
                       substitute(Gdot[i](r), list(i=paste(i))),
                       "Gdot")
  return(result)
}	

	
##########

"Gmulti" <- 	
function(X, I, J, r=NULL, breaks=NULL, ..., disjoint=NULL,
         correction=c("rs", "km", "han")) {
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
# choose correction(s)
  correction.given <- !missing(correction) && !is.null(correction)
  if(is.null(correction))
    correction <- c("rs", "km", "han")
  correction <- pickoption("correction", correction,
                           c(none="none",
                             border="rs",
                             rs="rs",
                             KM="km",
                             km="km",
                             Kaplan="km",
                             han="han",
                             Hanisch="han",
                             best="km"),
                           multi=TRUE)
#  determine breakpoints for r values
  lamJ <- nJ/area
  rmaxdefault <- rmax.rule("G", W, lamJ)
  breaks <- handle.r.b.args(r, breaks, W, rmaxdefault=rmaxdefault)
  brks <- breaks$val
  rmax <- breaks$max
  rvals <- breaks$r
  zeroes <- rep(0, length(rvals))
# initialise fv object
  df <- data.frame(r=rvals, theo=1-exp(-lamJ * pi * rvals^2))
  Z <- fv(df, "r", substitute(Gmulti(r), NULL), "theo", . ~ r,
          c(0,rmax),
          c("r", "%s[pois](r)"), 
          c("distance argument r", "theoretical Poisson %s"),
          fname="Gmulti")
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
# calculate estimates
  
  if("none" %in% correction) {
    #  UNCORRECTED e.d.f. of nearest neighbour distances: use with care
    if(npoints == 0)
      edf <- zeroes
    else {
      hh <- hist(nnd[nnd <= rmax],breaks=breaks$val,plot=FALSE)$counts
      edf <- cumsum(hh)/length(nnd)
    }
    Z <- bind.fv(Z, data.frame(raw=edf), "%s[raw](r)",
                 "uncorrected estimate of %s", "raw")
  }

  if("han" %in% correction) {
    # Hanisch style estimator
    if(npoints == 0)
      G <- zeroes
    else {
      #  uncensored distances
      x <- nnd[d]
      #  weights
      a <- eroded.areas(W, rvals)
      # calculate Hanisch estimator
      h <- hist(x[x <= rmax], breaks=breaks$val, plot=FALSE)$counts
      G <- cumsum(h/a)
      G <- G/max(G[is.finite(G)])
    }
    # add to fv object
    Z <- bind.fv(Z, data.frame(han=G),
                 "%s[han](r)", 
                 "Hanisch estimate of %s",
                 "han")
    # modify recommended plot range
    attr(Z, "alim") <- range(rvals[G <= 0.9])
  }
  
  if(any(correction %in% c("rs", "km"))) {
    # calculate Kaplan-Meier and border correction (Reduced Sample) estimators
    if(npoints == 0)
      result <- data.frame(rs=zeroes, km=zeroes, hazard=zeroes)
    else {
      result <- km.rs(o, bdry, d, breaks)
      result <- as.data.frame(result[c("rs", "km", "hazard")])
    }
    # add to fv object
    Z <- bind.fv(Z, result,
                 c("%s[bord](r)", "%s[km](r)", "hazard(r)"),
                 c("border corrected estimate of %s",
                   "Kaplan-Meier estimate of %s",
                   "Kaplan-Meier estimate of hazard function lambda(r)"),
                 "km")
    # modify recommended plot range
    attr(Z, "alim") <- range(rvals[result$km <= 0.9])
  }
  nama <- names(Z)
  fvnames(Z, ".") <- rev(nama[!(nama %in% c("r", "hazard"))])
  unitname(Z) <- unitname(X)
  return(Z)
}	


