#
#	Gest.S
#
#	Compute estimates of nearest neighbour distance distribution function G
#
#	$Revision: 4.17 $	$Date: 2009/04/07 06:46:50 $
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
        W <- X$window
        npoints <-X$n
# Intensity estimate  
        lambda <- npoints/area.owin(W)

#  determine breakpoints for r values
        rmaxdefault <- rmax.rule("G", W, lambda)
        breaks <- handle.r.b.args(r, breaks, W, rmaxdefault=rmaxdefault)

# handle special case of empty pattern
        if(X$n == 0) {
          rvalues <- breaks$r
          nr <- length(rvalues)
          zero <- rep(0, nr)
          result <- list(rs     = zero,
                         km     = zero,
                         hazard = zero,
                         r      = rvalues,
                         raw    = zero,
                         theo   = zero)
          unitname(result) <- unitname(X)
          return(result)
        }
        
#  compute nearest neighbour distances
	nnd <- nndist(X$x, X$y)
		
#  UNCORRECTED e.d.f. of nearest neighbour distances: use with care
        rightmost <- breaks$max
        hh <- hist(nnd[nnd <= rightmost],breaks=breaks$val,plot=FALSE)$counts
        edf <- cumsum(hh)/length(nnd)
	
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
        result$theo <- 1 - exp( - lambda * pi * result$r^2)

# neaten up and return        
        result$breaks <- NULL

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
        Z <- fv(Z, "r", substitute(G(r), NULL), "km", . ~ r, alim, labl, desc,
                fname="G")
        attr(Z, "dotnames") <- c("km", "rs", "theo")
        unitname(Z) <- unitname(X)
	return(Z)
}	

