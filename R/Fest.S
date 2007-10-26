#
#	Fest.S
#
#	S function empty.space()
#	Computes estimates of the empty space function
#
#	$Revision: 4.18 $	$Date: 2007/10/24 09:41:15 $
#
"Fest" <- 	
"empty.space" <-
function(X, eps = NULL, r=NULL, breaks=NULL) {
#
#	pp:		point pattern (an object of class 'ppp')
#	eps:		raster grid mesh size for distance transform
#				(unless specified by pp$window)
#       r:              (optional) values of argument r  
#	breaks:		(optional) breakpoints for argument r
#
# Intensity estimate  
        lambda <- X$n/area.owin(X$window)
# First discretise
	dwin <- as.mask(X$window, eps)
        dX <- ppp(X$x, X$y, window=dwin, check=FALSE)
#        
# histogram breakpoints 
#
        rmaxdefault <- rmax.rule("F", dwin, lambda)
        breaks <- handle.r.b.args(r, breaks, dwin, eps,
                                  rmaxdefault=rmaxdefault)
#
#  compute distances and censoring distances
	if(X$window$type == "rectangle") {
                # original data were in a rectangle
                # output of exactdt() is sufficient
		e <- exactdt(dX)
		dist <- e$d
		bdry <- e$b
	} else {
                # window is irregular..
          
                # Distance transform & boundary distance for all pixels
		e <- exactdt(dX)
		b <- bdist.pixels(dX$window, coords=FALSE)
                # select only those pixels inside mask
		dist <- e$d[dwin$m]
		bdry <- b[dwin$m]
	} 
# censoring indicators
	d <- (dist <= bdry)
#  observed distances
	o <- pmin(dist, bdry)
#        
#
# calculate Kaplan-Meier and border corrected estimates
	result <- km.rs(o, bdry, d, breaks)

# also calculate UNCORRECTED e.d.f. !!!! use with care
        rightmost <- breaks$max
        hh <- hist(dist[dist <= rightmost],breaks=breaks$val,plot=FALSE)$counts
        edf <- cumsum(hh)/length(dist)
        result$raw <- edf

# append theoretical value for Poisson
        result$theo <- 1 - exp( - lambda * pi * result$r^2)

# neaten up and return        
        result$breaks <- NULL

# convert to class "fv"
        result <- as.data.frame(result)
        Z <- result[, c("r", "theo", "rs", "km", "hazard", "raw")]
        alim <- range(result$r[result$km <= 0.9])
        labl <- c("r", "Fpois(r)", "Fbord(r)", "Fkm(r)",
                  "lambda(r)", "Fraw(r)")
        desc <- c("distance argument r",
                  "theoretical Poisson F(r)",
                  "border corrected estimate of F(r)",
                  "Kaplan-Meier estimate of F(r)",
                  "Kaplan-Meier estimate of hazard function lambda(r)",
                  "uncorrected estimate of F(r)")
        Z <- fv(Z, "r", substitute(F(r), NULL), "km", . ~ r, alim, labl, desc)
        attr(Z, "dotnames") <- c("km", "rs", "theo")
        unitname(Z) <- unitname(X)
	return(Z)
}

	
