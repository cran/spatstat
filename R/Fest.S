#
#	Fest.S
#
#	S function empty.space()
#	Computes estimates of the empty space function
#
#	$Revision: 4.5 $	$Date: 2001/11/06 06:59:55 $
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
# First discretise
	dwin <- as.mask(X$window, eps)
        dX <- ppp(X$x, X$y, window=dwin)
#        
# histogram breakpoints 
#
        breaks <- handle.r.b.args(r, breaks, dwin, eps)
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
		b <- bdist.pixels(dX$window, coords=F)
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
        hh <- hist(dist,breaks=breaks$val,plot=F)$counts
        edf <- cumsum(hh)/sum(hh)
        result$raw <- edf

# append theoretical value for Poisson
        lambda <- X$n/area.owin(X$window)
        result$theo <- 1 - exp( - lambda * pi * result$r^2)

# neaten up and return        
        result$breaks <- NULL
        
	return(data.frame(result))
}

	
