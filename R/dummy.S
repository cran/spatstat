#
#	dummy.S
#
#	Utilities for generating patterns of dummy points
#
#       $Revision: 4.8 $     $Date: 2003/07/22 19:13:44 $
#
#	corners()	corners of window
#	gridcenters()	points of a rectangular grid
#	stratrand()	random points in each tile of a rectangular grid
#	spokes()	Rolf's 'spokes' arrangement
#	
#	concatxy()	concatenate any lists of x, y coordinates
#
#	default.dummy()	Default action to create a dummy pattern
#		
	
corners <- function(window) {
	window <- as.owin(window)
	x <- window$xrange[c(1,2,1,2)]
	y <- window$yrange[c(1,1,2,2)]
	return(list(x=x, y=y))
}

gridcenters <-	
gridcentres <- function(window, nx, ny) {
	window <- as.owin(window)
	xr <- window$xrange
	yr <- window$yrange
	x <- seq(xr[1], xr[2], length = 2 * nx + 1)[2 * (1:nx)]
	y <- seq(yr[1], yr[2], length = 2 * ny + 1)[2 * (1:ny)]
	x <- rep(x, ny)
	y <- rep(y, rep(nx, ny))
	return(list(x=x, y=y))
}

stratrand <- function(window,nx,ny, k=1) {
	
	# divide window into an nx * ny grid of tiles
	# and place k points at random in each tile
	
	window <- as.owin(window)

	wide  <- diff(window$xrange)/(nx - 1)
	high  <- diff(window$yrange)/(ny - 1)
        cent <- gridcentres(window, nx, ny)
	cx <- rep(cent$x, k)
	cy <- rep(cent$y, k)
	n <- nx * ny * k
	x <- cx + runif(n, min = -wide/2, max = wide/2)
	y <- cy + runif(n, min = -high/2, max = high/2)
	return(list(x=x,y=y))
}

tilecentroids <- function (W, nx, ny)
{
  W <- as.owin(W)
  if(W$type == "rectangle")
    return(gridcentres(W, nx, ny))
  else {
    # approximate
    W   <- as.mask(W)
    xx  <- as.vector(raster.x(W)[W$m])
    yy  <- as.vector(raster.y(W)[W$m])
    pid <- gridindex(xx,yy,W$xrange,W$yrange,nx,nx)$index
    x   <- tapply(xx,pid,mean)
    y   <- tapply(yy,pid,mean)
    return(list(x=x,y=y))
  }
}

tilemiddles <- function (W, nx, ny)
{
    if(W$type == "rectangle")
      return(gridcentres(W, nx, ny))
    
    # pixel approximation to window
    # This matches the pixel approximation used to compute tile areas
    # and ensures that dummy points are generated only inside those tiles
    # that have nonzero digital area
    M   <- as.mask(W)
    xx  <- as.vector(raster.x(M))[M$m]
    yy  <- as.vector(raster.y(M))[M$m]
    # label the pixels by tile index
    pid <- gridindex(xx,yy,W$xrange,W$yrange,nx,ny)$index
    pidvalues <- levels(as.factor(pid))
    # Latter is for comparison with the output of tapply( ..., pid, ...)
    
    ######## 1st try: centroids of tiles
    # (usually inside tile, e.g. if it's convex, but not always inside)
    cx <- tapply(xx, pid, mean)
    cy <- tapply(yy, pid, mean)
    # centroids which are inside window
    cok <- inside.owin(cx, cy, W)
    x <- cx[cok]
    y <- cy[cok]
    # tiles for which this strategy did not work
    nbg <- pidvalues[!cok]
    
    ######### 2nd try: middle point in list of pixels in each tile
    # (always inside tile, by construction)
    todo <- pid %in% nbg
    xx <- xx[todo]
    yy <- yy[todo]
    pid <- pid[todo]
    middle <- function(v) { n <- length(v); mid <- ceiling(n/2); v[mid]}
    mx   <- tapply(xx,pid,middle)
    my   <- tapply(yy,pid,middle)
    
    x <- c( x, mx)
    y <- c( y, my)
    return(list(x=x,y=y))
}

spokes <- function(x, y, nrad = 3, nper = 3, fctr = 1.5, Mdefault=1) {
	#
	# Rolf Turner's "spokes" arrangement
	#
	# Places dummy points on radii of circles 
	# emanating from each data point x[i], y[i]
	#
	#       nrad:    number of radii from each data point
	#       nper:	 number of dummy points per radius
	#       fctr:	 length of largest radius = fctr * M
	#                where M is mean nearest neighbour distance in data
	#
	M <- if(length(x) > 1) mean(nndist(x,y)) else Mdefault
	lrad  <- fctr * M / nper
	theta <- 2 * pi * (1:nrad)/nrad
	cs    <- cos(theta)
	sn    <- sin(theta)
	xt    <- lrad * as.vector((1:nper) %o% cs)
	yt    <- lrad * as.vector((1:nper) %o% sn)
	xd    <- as.vector(outer(x, xt, "+"))
	yd    <- as.vector(outer(y, yt, "+"))
	
	return(list(x=xd,y=yd))
}
	
# concatenate any number of list(x,y) into a list(x,y)
		
concatxy <- function(...) {
	x <- unlist(lapply(list(...), function(w) {w$x}))
	y <- unlist(lapply(list(...), function(w) {w$y}))
	if(length(x) != length(y))
		stop("Internal error: lengths of x and y unequal")
	return(list(x=x,y=y))
}

#------------------------------------------------------------

default.dummy <- function(X, nx=NULL, ny=NULL, random=FALSE) {
	# default action to create dummy points.
	# regular grid of (order) nx * ny tiles
	# plus corner points of window frame.
	# 
	X <- as.ppp(X)
	win <- X$window
        # default dimensions
        if(is.null(nx)) nx <- default.ngrid(X)
        if(is.null(ny)) ny <- nx
        # make dummy points
        dummy <- if(random) 
                  stratrand(win, nx, ny, 1)
                else
                  tilemiddles(win, nx, ny)
        # add corner points
	dummy <- concatxy(dummy, corners(win))
	# restrict to window (if it's a mask)
	ok <- inside.owin(dummy$x, dummy$y, win)
	x <- dummy$x[ok]
	y <- dummy$y[ok]
	ppp(x, y, window=win)
}


default.ngrid <- function(X) {
	# default dimensions of rectangular grid of dummy points
        # for data pattern X
  X <- as.ppp(X)
  max(30, 10 * ceiling(2 * sqrt(X$n)/10))
}

