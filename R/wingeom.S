#
#	wingeom.S	Various geometrical computations in windows
#
#
#	$Revision: 4.15 $	$Date: 2006/09/26 07:44:09 $
#
#
#
#
#-------------------------------------
area.owin <- function(w) {
	verifyclass(w, "owin")
        switch(w$type,
               rectangle = {
		width <- abs(diff(w$xrange))
		height <- abs(diff(w$yrange))
		area <- width * height
               },
               polygonal = {
                 area <- sum(unlist(lapply(w$bdry, area.xypolygon)))
               },
               mask = {
                 pixelarea <- abs(w$xstep * w$ystep)
                 npixels <- sum(w$m)
                 area <- pixelarea * npixels
               },
               stop("Unrecognised window type")
        )
        return(area)
}

eroded.areas <- function(w, r) {
	verifyclass(w, "owin")
	
	switch(w$type,
               rectangle = {
                 width <- abs(diff(w$xrange))
                 height <- abs(diff(w$yrange))
                 areas <- pmax(width - 2 * r, 0) * pmax(height - 2 * r, 0)
               },
               polygonal = {
                 # warning("Approximating polygonal window by digital image")
                 w <- as.mask(w)
                 areas <- eroded.areas(w, r)
               },
               mask = {
                 # distances from each pixel to window boundary
                 b <- bdist.pixels(w, coords=FALSE)
                 # histogram breaks to satisfy hist()
                 Bmax <- max(b, r)
                 breaks <- c(-1,r,Bmax+1)
                 # histogram of boundary distances
                 h <- hist(b, breaks=breaks, plot=FALSE, probability=FALSE)$counts
                 # reverse cumulative histogram
                 H <- revcumsum(h)
                 # drop first entry corresponding to r=-1
                 H <- H[-1]
                 # convert count to area
                 pixarea <- w$xstep * w$ystep
                 areas <- pixarea * H
               },
 	       stop("unrecognised window type")
               )
	areas
}	

diameter <- function(w) {
	verifyclass(w, "owin")
	
        width <- abs(diff(w$xrange))
        height <- abs(diff(w$yrange))
        
        sqrt(width^2 + height^2)
}

even.breaks.owin <- function(w) {
	verifyclass(w, "owin")
        Rmax <- diameter(w)
        make.even.breaks(Rmax, Rmax/(100 * sqrt(2)))
}

unit.square <- function() { owin(c(0,1),c(0,1)) }

square <- function(r=1) {
  if(length(r) != 1 || !is.numeric(r))
    stop("argument r must be a single number")
  if(is.na(r) || !is.finite(r))
    stop("argument r is NA or infinite")
  if(r <= 0)
    stop("side of square must be positive")
  owin(c(0,r),c(0,r))
}

overlap.owin <- function(A, B) {
  # compute the area of overlap between two windows
  At <- A$type
  Bt <- B$type
  if(At=="rectangle" && Bt=="rectangle") {
    xmin <- max(A$xrange[1],B$xrange[1])
    xmax <- min(A$xrange[2],B$xrange[2])
    if(xmax <= xmin) return(0)
    ymin <- max(A$yrange[1],B$yrange[1])
    ymax <- min(A$yrange[2],B$yrange[2])
    if(ymax <= ymin) return(0)
    return((xmax-xmin) * (ymax-ymin))
  }
  if((At=="rectangle" && Bt=="polygonal")
     || (At=="polygonal" && Bt=="rectangle")
     || (At=="polygonal" && Bt=="polygonal"))
  {
    AA <- as.polygonal(A)$bdry
    BB <- as.polygonal(B)$bdry
    area <- 0
    for(i in seq(AA))
      for(j in seq(BB))
        area <- area + overlap.xypolygon(AA[[i]], BB[[j]])
    return(area)
  }
  if(At=="mask") {
    # count pixels in A that belong to B
    pixelarea <- abs(A$xstep * A$ystep)
    x <- as.vector(raster.x(A)[A$m])
    y <- as.vector(raster.y(A)[A$m])
    ok <- inside.owin(x, y, B) 
    return(pixelarea * sum(ok))
  }
  if(Bt== "mask") {
    # count pixels in B that belong to A
    pixelarea <- abs(B$xstep * B$ystep)
    x <- as.vector(raster.x(B)[B$m])
    y <- as.vector(raster.y(B)[B$m])
    ok <- inside.owin(x, y, A)
    return(pixelarea * sum(ok))
  }
  stop("Internal error")
}

#
#
#  Intersection and union of windows
#
#
intersect.owin <- function(A, B) {
  verifyclass(A, "owin")
  verifyclass(B, "owin")

  # chicken out 
  if(A$type == "polygonal")
    A <- as.mask(A)
  if(B$type == "polygonal")
    B <- as.mask(B)
  
  # determine intersection of x and y ranges
  intersect.ranges <- function(a, b) {
    lo <- max(a[1],b[1])
    hi <- min(a[2],b[2])
    if(lo >= hi) stop("Intersection is empty")
    return(c(lo, hi))
  }
  xr <- intersect.ranges(A$xrange, B$xrange)
  yr <- intersect.ranges(A$yrange, B$yrange)
  C <- owin(xr, yr)
  
  Arect <- (A$type == "rectangle")
  Brect <- (B$type == "rectangle")

  if(Arect && Brect)
    return(owin(xr, yr))

  # Otherwise, we'll need this
  trim.mask <- function(M, R) {
    # M is a mask,
    # R is a rectangle inside bounding.box(M)
    # Extract subset of image grid
    within.range <- function(u, v) { (u >= v[1]) & (u <= v[2]) }
    yrowok <- within.range(M$yrow, R$yrange)
    xcolok <- within.range(M$xcol, R$xrange)
    if(sum(yrowok) == 0 || sum(xcolok) == 0)
      stop("result is empty")
    Z <- M
    Z$xrange <- R$xrange
    Z$yrange <- R$yrange
    Z$yrow <- M$yrow[yrowok]
    Z$xcol <- M$xcol[xcolok]
    Z$m <- M$m[yrowok, xcolok]
    Z$dim <- dim(Z$m)
    return(Z)
  }

  if(!Arect && Brect)
    return(trim.mask(A, C))
  else if(Arect && !Brect) 
    return(trim.mask(B, C))
  else {
    # both are masks
    # First trim A 
    D <- trim.mask(A, C)
    # Then determine which pixels of D are inside B
    x <- raster.x(D)
    y <- raster.y(D)
    ok <- inside.owin(x, y, B)
    # 'and'
    if(!all(ok))
      D$m[!ok] <- FALSE
    return(D)
  }
          
  stop("Internal error")
}


union.owin <- function(A, B) {
  verifyclass(A, "owin")
  verifyclass(B, "owin")

  if(A$type == "rectangle" && B$type == "rectangle") {
    if(is.subset.owin(A, B))
      return(B)
    else if (is.subset.owin(B,A))
      return(A)
  }

  C <- owin(range(A$xrange, B$xrange),
            range(A$yrange, B$yrange))
  
  C <- as.mask(C)
  x <- raster.x(C)
  y <- raster.y(C)
  ok <- inside.owin(x, y, A) | inside.owin(x, y, B)

  if(!any(ok))
    stop("Internal error: union is empty")

  if(!all(ok))
    C$m[!ok] <- FALSE

  return(C)
}

expand.owin <- function(W, f=1) {
  # expand bounding box of 'win'
  # by factor 'f' in **area**
  if(f <= 0)
    stop("f must be > 0")
  if(f == 1)
    return(W)
  bb <- bounding.box(W)
  xr <- bb$xrange
  yr <- bb$yrange
  fff <- (sqrt(f) - 1)/2
  Wexp <- owin(xr + fff * c(-1,1) * diff(xr),
               yr + fff * c(-1,1) * diff(yr))
  return(Wexp)
}

trim.rectangle <- function(W, xmargin=0, ymargin=xmargin) {
  if(W$type != "rectangle")
    stop("Internal error: tried to trim margin off non-rectangular window")
  if(length(xmargin) == 1)
    xmargin <- rep(xmargin, 2)
  if(length(ymargin) == 1)
    ymargin <- rep(ymargin, 2)
  if(any(xmargin < 0) || any(ymargin < 0))
    stop("values of xmargin, ymargin must be nonnegative")
  if(sum(xmargin) > diff(W$xrange))
    stop("window is too small to cut off margins of the width specified")
  if(sum(ymargin) > diff(W$yrange))
    stop("window is too small to cut off margins of the height specified")
  owin(W$xrange + c(1,-1) * xmargin,
       W$yrange + c(1,-1) * ymargin)
}


