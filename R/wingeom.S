#
#	wingeom.S	Various geometrical computations in windows
#
#
#	$Revision: 4.7 $	$Date: 2002/05/13 12:41:10 $
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
                 H <- rev(cumsum(rev(h)))
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
  


  
