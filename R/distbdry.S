#
#	distbdry.S		Distance to boundary
#
#	$Revision: 4.2 $	$Date: 2001/08/07 09:22:52 $
#
# -------- functions ----------------------------------------
#
#	bdist.points()
#                       compute vector of distances 
#			from each point of point pattern
#                       to boundary of window
#
#       bdist.pixels()
#                       compute matrix of distances from each pixel
#                       to boundary of window
#
#       erode.mask()    erode the window mask by a distance r
#                       [yields a new window]
#
#       erode.owin()    erode the window by a distance r
#                       [yields a new window]
#
#
# 
"bdist.points"<-
function(X)
{
	verifyclass(X, "ppp") 
	x <- X$x
	y <- X$y
	window <- X$window
        switch(window$type,
               rectangle = {
		xmin <- min(window$xrange)
		xmax <- max(window$xrange)
		ymin <- min(window$yrange)
		ymax <- max(window$yrange)
		result <- pmin(x - xmin, xmax - x, y - ymin, ymax - y)
               },
               polygonal = {
                 xy <- cbind(x,y)
                 result <- rep(Inf, X$n)
                 bdry <- window$bdry
                 for(i in seq(bdry)) {
                   polly <- bdry[[i]]
                   nsegs <- length(polly$x)
                   for(j in 1:nsegs) {
                     j1 <- if(j < nsegs) j + 1 else 1
                     seg <- c(polly$x[j],  polly$y[j],
                              polly$x[j1], polly$y[j1])
                     result <- pmin(result, distppl(xy, seg))
                   }
                 }
               },
               mask = {
                 b <- bdist.pixels(window, coords=F)
                 loc <- nearest.raster.point(x,y,window)
                 result <- numeric(X$n)
                 for(i in 1:(X$n)) 
			result[i] <- b[loc$row[i],loc$col[i]]
               },
               stop("Unrecognised window type", window$type)
               )
        return(result)
}

"bdist.pixels" <- function(w, ..., coords=T) {
	verifyclass(w, "owin")

        masque <- as.mask(w, ...)
        
        switch(w$type,
               mask = {
                 neg <- complement.owin(masque)
                 m <- exactPdt(neg)
                 b <- pmin(m$d,m$b)
               },
               rectangle = {
                 x <- raster.x(masque)
                 y <- raster.y(masque)
                 xmin <- w$xrange[1]
                 xmax <- w$xrange[2]
                 ymin <- w$yrange[1]
                 ymax <- w$yrange[2]
                 b <- pmin(x - xmin, xmax - x, y - ymin, ymax - y)
               },
               polygonal = {
                 # set up pixel raster
                 x <- as.vector(raster.x(masque))
                 y <- as.vector(raster.y(masque))
                 b <- rep(0, length(x))
                 # test each pixel in/out, analytically
                 inside <- inside.owin(x, y, w)
                 # compute distances for these pixels
                 xy <- cbind(x[inside], y[inside])
                 dxy <- rep(Inf, sum(inside))
                 bdry <- w$bdry
                 for(i in seq(bdry)) {
                   polly <- bdry[[i]]
                   nsegs <- length(polly$x)
                   for(j in 1:nsegs) {
                     j1 <- if(j < nsegs) j + 1 else 1
                     seg <- c(polly$x[j],  polly$y[j],
                              polly$x[j1], polly$y[j1])
                     dxy <- pmin(dxy, distppl(xy, seg))
                   }
                 }
                 b[inside] <- dxy
               },
               stop("unrecognised window type", w$type)
               )

        # reshape it
        b <- matrix(b, nrow=masque$dim[1], ncol=masque$dim[2])
        
        if(coords)
          # return in a format which can be plotted by image(), persp() etc
          return(list(x=masque$xcol, y=masque$yrow, z=t(b)))
        else
          # return matrix (for internal use by package)
          return(b)
} 

erode.mask <- function(w, r) {
  # erode a binary image mask without changing any other entries
  
	verifyclass(w, "owin")
        if(w$type != "mask")
          stop("window w is not of type \'mask\'")
        
	bb <- bdist.pixels(w, coords=F)

        if(r > max(bb))
          warning("eroded mask is empty")

        w$m <- (bb >= r)
        return(w)
}

        
erode.owin <- function(w, r, shrink.frame=T, ...) {
	verifyclass(w, "owin")

        if(2 * r >= max(diff(w$xrange), diff(w$yrange)))
          stop("erosion distance r too large for frame of window")

        # compute the dimensions of the eroded frame
        exr <- w$xrange + c(r, -r)
        eyr <- w$yrange + c(r, -r)
        ebox <- list(x=exr[c(1,2,2,1)], y=eyr[c(1,1,2,2)])

        if(w$type == "rectangle") {
          # result is a smaller rectangle
          if(shrink.frame)
            return(owin(exr, eyr))  # type 'rectangle' 
          else
            return(owin(w$xrange, w$yrange, poly=ebox)) # type 'polygonal'
        }

        # otherwise erode the window in pixel image form
        if(w$type != "mask") 
          w <- as.mask(w, ...)
        
        wnew <- erode.mask(w, r)

        if(shrink.frame) {
          # trim off some rows & columns of pixel raster
          keepcol <- (wnew$xcol >= exr[1] & wnew$xcol <= exr[2])
          keeprow <- (wnew$yrow >= eyr[1] & wnew$yrow <= eyr[2])
          wnew$xcol <- w$xcol[keepcol]
          wnew$yrow <- w$yrow[keepcol]
          wnew$dim <- c(sum(keeprow), sum(keepcol))
          wnew$m <- wnew$m[keeprow, keepcol]
        }

        return(wnew)
}	
