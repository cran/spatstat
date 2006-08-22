#
#	distbdry.S		Distance to boundary
#
#	$Revision: 4.9 $	$Date: 2006/08/09 09:56:44 $
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
        if(X$n == 0)
          return(numeric(0))
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
                 b <- bdist.pixels(window, coords=FALSE)
                 loc <- nearest.raster.point(x,y,window)
                 result <- b[cbind(loc$row, loc$col)]
               },
               stop("Unrecognised window type", window$type)
               )
        return(result)
}

"bdist.pixels" <- function(w, ..., coords=TRUE) {
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
        
	bb <- bdist.pixels(w, coords=FALSE)

        if(r > max(bb))
          warning("eroded mask is empty")

        w$m <- (bb >= r)
        return(w)
}

        
erode.owin <- function(w, r, shrink.frame=TRUE, ...) {
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
          wnew$yrow <- w$yrow[keeprow]
          wnew$dim <- c(sum(keeprow), sum(keepcol))
          wnew$m <- wnew$m[keeprow, keepcol]
          wnew$xrange <- exr
          wnew$yrange <- eyr
        }

        return(wnew)
}	

rebound.owin <- function(w, rect) {
  verifyclass(w, "owin")
  verifyclass(rect, "owin")
  bb <- bounding.box(w)
  if(!is.subset.owin(bb, rect))
    stop("The new rectangle \`rect\' does not contain the window \`win\'")
  xr <- rect$xrange
  yr <- rect$yrange
  switch(w$type,
         rectangle={
           return(owin(xr, yr,
                       poly=list(x=w$xrange[c(1,2,2,1)],
                                 y=w$yrange[c(1,1,2,2)])))
         },
         polygonal={
           return(owin(xr, yr, poly=w$bdry))
         },
         mask={
           newseq <- function(oldseq, newrange, dstep) {
             oldends <- range(oldseq)
             nleft <- max(0, floor((oldends[1] - newrange[1])/dstep))
             nright <- max(0, floor((newrange[2] - oldends[2])/dstep))
             newstart <- max(oldends[1] - nleft * dstep, newrange[1])
             newend <- min(oldends[2] + nright * dstep, newrange[2])
             seq(from=newstart, by=dstep, to=newend)
           }
           xcol <- newseq(w$xcol, xr, mean(diff(w$xcol)))
           yrow <- newseq(w$yrow, yr, mean(diff(w$yrow)))
           newmask <- as.mask(xy=list(x=xcol, y=yrow))
           xx <- raster.x(newmask)
           yy <- raster.y(newmask)
           newmask$m <- inside.owin(xx, yy, w)
           return(newmask)
         }
         )
}
  
    
dilate.owin <- function(w, r, ...) {
  verifyclass(w, "owin")
  if(w$type != "mask")
    w <- as.mask(w, ...)
  bb <- bounding.box(w)
  newbox <- owin(bb$xrange + c(-r,r), bb$yrange + c(-r,r))
  w <- rebound.owin(w, newbox)
  distant <- distmap(w)
  dil <- w
  dil$m <- (distant$v <= r)
  return(dil)
}

  
  
  
