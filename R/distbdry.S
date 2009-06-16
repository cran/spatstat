#
#	distbdry.S		Distance to boundary
#
#	$Revision: 4.26 $	$Date: 2009/06/12 23:33:12 $
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
#       erosion.owin
#       dilation.owin
#       opening.owin
#       closing.owin
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

erode.mask <- function(w, r, strict=FALSE) {
  # erode a binary image mask without changing any other entries
  
  verifyclass(w, "owin")
  if(w$type != "mask")
    stop(paste("window w is not of type", sQuote("mask")))
  if(!is.numeric(r) || length(r) != 1)
    stop("r must be a single number")
  if(r < 0)
    stop("r must be nonnegative")
        
  bb <- bdist.pixels(w, coords=FALSE)

  if(r > max(bb))
    warning("eroded mask is empty")

  if(identical(strict, TRUE))
    w$m <- (bb > r)
  else 
    w$m <- (bb >= r)
  return(w)
}

        
erosion.owin <- erode.owin <-
  function(w, r, shrink.frame=TRUE, ..., strict=FALSE, polygonal=TRUE) {
  verifyclass(w, "owin")
  if(!is.numeric(r) || length(r) != 1)
    stop("r must be a single number")
  if(r < 0)
    stop("r must be nonnegative")

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

  if(w$type == "polygonal" && polygonal) {
    # compute polygonal region
    cw <- complement.owin(w)
    dcw <- dilate.owin(cw, r, polygonal=TRUE)
    cdcw <- complement.owin(dcw)
    wnew <- intersect.owin(cdcw, w)
    return(wnew)
  }
  
  # otherwise erode the window in pixel image form
  if(w$type == "mask") 
    wnew <- erode.mask(w, r, strict=strict)
  else {
    D <- distmap(w, invert=TRUE, ...)
    wnew <- levelset(D, r, if(strict) ">" else ">=")
  }
        
  if(shrink.frame) {
    # trim off some rows & columns of pixel raster
    keepcol <- (wnew$xcol >= exr[1] & wnew$xcol <= exr[2])
    keeprow <- (wnew$yrow >= eyr[1] & wnew$yrow <= eyr[2])
    wnew$xcol <- wnew$xcol[keepcol]
    wnew$yrow <- wnew$yrow[keeprow]
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
    stop(paste("The new rectangle",
               sQuote("rect"),
               "does not contain the window",
               sQuote("win")))
  xr <- rect$xrange
  yr <- rect$yrange
  switch(w$type,
         rectangle={
           return(owin(xr, yr,
                       poly=list(x=w$xrange[c(1,2,2,1)],
                                 y=w$yrange[c(1,1,2,2)])))
         },
         polygonal={
           return(owin(xr, yr, poly=w$bdry, check=FALSE))
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
  
    
dilation.owin <- dilate.owin <-
  function(w, r, ..., polygonal=TRUE, tight=TRUE) {
  verifyclass(w, "owin")
  if(!is.numeric(r) || length(r) != 1)
    stop("r must be a single number")
  if(r < 0)
    stop("r must be nonnegative")
  if(r == 0)
    return(w)

  ismask <- (w$type == "mask")
  polygonal <- polygonal && !ismask

  if(ismask && !any(w$m)) {
    # w is empty
    return(w)
  }
  
  # bounding frame
  bb <- if(tight) bounding.box(w) else as.rectangle(w)

  # compute dilation
  if(!polygonal) {
    # compute pixel approximation
    epsilon <- sqrt(w$xstep^2 + w$ystep^2)
    r <- max(r, epsilon)
    newbox <- grow.rectangle(bb, r)
    w <- rebound.owin(w, newbox)
    distant <- distmap(w, ...)
    dil <- levelset(distant, r, "<=")
    return(dil)
  } else {
    # compute polygonal approximation
    newbox <- grow.rectangle(bb, r)
    # extract individual edges
    w <- as.polygonal(w)
    edges <- as.psp(w)
    ends   <- edges$ends
    angles <- angles.psp(edges, directed=TRUE)
    lengths <- lengths.psp(edges)

    # initialise dilation = interior of w
    out <- w
    # dilate individual edges
    halfcircle <- seq(0, pi, length=128)[-c(1,128)]
    for(i in seq(edges$n)) {
      seg <- ends[i,]
      co <- cos(angles[i])
      si <- sin(angles[i])
      # draw sausage around i-th segment
      xx <- c(seg$x0, seg$x1) + r * si
      yy <- c(seg$y0, seg$y1) - r * co
      rightcircle <- angles[i] - pi/2 + halfcircle
      xx <- c(xx, seg$x1 + r * cos(rightcircle))
      yy <- c(yy, seg$y1 + r * sin(rightcircle))
      xx <- c(xx, c(seg$x1, seg$x0) - r * si)
      yy <- c(yy, c(seg$y1, seg$y0) + r * co)
      leftcircle <- angles[i] + pi/2 + halfcircle
      xx <- c(xx, seg$x0 + r * cos(leftcircle))
      yy <- c(yy, seg$y0 + r * sin(leftcircle))
      sausage <- owin(newbox$xrange, newbox$yrange, poly=list(x=xx, y=yy))
      # add to set
      out <- union.owin(out, sausage)
    }
    return(out)
  }
}

closing.owin <- function(w, r, ..., polygonal=TRUE) {
  if(missing(r))
    stop("r is required")
  wplus <- dilate.owin(w, r, ..., polygonal=polygonal, tight=FALSE)
  if(wplus$type == "polygonal")
    wplus <- simplify.owin(wplus, r/10)
  wclose <- erode.owin(wplus, r, strict=TRUE)
  wclose <- rebound.owin(wclose, as.rectangle(w))
  return(wclose)
}

opening.owin <- function(w, r, ..., polygonal=TRUE) {
  if(missing(r))
    stop("r is required")
  wminus <- erode.owin(w, r, ..., polygonal=polygonal, shrink.frame=FALSE)
  if(wminus$type == "polygonal")
    wminus <- simplify.owin(wminus, r/10)
  wopen <- dilate.owin(wminus, r, tight=FALSE)
  wopen <- rebound.owin(wopen, as.rectangle(w))
  return(wopen)
}

  
