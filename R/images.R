#
#       images.R
#
#         $Revision: 1.13 $     $Date: 2005/01/25 23:58:22 $
#
#      The class "im" of raster images
#
# Temporary code until we sort out the class structure
#
#     im()     object creator
#
#     is.im()   tests class membership
#
#     rasterx.im(), rastery.im()    
#                      raster X and Y coordinates
#
#     nearest.pixel()   
#     lookup.im()
#                      facilities for looking up pixel values
#
################################################################
########   basic support for class "im"
################################################################
#
#   creator 

im <- function(mat, xcol=seq(ncol(mat)), yrow=seq(nrow(mat))) {
  nr <- nrow(mat)
  nc <- ncol(mat)
  if(length(xcol) != nc)
    stop("Length of xcol does not match ncol(x)")
  if(length(yrow) != nr)
    stop("Length of yrow does not match nrow(x)")
  xstep <- diff(xcol)[1]
  ystep <- diff(yrow)[1]
  xrange <- range(xcol) + c(-1,1) * xstep/2
  yrange <- range(yrow) + c(-1,1) * ystep/2
  out <- list(v   = mat,
              dim = c(nr, nc),
              xrange   = xrange,
              yrange   = yrange,
              xstep    = xstep,
              ystep    = ystep,
              xcol    = xcol,
              yrow    = yrow)
  class(out) <- "im"
  return(out)
}

is.im <- function(x) {
inherits(x,"im")
}

################################################################
########   methods for class "im"
################################################################

shift.im <- function(X, vec=c(0,0), ...) {
  verifyclass(X, "im")
  X$xrange <- X$xrange + vec[1]
  X$yrange <- X$yrange + vec[2]
  X$xcol <- X$xcol + vec[1]
  X$yrow <- X$yrow + vec[2]
  return(X)
}

"[.im" <- subset.im <-
function(x, i, j, drop=TRUE, ...) {
  if(verifyclass(i, "ppp", fatal=FALSE) && missing(j)) {
    # 'i' is a point pattern
    # Look up the greyscale values for the points of the pattern
    values <- lookup.im(x, i$x, i$y, naok=TRUE)
    if(drop) return(values[!is.na(values)]) else return(values)
  }
  if(verifyclass(i, "owin", fatal=FALSE) && missing(j)) {
    # 'i' is a window
    # if drop = FALSE, just set values outside window to NA
    # if drop = TRUE, extract values for all pixels inside window
    #                 as an image (if 'i' is a rectangle)
    #                 or as a vector (otherwise)

    xy <- expand.grid(y=x$yrow,x=x$xcol)
    inside <- inside.owin(xy$x, xy$y, i)
    if(!drop) { 
      x$v[!inside] <- NA
      return(x)
    } else if(i$type != "rectangle") {
      return(x$v[inside])
    } else {
      disjoint <- function(r, s) { (r[2] < s[1]) || (r[1] > s[2])  }
      clip <- function(r, s) { c(max(r[1],s[1]), min(r[2],s[2])) }
      inrange <- function(x, r) { (x >= r[1]) & (x <= r[2]) }
      if(disjoint(i$xrange, x$xrange) || disjoint(i$yrange, x$yrange))
        # empty intersection
        return(numeric(0))
      xr <- clip(i$xrange, x$xrange)
      yr <- clip(i$yrange, x$yrange)
      colsub <- inrange(x$xcol, xr)
      rowsub <- inrange(x$yrow, yr)
      return(im(x$v[rowsub,colsub], x$xcol[colsub], x$yrow[rowsub]))
    } 
  }
  stop("The subset operation is undefined for this type of index")
}



################################################################
########   other tools
################################################################

#
# This function is similar to nearest.raster.point except for
# the third argument 'im' and the different idiom for calculating
# row & column - which could be used in nearest.raster.point()

nearest.pixel <- function(x,y,im) {
  verifyclass(im, "im")
  nr <- im$dim[1]
  nc <- im$dim[2]
  cc <- round(1 + (x - im$xcol[1])/im$xstep)
  rr <- round(1 + (y - im$yrow[1])/im$ystep)
  cc <- pmax(1,pmin(cc, nc))
  rr <- pmax(1,pmin(rr, nr))
  return(list(row=rr, col=cc))
}

# This function is a generalisation of inside.owin()
# to images other than binary-valued images.

lookup.im <- function(im, x, y, naok=FALSE) {
  verifyclass(im, "im")

  if(length(x) != length(y))
    stop("x and y must be numeric vectors of equal length")
  value <- rep(NA, length(x))
               
  # test whether inside bounding rectangle
  xr <- im$xrange
  yr <- im$yrange
  eps <- sqrt(.Machine$double.eps)
  frameok <- (x >= xr[1] - eps) & (x <= xr[2] + eps) & 
             (y >= yr[1] - eps) & (y <= yr[2] + eps)
  
  if(!any(frameok))  # all points OUTSIDE range - no further work needed
    return(value)  # all zero

  # consider only those points which are inside the frame
  xf <- x[frameok]
  yf <- y[frameok]
  # map locations to raster (row,col) coordinates
  loc <- nearest.pixel(xf,yf,im)
  # look up image values
  vf <- im$v[cbind(loc$row, loc$col)]
  
  # insert into 'ok' vector
  value[frameok] <- vf

  if(!naok && any(is.na(value)))
    warning("Internal error: NA's generated")
  
  return(value)
}
  

rasterx.im <- function(x) {
  verifyclass(x, "im")
  v <- x$v
  xx <- x$xcol
  matrix(xx[col(v)], ncol=ncol(v), nrow=nrow(v))
}

rastery.im <- function(x) {
  verifyclass(x, "im")
  v <- x$v
  yy <- x$yrow
  matrix(yy[row(v)], ncol=ncol(v), nrow=nrow(v))
}

##############

