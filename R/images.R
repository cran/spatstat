#
#       images.R
#
#         $Revision: 1.1 $     $Date: 2002/07/18 10:54:34 $
#
#      The class "im" of raster images
#
# Temporary code until we sort out the class structure
#
#     im()     object creator
#
#     plot.im(), image.im(), contour.im(), persp.im()
#                      plotting functions
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

################################################################
########   methods for class "im"
################################################################

image.im <- function(x, ...) {
  image.default(x$xcol, x$yrow, t(x$v), ..., asp=1.0)
}

persp.im <- function(x, ...) {
  persp.default(x$xcol, x$yrow, t(x$v), ..., asp=1.0)
}

contour.im <- function(x, ...) {
  plot.default(x$xcol, x$yrow, type="n", asp=1.0, ...)
  contour.default(x$xcol, x$yrow, t(x$v), ..., add=TRUE)
}

plot.im <- image.im


################################################################
########   other stuff
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

lookup.im <- function(im, x, y) {
  verifyclass(im, "im")

  if(length(x) != length(y))
    stop("x and y must be numeric vectors of equal length")
  value <- rep(NA, length(x))
               
  # test whether inside bounding rectangle
  xr <- range(im$xcol)
  yr <- range(im$yrow)
  frameok <- (xr[1] <= x) & (x <= xr[2]) & (yr[1] <= y) & (y <= yr[2])
  value[!frameok] <- 0
  
  if(all(!frameok))  # all points OUTSIDE range - no further work needed
    return(value)  # all zero

  # consider only those points which are inside the frame
  xf <- x[frameok]
  yf <- y[frameok]
  # map locations to raster (row,col) coordinates
  loc <- nearest.pixel(xf,yf,im)
  # look up image values
  v <- im$v
  nf <- sum(frameok)
  vf <- logical(nf)
  for(i in 1:nf) 
    vf[i] <- v[loc$row[i],loc$col[i]]
  
  # insert into 'ok' vector
  value[frameok] <- vf

  if(any(is.na(value)))
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

