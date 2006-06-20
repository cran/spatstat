#
#       images.R
#
#         $Revision: 1.20 $     $Date: 2006/06/21 04:00:14 $
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

im <- function(mat, xcol=seq(ncol(mat)), yrow=seq(nrow(mat)),
               lev=levels(mat)) {

  typ <- typeof(mat)
  if(typ == "double")
    typ <- "real"
  
  if(is.matrix(mat)) {
    nr <- nrow(mat)
    nc <- ncol(mat)
    if(length(xcol) != nc)
      stop("Length of xcol does not match ncol(mat)")
    if(length(yrow) != nr)
      stop("Length of yrow does not match nrow(mat)")
  } else {
    if(missing(xcol) || missing(yrow))
      stop(paste(sQuote("mat"),
                 "is not a matrix and I can't guess its dimensions"))
    stopifnot(length(mat) == length(xcol) * length(yrow))
    nc <- length(xcol)
    nr <- length(yrow)
  }

  # deal with factor case
  if(!is.null(lev)) {
    # convert to integer codes
    mat <- as.integer(factor(mat, levels=lev))
    typ <- "factor"
  }
  # finally coerce 'mat' to a matrix
  if(!is.matrix(mat))
    mat <- matrix(mat, nrow=nr, ncol=nc)

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
              yrow    = yrow,
              lev     = lev,
              type    = typ)
  class(out) <- "im"
  attr(out, "levels") <- lev 
  return(out)
}

is.im <- function(x) {
  inherits(x,"im")
}

"levels<-.im" <- function(x, value) {
  if(x$type != "factor") 
    stop("image is not factor-valued")
  attr(x, "levels") <- x$lev <- value
  x
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
  lev <- x$lev
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
      values <- x$v[inside]
      if(!is.null(lev))
        values <- factor(values, levels=seq(lev), labels=lev)
      return(values)
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
      return(im(x$v[rowsub,colsub], x$xcol[colsub], x$yrow[rowsub], lev))
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

lookup.im <- function(Z, x, y, naok=FALSE) {
  verifyclass(Z, "im")

  if(length(x) != length(y))
    stop("x and y must be numeric vectors of equal length")
  value <- rep(NA, length(x))
               
  # test whether inside bounding rectangle
  xr <- Z$xrange
  yr <- Z$yrange
  eps <- sqrt(.Machine$double.eps)
  frameok <- (x >= xr[1] - eps) & (x <= xr[2] + eps) & 
             (y >= yr[1] - eps) & (y <= yr[2] + eps)
  
  if(!any(frameok))  # all points OUTSIDE range - no further work needed
    return(value)  # all zero

  # consider only those points which are inside the frame
  xf <- x[frameok]
  yf <- y[frameok]
  # map locations to raster (row,col) coordinates
  loc <- nearest.pixel(xf,yf,Z)
  # look up image values
  vf <- Z$v[cbind(loc$row, loc$col)]
  
  # insert into 'ok' vector
  value[frameok] <- vf

  if(!naok && any(is.na(value)))
    warning("Internal error: NA's generated")

  # return factor, if it's a factor valued image
  if(!is.null(lev <- Z$lev))
      value <- factor(value, levels=seq(lev), labels=lev)
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

# methods for other functions

as.matrix.im <- function(x) {
  return(x$v)
}

mean.im <- function(x, ...) {
  verifyclass(x, "im")
  return(mean.default(as.matrix(x), na.rm=TRUE, ...))
}

hist.im <- function(x, ...) {
  verifyclass(x, "im")
  v <- as.numeric(as.matrix(x))
  v <- v[!is.na(v)]
  xname <- paste(deparse(substitute(x), 500), collapse="\n")
  out <- do.call("hist.default",
                 resolve.defaults(list(v),
                                  list(...),
                                  list(xlab=paste("Pixel value"),
                                       main = paste("Histogram of", xname))))
  return(invisible(out))
}


cut.im <- function(x, ...) {
  verifyclass(x, "im")
  vcut <- cut(as.numeric(as.matrix(x)), ...)
  lev <- if(is.factor(vcut)) levels(vcut) else NULL
  return(im(vcut, xcol=x$xcol, yrow=x$yrow, lev=lev))
}

quantile.im <- function(x, ...) {
  verifyclass(x, "im")
  q <- do.call("quantile",
               resolve.defaults(list(as.numeric(as.matrix(x))),
                                list(...),
                                list(na.rm=TRUE)))
  return(q)
}



