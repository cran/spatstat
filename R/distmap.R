#
#
#      distmap.R
#
#      $Revision: 1.8 $     $Date: 2007/05/10 17:35:00 $
#
#
#     Distance transforms
#
#
distmap <- function(X, ...) {
  UseMethod("distmap")
}

distmap.ppp <- function(X, ...) {
  verifyclass(X, "ppp")
  e <- exactdt(X, ...)
  W <- e$w
  uni <- units(W)
  dmat <- e$d
  imat <- e$i
  V <- im(dmat, W$xcol, W$yrow, units=uni)
  I <- im(imat, W$xcol, W$yrow, units=uni)
  if(X$window$type == "rectangle") {
    # distance to frame boundary
    bmat <- e$b
    B <- im(bmat, W$xcol, W$yrow, units=uni)
  } else {
    # distance to window boundary, not frame boundary
    bmat <- bdist.pixels(W, coords=FALSE)
    B <- im(bmat, W$xcol, W$yrow, units=uni)
    # clip all to window
    V <- V[W, drop=FALSE]
    I <- I[W, drop=FALSE]
    B <- B[W, drop=FALSE]
  }
  attr(V, "index") <- I
  attr(V, "bdry")  <- B
  return(V)
}

distmap.owin <- function(X, ...) {
  verifyclass(X, "owin")
  X <- as.mask(X, ...)
  uni <- units(X)
  xc <- X$xcol
  yr <- X$yrow
  nr <- X$dim[1]
  nc <- X$dim[2]
# pad out the input image with a margin of width 1 on all sides
  mat <- X$m
  mat <- cbind(FALSE, mat, FALSE)
  mat <- rbind(FALSE, mat, FALSE)
# call C routine
  res <- .C("distmapbin",
            as.double(X$xrange[1]),
            as.double(X$yrange[1]),
            as.double(X$xrange[2]),
            as.double(X$yrange[2]),
            nr = as.integer(nr),
            nc = as.integer(nc),
            as.logical(t(mat)),
            distances = as.double (matrix(0, ncol = nc + 2, nrow = nr + 2)),
            boundary = as.double (matrix(0, ncol = nc + 2, nrow = nr + 2)),
            PACKAGE="spatstat"
            )
  # strip off margins again
  dist <- matrix(res$distances,
                 ncol = nc + 2, byrow = TRUE)[2:(nr + 1), 2:(nc +1)]
  bdist <- matrix(res$boundary,
                  ncol = nc + 2, byrow = TRUE)[2:(nr + 1), 2:(nc +1)]
  # cast as image objects
  V <- im(dist,  xc, yr, units=uni)
  B <- im(bdist, xc, yr, units=uni)
  attr(V, "bdry") <- B
  return(V)
}

distmap.psp <- function(X, ...) {
  verifyclass(X, "psp")
  W <- as.mask(X$window, ...)
  uni <- units(W)
  xp <- as.vector(raster.x(W))
  yp <- as.vector(raster.y(W))
  np <- length(xp)
  E <- X$ends
  z <- .C("distmap2segs",
          xp=as.double(xp),
          yp=as.double(yp),
          npoints=as.integer(np),
          x0=as.double(E$x0),
          y0=as.double(E$y0),
          x1=as.double(E$x1),
          y1=as.double(E$y1),
          nsegments=as.integer(nrow(E)),
          epsilon=as.double(.Machine$double.eps),
          dist2=as.double(numeric(np)),
          index=as.integer(integer(np)),
          PACKAGE="spatstat")
  xc <- W$xcol
  yr <- W$yrow
  Dist <- im(array(sqrt(z$dist2), dim=W$dim), xc, yr, units=uni)
  Indx <- im(array(z$index + 1, dim=W$dim), xc, yr, units=uni)
  Bdry <- im(bdist.pixels(W, coords=FALSE), xc, yr, units=uni)
  attr(Dist, "index") <- Indx
  attr(Dist, "bdry")  <- Bdry
  return(Dist)
}

