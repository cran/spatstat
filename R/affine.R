#
#	affine.R
#
#	$Revision: 1.35 $	$Date: 2011/10/04 10:00:44 $
#

affinexy <- function(X, mat=diag(c(1,1)), vec=c(0,0), invert=FALSE) {
  if(length(X$x) == 0 && length(X$y) == 0)
    return(list(x=numeric(0),y=numeric(0)))
  if(invert) {
    mat <- invmat <- solve(mat)
    vec <- - as.numeric(invmat %*% vec)
  }
  # Y = M X + V
  ans <- mat %*% rbind(X$x, X$y) + matrix(vec, nrow=2, ncol=length(X$x))
  return(list(x = ans[1,],
              y = ans[2,]))
}

affinexypolygon <- function(p, mat=diag(c(1,1)), vec=c(0,0),
                             detmat=det(mat)) {
  # transform (x,y)
  p[c("x","y")] <- affinexy(p, mat=mat, vec=vec)
  # transform area
  if(!is.null(p$area))
    p$area <- p$area * detmat
  # if map has negative sign, cyclic order was reversed; correct it
  if(detmat < 0)
    p <- reverse.xypolygon(p, adjust=TRUE)
  return(p)
}
       
"affine.owin" <- function(X,  mat=diag(c(1,1)), vec=c(0,0), ...,
                          rescue=TRUE) {
  verifyclass(X, "owin")
  if(!is.vector(vec) || length(vec) != 2)
    stop(paste(sQuote("vec"), "should be a vector of length 2"))
  if(!is.matrix(mat) || any(dim(mat) != c(2,2)))
    stop(paste(sQuote("mat"), "should be a 2 x 2 matrix"))
  # Inspect the determinant
  detmat <- det(mat)
  if(abs(detmat) < .Machine$double.eps)
    stop("Matrix of linear transformation is singular")
  #
  diagonalmatrix <- all(mat == diag(diag(mat)))
  scaletransform <- diagonalmatrix && (length(unique(diag(mat))) == 1)
  newunits <- if(scaletransform) unitname(X) else as.units(NULL)
  #
  switch(X$type,
         rectangle={
           if(diagonalmatrix) {
             # result is a rectangle
             Y <- owin(range(mat[1,1] * X$xrange + vec[1]),
                       range(mat[2,2] * X$yrange + vec[2]))
             unitname(Y) <- newunits
             return(Y)
           } else {
             # convert rectangle to polygon
             P <- as.polygonal(X)
             # call polygonal case
             return(affine.owin(P, mat, vec, rescue=rescue))
           }
         },
         polygonal={
           # Transform the polygonal boundaries
           bdry <- lapply(X$bdry, affinexypolygon, mat=mat, vec=vec,
                          detmat=detmat)
           # Compile result
           W <- owin(poly=bdry, check=FALSE, unitname=newunits)
           # Result might be a rectangle: if so, convert to rectangle type
           if(rescue)
             W <- rescue.rectangle(W)
           return(W)
         },
         mask={
           # binary mask
           newframe <- bounding.box.xy(affinexy(corners(X), mat, vec))
           W <- if(length(list(...)) > 0) as.mask(newframe, ...) else 
                   as.mask(newframe, eps=with(X, min(xstep, ystep)))
           pixelxy <- raster.xy(W)
           xybefore <- affinexy(pixelxy, mat, vec, invert=TRUE)
           W$m[] <- with(xybefore, inside.owin(x, y, X))
           W <- intersect.owin(W, bounding.box(W))
           if(rescue)
             W <- rescue.rectangle(W)
           return(W)
         },
         stop("Unrecognised window type")
         )
}

"affine.ppp" <- function(X, mat=diag(c(1,1)), vec=c(0,0), ...) {
  verifyclass(X, "ppp")
  r <- affinexy(X, mat, vec)
  w <- affine.owin(X$window, mat, vec, ...)
  return(ppp(r$x, r$y, window=w, marks=marks(X, dfok=TRUE), check=FALSE))
}


"affine" <- function(X, ...) {
  UseMethod("affine")
}

### ---------------------- shift ----------------------------------

"shift" <- function(X, ...) {
  UseMethod("shift")
}

shiftxy <- function(X, vec=c(0,0)) {
  list(x = X$x + vec[1],
       y = X$y + vec[2])
}

shiftxypolygon <- function(p, vec=c(0,0)) {
  # transform (x,y), retaining other data
  p[c("x","y")] <- shiftxy(p, vec=vec)
  return(p)
}

"shift.owin" <- function(X,  vec=c(0,0), ..., origin=NULL) {
  verifyclass(X, "owin")
  if(!is.null(origin)) {
    stopifnot(is.character(origin))
    if(!missing(vec))
      warning("argument vec ignored; overruled by argument origin")
    origin <- pickoption("origin", origin, c(centroid="centroid",
                                             midpoint="midpoint",
                                             bottomleft="bottomleft"))
    locn <- switch(origin,
                   centroid={ unlist(centroid.owin(X)) },
                   midpoint={ c(mean(X$xrange), mean(X$yrange)) },
                   bottomleft={ c(X$xrange[1], X$yrange[1]) })
    return(shift(X, -locn))
  }
  # Shift the bounding box
  X$xrange <- X$xrange + vec[1]
  X$yrange <- X$yrange + vec[2]
  switch(X$type,
         rectangle={
         },
         polygonal={
           # Shift the polygonal boundaries
           X$bdry <- lapply(X$bdry, shiftxypolygon, vec=vec)
         },
         mask={
           # Shift the pixel coordinates
           X$xcol <- X$xcol + vec[1]
           X$yrow <- X$yrow + vec[2]
           # That's all --- the mask entries are unchanged
         },
         stop("Unrecognised window type")
         )
  # tack on shift vector
  attr(X, "lastshift") <- vec
  # units are unchanged
  return(X)
}

"shift.ppp" <- function(X, vec=c(0,0), ..., origin=NULL) {
  verifyclass(X, "ppp")
  if(!is.null(origin)) {
    stopifnot(is.character(origin))
    if(!missing(vec))
      warning("argument vec ignored; overruled by argument origin")
    origin <- pickoption("origin", origin, c(centroid="centroid",
                                             midpoint="midpoint",
                                             bottomleft="bottomleft"))
    W <- X$window
    locn <- switch(origin,
                   centroid={ unlist(centroid.owin(W)) },
                   midpoint={ c(mean(W$xrange), mean(W$yrange)) },
                   bottomleft={ c(W$xrange[1], W$yrange[1]) })
    vec <- -locn
  }
  # perform shift
  r <- shiftxy(X, vec)
  w <- shift.owin(X$window, vec)
  Y <- ppp(r$x, r$y, window=w, marks=marks(X, dfok=TRUE), check=FALSE)
  # tack on shift vector
  attr(Y, "lastshift") <- vec
  return(Y)
}




