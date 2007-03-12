#
#	affine.S
#
#	$Revision: 1.13 $	$Date: 2007/03/06 09:43:38 $
#

affinexy <- function(X, mat=diag(c(1,1)), vec=c(0,0)) {
  if(length(X$x) == 0 && length(X$y) == 0)
    return(list(x=c(),y=c()))
  # Y = M X + V
  ans <- mat %*% rbind(X$x, X$y) + matrix(vec, nrow=2, ncol=length(X$x))
  return(list(x = ans[1,],
              y = ans[2,]))
}

"affine.owin" <- function(X,  mat=diag(c(1,1)), vec=c(0,0), ...) {
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
  newunits <- if(scaletransform) units(X) else as.units(NULL)
  #
  switch(X$type,
         rectangle={
           if(diagonalmatrix) {
             # result is a rectangle
             Y <- owin(range(mat[1,1] * X$xrange + vec[1]),
                       range(mat[2,2] * X$yrange + vec[2]))
             units(Y) <- newunits
             return(Y)
           } else {
             # convert rectangle to polygon
             P <- as.polygonal(X)
             # call polygonal case
             return(affine.owin(P, mat, vec))
           }
         },
         polygonal={
           # First transform the polygonal boundaries
           bdry <- lapply(X$bdry, affinexy, mat=mat, vec=vec)
           # If determinant < 0, traverse polygons in reverse direction
           if(detmat < 0)
             bdry <- lapply(bdry, reverse.xypolygon)
           # Compute bounding box of new polygons
           xr <- range(unlist(lapply(bdry, function(a) a$x)))
           yr <- range(unlist(lapply(bdry, function(a) a$y)))
           # Compile result
           W <- owin(xr, yr, poly=bdry)
           # Result might be a rectangle: if so, convert to rectangle type
           W <- rescue.rectangle(W)
           units(W) <- newunits
           return(W)
         },
         mask={
           stop(paste("Sorry,", sQuote("affine.owin"),
                      "is not yet implemented for masks"))
         },
         stop("Unrecognised window type")
         )
}

"affine.ppp" <- function(X, mat=diag(c(1,1)), vec=c(0,0), ...) {
  verifyclass(X, "ppp")
  r <- affinexy(X, mat, vec)
  w <- affine.owin(X$window, mat, vec)
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

"shift.owin" <- function(X,  vec=c(0,0), ...) {
  verifyclass(X, "owin")
  # Shift the bounding box
  X$xrange <- X$xrange + vec[1]
  X$yrange <- X$yrange + vec[2]
  switch(X$type,
         rectangle={
         },
         polygonal={
           # Shift the polygonal boundaries
           X$bdry <- lapply(X$bdry, shiftxy, vec=vec)
         },
         mask={
           # Shift the pixel coordinates
           X$xcol <- X$xcol + vec[1]
           X$yrow <- X$yrow + vec[2]
           # That's all --- the mask entries are unchanged
         },
         stop("Unrecognised window type")
         )
  # units are unchanged
  return(X)
}

"shift.ppp" <- function(X, vec=c(0,0), ...) {
  verifyclass(X, "ppp")
  r <- shiftxy(X, vec)
  w <- shift.owin(X$window, vec)
  return(ppp(r$x, r$y, window=w, marks=marks(X), check=FALSE))
}



