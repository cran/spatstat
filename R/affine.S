#
#	affine.S
#
#	$Revision: 1.1 $	$Date: 2001/11/21 09:07:22 $
#

affinexy <- function(X, mat=diag(c(1,1)), vec=c(0,0)) {
  # Y = M X + V
  ans <- mat %*% rbind(X$x, X$y) + matrix(vec, nrow=2, ncol=length(X$x))
  list(x = ans[1,],
       y = ans[2,])
}

"affine.owin" <- function(W,  mat=diag(c(1,1)), vec=c(0,0)) {
  verifyclass(W, "owin")
  # Inspect the determinant
  detmat <- det(mat)
  if(abs(detmat) < .Machine$double.eps)
    stop("Matrix of linear transformation is singular")
  #
  switch(W$type,
         rectangle={
           # convert rectangle to polygon
           P <- owin(W$xrange, W$yrange, poly=
                     list(x=W$xrange[c(1,2,2,1)],
                          y=W$yrange[c(1,1,2,2)]))
           # call polygonal case
           return(affine.owin(P, mat, vec))
         },
         polygonal={
           # First transform the polygonal boundaries
           bdry <- lapply(W$bdry, affinexy, mat=mat, vec=vec)
           # If determinant < 0, traverse polygons in reverse direction
           if(detmat < 0)
             bdry <- lapply(bdry, reverse.xypolygon)
           # Compute bounding box of new polygons
           xr <- range(unlist(lapply(bdry, function(a) a$x)))
           yr <- range(unlist(lapply(bdry, function(a) a$y)))
           # wrap up
           return(owin(xr, yr, poly=bdry))
         },
         mask={
           stop("Sorry, \'affine.owin\' is not yet implemented for masks")
         },
         stop("Unrecognised window type")
         )
}

"affine.ppp" <- function(X, mat=diag(c(1,1)), vec=c(0,0)) {
  verifyclass(X, "ppp")
  r <- affinexy(X, mat, vec)
  w <- affine.owin(X$window, mat, vec)
  return(ppp(r$x, r$y, window=w, marks=X$marks))
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

"shift.owin" <- function(W,  vec=c(0,0)) {
  verifyclass(W, "owin")
  # Shift the bounding box
  W$xrange <- W$xrange + vec[1]
  W$yrange <- W$yrange + vec[2]
  switch(W$type,
         rectangle={
         },
         polygonal={
           # Shift the polygonal boundaries
           W$bdry <- lapply(W$bdry, shiftxy, vec=vec)
         },
         mask={
           # Shift the pixel coordinates
           W$xcol <- W$xcol + vec[1]
           W$yrow <- W$yrow + vec[2]
           # That's all --- the mask entries are unchanged
         },
         stop("Unrecognised window type")
         )
  return(W)
}

"shift.ppp" <- function(X, vec=c(0,0)) {
  verifyclass(X, "ppp")
  r <- shiftxy(X, vec)
  w <- shift.owin(X$window, vec)
  return(ppp(r$x, r$y, window=w, marks=X$marks))
}


