#
#	rotate.S
#
#	$Revision: 1.1 $	$Date: 2001/11/21 09:07:22 $
#

rotxy <- function(X, angle=pi/2) {
  co <- cos(angle)
  si <- sin(angle)
  list(x = co * X$x - si * X$y,
       y = si * X$x + co * X$y)
}

"rotate.owin" <- function(W, angle=pi/2) {
  verifyclass(W, "owin")
  switch(W$type,
         rectangle={
           # convert rectangle to polygon
           P <- owin(W$xrange, W$yrange, poly=
                     list(x=W$xrange[c(1,2,2,1)],
                          y=W$yrange[c(1,1,2,2)]))
           # call polygonal case
           return(rotate.owin(P, angle))
         },
         polygonal={
           # First rotate the polygonal boundaries
           bdry <- lapply(W$bdry, rotxy, angle=angle)
           # Compute bounding box of new polygons
           xr <- range(unlist(lapply(bdry, function(a) a$x)))
           yr <- range(unlist(lapply(bdry, function(a) a$y)))
           # wrap up
           return(owin(xr, yr, poly=bdry))
         },
         mask={
           stop("Sorry, \'rotate.owin\' is not yet implemented for masks")
         },
         stop("Unrecognised window type")
         )
}

"rotate.ppp" <- function(X, angle=pi/2) {
  verifyclass(X, "ppp")
  r <- rotxy(X, angle)
  w <- rotate.owin(X$window, angle)
  return(ppp(r$x, r$y, window=w, marks=X$marks))
}


"rotate" <- function(X, ...) {
  UseMethod("rotate")
}

  
