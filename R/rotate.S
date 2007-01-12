#
#	rotate.S
#
#	$Revision: 1.7 $	$Date: 2006/12/12 08:57:55 $
#

rotxy <- function(X, angle=pi/2) {
  co <- cos(angle)
  si <- sin(angle)
  list(x = co * X$x - si * X$y,
       y = si * X$x + co * X$y)
}

"rotate.owin" <- function(X, angle=pi/2, ...) {
  verifyclass(X, "owin")
  switch(X$type,
         rectangle={
           # convert rectangle to polygon
           P <- owin(X$xrange, X$yrange, poly=
                     list(x=X$xrange[c(1,2,2,1)],
                          y=X$yrange[c(1,1,2,2)]),
                     units=units(X))
           # call polygonal case
           return(rotate.owin(P, angle))
         },
         polygonal={
           # First rotate the polygonal boundaries
           bdry <- lapply(X$bdry, rotxy, angle=angle)
           # Compute bounding box of new polygons
           xr <- range(unlist(lapply(bdry, function(a) a$x)))
           yr <- range(unlist(lapply(bdry, function(a) a$y)))
           # wrap up
           W <- owin(xr, yr, poly=bdry, units=units(X))
           W <- rescue.rectangle(W)
           return(W)
         },
         mask={
           stop(paste("Sorry,", sQuote("rotate.owin"),
                      "is not yet implemented for masks"))
         },
         stop("Unrecognised window type")
         )
}

"rotate.ppp" <- function(X, angle=pi/2, ...) {
  verifyclass(X, "ppp")
  r <- rotxy(X, angle)
  w <- rotate.owin(X$window, angle)
  return(ppp(r$x, r$y, window=w, marks=marks(X, dfok=TRUE), check=FALSE))
}


"rotate" <- function(X, ...) {
  UseMethod("rotate")
}

  
