#
#	ripras.S	Ripley-Rasson estimator of domain
#
#
#	$Revision: 1.5 $	$Date: 2007/03/07 03:48:39 $
#
#
#
#
#-------------------------------------
bounding.box.xy <- function(x, y=NULL) {
  xy <- xy.coords(x,y)
  owin(range(xy$x), range(xy$y))
}

convexhull.xy <- function(x, y=NULL) {
  xy <- xy.coords(x, y)
  x <- xy$x
  y <- xy$y
  h <- rev(chull(x, y))  # must be anticlockwise
  w <- owin(poly=list(x=x[h], y=y[h]))
  return(w)
}

ripras <- function(x, y=NULL, shape="convex") {
  xy <- xy.coords(x, y)
  n <- length(xy$x)
  switch(shape,
         convex = {
           w <- convexhull.xy(xy)
           # expansion factor
           m <- summary(w)$nvertices
           f <- 1/sqrt(1 - m/n)
         },
         rectangle={
           # bounding rectangle
           w <- bounding.box.xy(xy)
           # expansion factor
           f <- 1/sqrt(1 - 4/n)
         },
         stop(paste("Unrecognised option: shape=", dQuote(shape)))
         )
  # centroid
  ce <- unlist(centroid.owin(w))
  # shift centroid to origin
  W <- shift(w, -ce)
  # rescale
  W <- affine(W, mat=diag(c(f,f)))
  # shift origin to centroid
  W <- shift(W, ce)
  return(W)
}

