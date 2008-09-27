#
#  cut.ppp.R
#
#  cut method for ppp objects
#
#  $Revision: 1.3 $   $Date: 2008/09/24 17:15:26 $
#

cut.ppp <- function(x, z=marks(x), ...) {
  x <- as.ppp(x)
  if(missing(z)) {
    # cut the marks of x
    if(!is.marked(x))
      stop("x has no marks to cut")
    z <- marks(x, dfok=FALSE)
    m <- if(is.numeric(z)) cut(z, ...) else factor(z)
    marks(x) <- m
    return(x)
  }
  if(inherits(z, "im"))
    return(cut(x, z[x, drop=FALSE], ...))
  if(!inherits(z, "tess")) {
    # z should be a vector
    stopifnot(length(z) == x$n)
    m <- if(is.numeric(z)) cut(z, ...) else factor(z)
    marks(x) <- m
    return(x)
  }
  # z is a tessellation 
  switch(z$type,
         rect={
           jx <- findInterval(x$x, z$xgrid,
                              rightmost.closed=TRUE, all.inside=TRUE)
           iy <- findInterval(x$y, z$ygrid, 
                              rightmost.closed=TRUE, all.inside=TRUE)
           nrows    <- length(z$ygrid) - 1
           ncols <- length(z$xgrid) - 1
           jcol <- jx
           irow <- nrows - iy + 1
           ktile <- jcol + ncols * (irow - 1)
           m <- factor(ktile, levels=seq(nrows*ncols))
           ij <- expand.grid(j=seq(ncols),i=seq(nrows))
           levels(m) <- paste("Tile row ", ij$i, ", col ", ij$j, sep="")
         },
         tiled={
           todo <- seq(x$n)
           nt <- length(z$tiles)
           m <- integer(nt)
           for(i in 1:nt) {
             ti <- z$tiles[[i]]
             hit <- inside.owin(x$x[todo], x$y[todo], ti)
             if(any(hit)) {
               m[todo[hit]] <- i
               todo <- todo[!hit]
             }
             if(length(todo) == 0)
               break
           }
           m[m == 0] <- NA
           m <- factor(m, levels=seq(nt))
         },
         image={
           zim <- z$image
           m <- factor(zim[x], levels=levels(zim))
         }
         )
  marks(x) <- m
  return(x)
} 

