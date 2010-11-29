#
#  cut.ppp.R
#
#  cut method for ppp objects
#
#  $Revision: 1.8 $   $Date: 2010/11/25 06:06:05 $
#

cut.ppp <- function(x, z=marks(x), ...) {
  x <- as.ppp(x)
  if(missing(z) || is.null(z)) {
    z <- marks(x, dfok=TRUE)
    if(is.null(z))
      stop("x has no marks to cut")
  }
  if(is.character(z)) {
    # interpret as the name of a column of marks
    zname <- z
    m <- marks(x, dfok=TRUE)
    if(!(zname %in% colnames(m)))
      stop(paste("Unrecognised mark variable", sQuote(zname)))
    z <- m[, zname]
  }
  if(is.vector(z)) {
    stopifnot(length(z) == npoints(x))
    g <- if(is.numeric(z)) cut(z, ...) else factor(z)
    marks(x) <- g
    return(x)
  }
  if(is.data.frame(z) || is.matrix(z)) {
    stopifnot(nrow(z) == npoints(x))
    # take first column 
    z <- z[,1]
    g <- if(is.numeric(z)) cut(z, ...) else factor(z)
    marks(x) <- g
    return(x)
  }
  if(is.im(z)) 
    return(cut(x, z[x, drop=FALSE], ...))

  if(is.tess(z)) {
    switch(z$type,
           rect={
             jx <- findInterval(x$x, z$xgrid, rightmost.closed=TRUE)
             iy <- findInterval(x$y, z$ygrid, rightmost.closed=TRUE)
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
             m <- integer(x$n)
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
             m <- factor(zim[x, drop=FALSE], levels=levels(zim))
           }
           )
    marks(x) <- m
    return(x)
  }
  stop("Format of z not understood")
} 

