# replaces part of ppp.S

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
           fx <- cut(x$x, breaks=z$xgrid, include.lowest=TRUE)
           fy <- cut(x$y, breaks=z$ygrid, include.lowest=TRUE)
           m <- factor(paste(fx, fy, sep="x"))
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

