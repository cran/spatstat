#
#    summary.im.R
#
#    summary() method for class "im"
#
#    $Revision: 1.1 $   $Date: 2004/01/06 10:13:58 $
#
#    summary.im()
#    print.summary.im()
#    print.im()
#
summary.im <- function(object, ...) {
  verifyclass(object, "im")

  x <- object

  y <- unclass(x)[c("dim", "xstep", "ystep")]
  pixelarea <- y$xstep * y$ystep

  # extract image values
  v <- x$v
  inside <- !is.na(v)
  v <- v[inside]

  # summarise image values
  y$integral <- sum(v) * pixelarea
  y$mean <- mean(v)
  y$range <- range(v)
  y$min <- y$range[1]  
  y$max <- y$range[2]  
  
  # summarise pixel raster
  win <- as.owin(x)
  y$window <- summary.owin(win)

  y$fullgrid <- all(inside & win$m)

  class(y) <- "summary.im"
  return(y)
}

print.summary.im <- function(x, ...) {
  verifyclass(x, "summary.im")
  cat("Pixel image\n")
  di <- x$dim
  win <- x$window
  cat(paste(di[1], "x", di[2], "pixel array\n"))
  cat("enclosing rectangle: ")
  cat(paste("[",
            paste(win$xrange, collapse=", "),
            "] x [",
            paste(win$yrange, collapse=", "),
            "]\n"))
  cat(paste("dimensions of each pixel:", x$xstep, "x", x$ystep, "\n"))
  if(x$fullgrid) {
    cat("Image is defined on the full rectangular grid\n")
    cat(paste("Frame area = ", win$area, "\n"))
  } else {
    cat("Image is defined on a subset of the rectangular grid\n")
    cat(paste("Subset area = ", win$area, "\n"))
  }
  cat(paste("Pixel values ",
            if(x$fullgrid) "" else "(inside window)",
            ":\n", sep=""))
  cat(paste(
      "\trange = [",
      paste(x$range, collapse=","),
      "]\n",
      "\tintegral = ",
      x$integral,
      "\n",
      "\tmean = ",
      x$mean,
      "\n",
      sep=""))

  
  return(invisible(NULL))
}



print.im <- function(x, ...) {
  cat("Pixel image\n")
  di <- x$dim
  cat(paste(di[1], "x", di[2], "pixel array\n"))
  cat("enclosing rectangle: ")
  cat(paste("[",
            paste(x$xrange, collapse=", "),
            "] x [",
            paste(x$yrange, collapse=", "),
            "]\n"))
  return(invisible(NULL))
}
