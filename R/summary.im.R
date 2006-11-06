#
#    summary.im.R
#
#    summary() method for class "im"
#
#    $Revision: 1.10 $   $Date: 2006/11/03 09:55:42 $
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

  # type of values?
  y$type <- x$type
  
  # factor-valued?
  lev <- x$lev
  if(fak <- !is.null(lev))
    v <- factor(v, levels=seq(lev), labels=lev)

  switch(x$type,
         integer=,
         real={
           y$integral <- sum(v) * pixelarea
           y$mean <- mean(v)
           y$range <- range(v)
           y$min <- y$range[1]  
           y$max <- y$range[2]
         },
         factor={
           y$levels <- lev
           y$table <- table(v, dnn="")
         },
         complex={
           y$integral <- sum(v) * pixelarea
           y$mean <- mean(v)
           rr <- range(Re(v))
           y$Re <- list(range=rr, min=rr[1], max=rr[2])
           ri <- range(Im(v))
           y$Im <- list(range=ri, min=ri[1], max=ri[2])
         },
         {
           # another unknown type
           pixelvalues <- v
           y$summary <- summary(pixelvalues)
         })
    
  # summarise pixel raster
  win <- as.owin(x)
  y$window <- summary.owin(win)

  y$fullgrid <- (rescue.rectangle(win)$type == "rectangle")

  class(y) <- "summary.im"
  return(y)
}

print.summary.im <- function(x, ...) {
  verifyclass(x, "summary.im")
  cat(paste(x$type, "-valued pixel image\n", sep=""))
  di <- x$dim
  win <- x$window
  cat(paste(di[1], "x", di[2], "pixel array (ny, nx)\n"))
  cat("enclosing rectangle: ")
  cat(paste("[",
            paste(win$xrange, collapse=", "),
            "] x [",
            paste(win$yrange, collapse=", "),
            "] ",
            win$units[2], "\n", sep=""))
  cat(paste("dimensions of each pixel:",
            signif(x$xstep, 3), "x", signif(x$ystep, 3),
            win$units[2], "\n"))
  if(x$fullgrid) {
    cat("Image is defined on the full rectangular grid\n")
    whatpart <- "Frame"
  } else {
    cat("Image is defined on a subset of the rectangular grid\n")
    whatpart <- "Subset"
  }
  cat(paste(whatpart, "area = ", win$area, "square", win$units[2], "\n"))
  cat(paste("Pixel values ",
            if(x$fullgrid) "" else "(inside window)",
            ":\n", sep=""))
  switch(x$type,
         integer=,
         real={
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
         },
         factor={
           print(x$table)
         },
         complex={
           cat(paste(
                     "\trange: Real [",
                     paste(x$Re$range, collapse=","),
                     "], Imaginary [",
                     paste(x$Im$range, collapse=","),
                     "]\n",
                     "\tintegral = ",
                     x$integral,
                     "\n",
                     "\tmean = ",
                     x$mean,
                     "\n",
                     sep=""))
         },
         {
           print(x$summary)
         })

  return(invisible(NULL))
}



print.im <- function(x, ...) {
  cat(paste(x$type, "-valued pixel image\n", sep=""))
  if(x$type == "factor") {
    cat("factor levels:\n")
    print(levels(x))
  }
  di <- x$dim
  cat(paste(di[1], "x", di[2], "pixel array (ny, nx)\n"))
  cat("enclosing rectangle: ")
  cat(paste("[",
            paste(x$xrange, collapse=", "),
            "] x [",
            paste(x$yrange, collapse=", "),
            "] ",
            units(x)[2], "\n", sep=""))
  return(invisible(NULL))
}
