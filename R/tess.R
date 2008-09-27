#
# tess.R
#
# support for tessellations
#
#   $Revision: 1.4 $ $Date: 2008/09/24 17:15:26 $
#
tess <- function(..., xgrid=NULL, ygrid=NULL, tiles=NULL, image=NULL) {
  isrect <- !is.null(xgrid) && !is.null(ygrid)
  istiled <- !is.null(tiles)
  isimage <- !is.null(image)
  if(isrect + istiled + isimage != 1)
    stop("Must specify either (xgrid, ygrid) or tiles or img")
  if(isrect) {
    stopifnot(is.numeric(xgrid) && all(diff(xgrid) > 0))
    stopifnot(is.numeric(ygrid) && all(diff(ygrid) > 0))
    win <- owin(range(xgrid), range(ygrid))
    out <- list(type="rect", window=win, xgrid=xgrid, ygrid=ygrid)
  } else if(istiled) {
    stopifnot(is.list(tiles))
    if(!all(unlist(lapply(tiles, is.owin))))
      stop("tiles must be a list of owin objects")
    for(i in seq(along=tiles)) {
      if(i == 1)
        win <- tiles[[1]]
      else
        win <- union.owin(win, tiles[[i]])
    }
    win <- rescue.rectangle(win)
    out <- list(type="tiled", window=win, tiles=tiles)
  } else if(isimage) {
    image <- eval.im(factor(image))
    win <- as.owin(image)
    out <- list(type="image", window=win, image=image)
  } else stop("Internal error: unrecognised format")
  class(out) <- c("tess", class(out))
  return(out)
}

print.tess <- function(x, ...) {
  cat("Tessellation\n")
  win <- x$window
  unitinfo <- summary(unitname(win))
  switch(x$type,
         rect={
           equispaced <- function(z) {
             dz <- diff(z)
             diff(range(dz))/mean(dz) < 0.01
           }
           if(equispaced(x$xgrid) && equispaced(x$ygrid)) 
             cat(paste("Tiles are equal rectangles, of dimension",
                       signif(mean(diff(x$xgrid)), 5),
                       "x",
                       signif(mean(diff(x$ygrid)), 5),
                       unitinfo$plural, " ", unitinfo$explain,
                       "\n"))
           else
             cat(paste("Tiles are unequal rectangles\n"))
           cat(paste(paren(paste(length(x$xgrid)-1, "by",
                                 length(x$ygrid)-1,
                                 "array of tiles")),
                     "\n"))
         },
         tiled={
           if(win$type == "polygonal")
             cat("Tiles are irregular polygons\n")
           else
             cat("Tiles are windows of general type\n")
           cat(paste(length(x$tiles), "tiles\n"))
         },
         image={
           cat(paste("Tessellation is determined by",
                     "a factor-valued image",
                     "with", length(levels(x$image)), "levels\n"))
         })
  print(win)
  invisible(NULL)
}

plot.tess <- function(x, ..., main) {
  xname <- deparse(substitute(x))
  if(missing(main))
    main <- xname
  switch(x$type,
         rect={
           win <- x$window
           do.call.matched("plot.owin",
                           resolve.defaults(list(x=win, main=main),
                                            list(...)),
                           extrargs=c("sub", "lty", "lwd"))
           xg <- x$xgrid
           yg <- x$ygrid
           do.call.matched("segments",
                           resolve.defaults(list(x0=xg, y0=win$yrange[1],
                                                 x1=xg, y1=win$yrange[2]),
                                            list(...)))
           do.call.matched("segments",
                           resolve.defaults(list(x0=win$xrange[1], y0=yg,
                                                 x1=win$xrange[2], y1=yg),
                                            list(...)))
         },
         tiled={
           plotme <- function(z, ..., hatch) { plot(z, ..., hatch=FALSE) }
           plotme(x$window, main=main, ...)
           til <- tiles(x)
           lapply(til,
                  function(z, add=TRUE, ...) { plot(z, add=add, ...) },
                  ...)
         },
         image={
           plot(x$image, main=main, ...)
         })
  return(invisible(NULL))
}

"[<-.tess" <- function(x, ..., value) {
  switch(x$type,
         rect=,
         tiled={
           til <- tiles(x)
           til[...] <- value
           ok <- !unlist(lapply(til, is.null))
           x <- tess(tiles=til[ok])
         },
         image={
           stop("Cannot assign new values to subsets of a pixel image")
         })
  return(x)
}
  
"[.tess" <- function(x, ...) {
  switch(x$type,
         rect=,
         tiled={
           til <- tiles(x)[...]
           return(tess(tiles=til))
         },
         image={
           img <- x$image
           oldlev <- levels(img)
           newlev <- unique(oldlev[...])
           img <- eval.im(factor(ifelse(img %in% newlev, img, NA)))
           levels(img) <- newlev
           return(tess(image=img))
         })
}

tiles <- function(x) {
  switch(x$type,
         rect={
           out <- list()
           xg <- x$xgrid
           yg <- x$ygrid
           nx <- length(xg) - 1
           ny <- length(yg) - 1
           for(j in rev(seq(ny)))
             for(i in seq(nx)) {
               winij <- owin(xg[c(i,i+1)], yg[c(j,j+1)])
               dout <- list(winij)
               names(dout) <- paste("Tile row ", ny-j+1, ", col ", i,
                                    sep="")
               out <- append(out, dout)
             }
         },
         tiled={
           out <- x$tiles
           if(is.null(names(out)))
             names(out) <- paste("Tile", seq(length(out)))
         },
         image={
           out <- list()
           ima <- x$image
           lev <- levels(ima)
           for(i in seq(lev))
             out[[i]] <- solutionset(ima == lev[i])
           names(out) <- paste(lev)
         })
  class(out) <- c("listof", class(out))
  return(out)
}
