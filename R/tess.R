#
# tess.R
#
# support for tessellations
#
#   $Revision: 1.30 $ $Date: 2009/08/26 02:53:15 $
#
tess <- function(..., xgrid=NULL, ygrid=NULL, tiles=NULL, image=NULL,
                 window=NULL) {
  if(!is.null(window))
    win <- as.owin(window)
  else win <- NULL
  isrect <- !is.null(xgrid) && !is.null(ygrid)
  istiled <- !is.null(tiles)
  isimage <- !is.null(image)
  if(isrect + istiled + isimage != 1)
    stop("Must specify either (xgrid, ygrid) or tiles or img")
  if(isrect) {
    stopifnot(is.numeric(xgrid) && all(diff(xgrid) > 0))
    stopifnot(is.numeric(ygrid) && all(diff(ygrid) > 0))
    if(is.null(win)) win <- owin(range(xgrid), range(ygrid))
    ntiles <- (length(xgrid)-1) * (length(ygrid)-1)
    out <- list(type="rect", window=win, xgrid=xgrid, ygrid=ygrid, n=ntiles)
  } else if(istiled) {
    stopifnot(is.list(tiles))
    if(!all(unlist(lapply(tiles, is.owin))))
      stop("tiles must be a list of owin objects")
    if(is.null(win)) {
      for(i in seq(along=tiles)) {
        if(i == 1)
          win <- tiles[[1]]
        else
          win <- union.owin(win, tiles[[i]])
      }
    }
    ismask <- function(x) {x$type == "mask"}
    if(ismask(win) || any(unlist(lapply(tiles, ismask)))) {
      # convert to pixel image tessellation
      win <- as.mask(win)
      ima <- as.im(win)
      for(i in seq(along=tiles))
        ima[tiles[[i]]] <- i
      ima <- ima[win, drop=FALSE]
      ima <- eval.im(factor(ima))
      out <- list(type="image",
                  window=win, image=ima, n=length(levels(ima)))
    } else {
      # tile list
      win <- rescue.rectangle(win)
      out <- list(type="tiled", window=win, tiles=tiles, n=length(tiles))
    }
  } else if(isimage) {
    image <- eval.im(factor(image))
    if(is.null(win)) win <- as.owin(image)
    out <- list(type="image", window=win, image=image, n=length(levels(image)))
  } else stop("Internal error: unrecognised format")
  class(out) <- c("tess", class(out))
  return(out)
}

is.tess <- function(x) { inherits(x, "tess") }

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

plot.tess <- function(x, ..., main, add=FALSE, col=NULL) {
  xname <- deparse(substitute(x))
  if(missing(main))
    main <- xname
  switch(x$type,
         rect={
           win <- x$window
           if(!add)
             do.call.matched("plot.owin",
                             resolve.defaults(list(x=win, main=main),
                                              list(...)),
                             extrargs=c("sub", "lty", "lwd"))
           xg <- x$xgrid
           yg <- x$ygrid
           do.call.matched("segments",
                           resolve.defaults(list(x0=xg, y0=win$yrange[1],
                                                 x1=xg, y1=win$yrange[2]),
                                            list(col=col),
                                            list(...),
                                            .StripNull=TRUE))
           do.call.matched("segments",
                           resolve.defaults(list(x0=win$xrange[1], y0=yg,
                                                 x1=win$xrange[2], y1=yg),
                                            list(col=col),
                                            list(...),
                                            .StripNull=TRUE))
         },
         tiled={
           if(!add)
             do.call.matched("plot.owin",
                             resolve.defaults(list(x=x$window, main=main),
                                              list(...)))
           til <- tiles(x)
           plotem <- function(z, ..., col=NULL) {
             if(is.null(col))
               plot(z, ..., add=TRUE)
             else if(z$type != "mask")
               plot(z, ..., border=col, add=TRUE)
             else plot(z, ..., col=col, add=TRUE)
           }
           lapply(til, plotem, ..., col=col)
         },
         image={
           plot(x$image, main=main, ..., add=add)
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
  out <- as.listof(out)
  return(out)
}


as.im.tess <- function(X, W=NULL, ...,
                       eps=NULL, dimyx=NULL, xy=NULL,
                       na.replace=NULL) {
  # if W is present, it may have to be converted
  if(!is.null(W)) {
    stopifnot(is.owin(W))
    if(W$type != "mask")
      W <- as.mask(W, eps=eps, dimyx=dimyx, xy=xy)
  } 
  switch(X$type,
         tess={
           out <- as.im(X$image, W=W, eps=eps, dimyx=dimyx, xy=xy,
                        na.replace=na.replace)
         },
         tiled={
           if(is.null(W))
             W <- as.mask(as.owin(X), eps=eps, dimyx=dimyx, xy=xy)
           til <- X$tiles
           ntil <- length(til)
           xy <- list(x=W$xcol, y=W$yrow)
           for(i in seq(ntil)) {
             indic <- as.mask(til[[i]], xy=xy)
             tag <- as.im(indic, value=i)
             if(i == 1) {
               out <- tag
               outv <- out$v
             } else {
               outv <- pmin(outv, tag$v, na.rm=TRUE)
             }
           }
           out$v    <- outv
           out$type <- "factor"
           nama <- names(til)
           if(is.null(nama) || !all(nzchar(nama)))
             nama <- paste(seq(ntil))
           levels(out) <- nama
           unitname(out) <- unitname(W)
         },
         rect={
           if(is.null(W))
             out <- as.im(as.rectangle(X), eps=eps, dimyx=dimyx, xy=xy)
           else
             out <- as.im(W)
           xg <- X$xgrid
           yg <- X$ygrid
           nrows <- length(yg) - 1
           ncols <- length(xg) - 1
           jx <- findInterval(out$xcol, xg, rightmost.closed=TRUE)
           iy <- findInterval(out$yrow, yg, rightmost.closed=TRUE)
           M <- as.matrix(out)
           Jcol <- jx[col(M)]
           Irow <- nrows - iy[row(M)] + 1
           Ktile <- Jcol + ncols * (Irow - 1)
           out <- im(Ktile, xcol=out$xcol, yrow=out$yrow,
                     lev=seq(nrows * ncols), unitname=unitname(W))
         }
         )
  return(out)
}

as.tess <- function(X) {
  UseMethod("as.tess")
}

as.tess.tess <- function(X) {
  fields <- 
    switch(X$type,
           rect={ c("xgrid", "ygrid") },
           tiled={ "tiles" },
           image={ "image" },
           stop(paste("Unrecognised tessellation type", sQuote(X$type))))
  fields <- c(c("type", "window"), fields)
  X <- unclass(X)[fields]
  class(X) <- c("tess", class(X))
  return(X)
}

as.tess.im <- function(X) {
  return(tess(image = X))
}

as.tess.quadratcount <- function(X) {
  return(attr(X, "tess"))
}

as.tess.quadrattest <- function(X) {
  Y <- attr(X, "quadratcount")
  Z <- attr(Y, "tess")
  return(Z)
}

as.tess.list <- function(X) {
  W <- lapply(X, as.owin)
  return(tess(tiles=W))
}

as.tess.owin <- function(X) {
  return(tess(tiles=list(X)))
}

intersect.tess <- function(X, Y, ...) {
  X <- as.tess(X)
  if(is.owin(Y) && Y$type == "mask") {
    # special case
    # convert to pixel image 
    result <- as.im(Y)
    Xtiles <- tiles(X)
    for(i in seq(Xtiles)) {
      tilei <- Xtiles[[i]]
      result[tilei] <- i
    }
    result <- result[Y, drop=FALSE]
    return(tess(image=result, window=Y))
  }
  Y <- as.tess(Y)
  Xtiles <- tiles(X)
  Ytiles <- tiles(Y)
  Ztiles <- list()
  for(i in seq(Xtiles))
    for(j in seq(Ytiles)) {
      Tij <- intersect.owin(Xtiles[[i]], Ytiles[[j]], ..., fatal=FALSE)
      if(!is.null(Tij) && !is.empty(Tij))
        Ztiles <- append(Ztiles, list(Tij))
    }
  Xwin <- as.owin(X)
  Ywin <- as.owin(Y)
  Zwin <- intersect.owin(Xwin, Ywin)
  return(tess(tiles=Ztiles, window=Zwin))
}
