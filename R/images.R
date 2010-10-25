#
#       images.R
#
#         $Revision: 1.75 $     $Date: 2010/10/21 03:26:05 $
#
#      The class "im" of raster images
#
#     im()     object creator
#
#     is.im()   tests class membership
#
#     rasterx.im(), rastery.im()    
#                      raster X and Y coordinates
#
#     nearest.pixel()   
#     lookup.im()
#                      facilities for looking up pixel values
#
################################################################
########   basic support for class "im"
################################################################
#
#   creator 

im <- function(mat, xcol=seq(ncol(mat)), yrow=seq(nrow(mat)), 
               xrange=NULL, yrange=NULL,
               unitname=NULL) {

  typ <- typeof(mat)
  if(typ == "double")
    typ <- "real"

  miss.xcol <- missing(xcol)
  miss.yrow <- missing(yrow)
  
  # determine dimensions
  if(is.matrix(mat)) {
    nr <- nrow(mat)
    nc <- ncol(mat)
    if(length(xcol) != nc)
      stop("Length of xcol does not match ncol(mat)")
    if(length(yrow) != nr)
      stop("Length of yrow does not match nrow(mat)")
  } else {
    if(miss.xcol || miss.yrow)
      stop(paste(sQuote("mat"),
                 "is not a matrix and I can't guess its dimensions"))
    stopifnot(length(mat) == length(xcol) * length(yrow))
    nc <- length(xcol)
    nr <- length(yrow)
  }

  # deal with factor case
  if(is.factor(mat)) {
    typ <- "factor"
  } else if(!is.null(lev <- levels(mat))) {
    typ <- "factor"
    mat <- factor(mat, levels=lev)
  }

  # Ensure 'mat' is a matrix (without destroying factor information)
  if(!is.matrix(mat))
    dim(mat) <- c(nr, nc)

  # set up coordinates
  if(miss.xcol && !is.null(xrange)) {
    # use 'xrange' 
    xstep <- diff(xrange)/nc
    xcol <- xrange[1] - xstep/2 + xstep * seq(nc)
  } else {
    # use 'xcol'
    xstep <- diff(xcol)[1]
    xrange <- range(xcol) + c(-1,1) * xstep/2
  }
  if(miss.yrow && !is.null(yrange)) {
    # use 'yrange'
    ystep <- diff(yrange)/nr
    yrow <- yrange[1] - ystep/2 + ystep * seq(nr)
  } else {
    # use 'yrow'
    ystep <- diff(yrow)[1]
    yrange <- range(yrow) + c(-1,1) * ystep/2
  }  
  unitname <- as.units(unitname)

  # get rid of those annoying 8.67e-19 printouts
  swat <- function(x) {ifelse(abs(x) < .Machine$double.eps, 0, x)}
  xrange <- swat(xrange)
  yrange <- swat(yrange)
  
  out <- list(v   = mat,
              dim = c(nr, nc),
              xrange   = xrange,
              yrange   = yrange,
              xstep    = xstep,
              ystep    = ystep,
              xcol    = xcol,
              yrow    = yrow,
              type    = typ,
              units   = unitname)
  class(out) <- "im"
  return(out)
}

is.im <- function(x) {
  inherits(x,"im")
}

levels.im <- function(x) {
  levels(x$v)
}

"levels<-.im" <- function(x, value) {
  if(x$type != "factor") 
    stop("image is not factor-valued")
  levels(x$v) <- value
  x
}

################################################################
########   methods for class "im"
################################################################

shift.im <- function(X, vec=c(0,0), ..., origin=NULL) {
  verifyclass(X, "im")
  if(!is.null(origin)) {
    stopifnot(is.character(origin))
    if(!missing(vec))
      warning("argument vec ignored; overruled by argument origin")
    origin <- pickoption("origin", origin, c(centroid="centroid",
                                             midpoint="midpoint",
                                             bottomleft="bottomleft"))
    W <- as.owin(X)
    locn <- switch(origin,
                   centroid={ unlist(centroid.owin(W)) },
                   midpoint={ c(mean(W$xrange), mean(W$yrange)) },
                   bottomleft={ c(W$xrange[1], W$yrange[1]) })
    return(shift(X, -locn))
  }
  X$xrange <- X$xrange + vec[1]
  X$yrange <- X$yrange + vec[2]
  X$xcol <- X$xcol + vec[1]
  X$yrow <- X$yrow + vec[2]
  return(X)
}

"[.im" <- 
function(x, i, drop=TRUE, ..., raster=NULL) {
  
  if(missing(i)) {
    # entire image 
    out <- if(is.null(raster)) x else as.im(raster)
    xy <- expand.grid(y=out$yrow,x=out$xcol)
    if(!is.null(raster)) {
      # resample image on new pixel raster
      values <- lookup.im(x, xy$x, xy$y, naok=TRUE)
      out <- im(values, out$xcol, out$yrow, unitname=unitname(out))
    }
    if(!drop)
      return(out)
    else {
      v <- out$v
      return(v[!is.na(v)])
    }
  }
  if(verifyclass(i, "owin", fatal=FALSE)) {
    # 'i' is a window
    # if drop = FALSE, just set values outside window to NA
    # if drop = TRUE, extract values for all pixels inside window
    #                 as an image (if 'i' is a rectangle)
    #                 or as a vector (otherwise)

    out <- if(is.null(raster)) x else as.im(raster)
    xy <- expand.grid(y=out$yrow,x=out$xcol)
    if(!is.null(raster)) {
      # resample image on new pixel raster
      values <- lookup.im(x, xy$x, xy$y, naok=TRUE)
      out <- im(values, out$xcol, out$yrow, unitname=unitname(out))
    }
    inside <- inside.owin(xy$x, xy$y, i)
    if(!drop) { 
      out$v[!inside] <- NA
      return(out)
    } else if(i$type != "rectangle") {
      values <- out$v[inside]
      return(values)
    } else {
      disjoint <- function(r, s) { (r[2] < s[1]) || (r[1] > s[2])  }
      clip <- function(r, s) { c(max(r[1],s[1]), min(r[2],s[2])) }
      inrange <- function(x, r) { (x >= r[1]) & (x <= r[2]) }
      if(disjoint(i$xrange, x$xrange) || disjoint(i$yrange, x$yrange))
        # empty intersection
        return(numeric(0))
      xr <- clip(i$xrange, x$xrange)
      yr <- clip(i$yrange, x$yrange)
      colsub <- inrange(out$xcol, xr)
      rowsub <- inrange(out$yrow, yr)
      return(im(out$v[rowsub,colsub], out$xcol[colsub], out$yrow[rowsub],
                unitname=unitname(x)))
    } 
  }
  if(verifyclass(i, "im", fatal=FALSE)) {
    # logical images OK
    if(i$type == "logical") {
      # convert to window
      w <- as.owin(eval.im(ifelse(i, 1, NA)))
      return(x[w, drop=drop, ..., raster=raster])
    } else stop("Subset argument \'i\' is an image, but not of logical type")
  }

  # Try indexing as a matrix
  y <- try(as.matrix(x)[i, drop=FALSE], silent=TRUE)
  if(!inherits(y, "try-error")) {
    # valid subset index for a matrix
    # check whether it's a rectangular block, in correct order
    RR <- row(x$v)
    CC <- col(x$v)
    rr <- RR[i]
    cc <- CC[i]
    rseq <- sort(unique(rr))
    cseq <- sort(unique(cc))
    if(all(diff(rseq) == 1) && all(diff(cseq) == 1) &&
       all(rr == RR[rseq, cseq]) && all(cc == CC[rseq,cseq])) {
      # yes - make image
      dim(y) <- c(length(rseq), length(cseq))
      Y <- x
      Y$v <- y
      Y$dim <- dim(y)
      Y$xcol <- x$xcol[cseq]
      Y$yrow <- x$yrow[rseq]
      Y$xrange <- range(Y$xcol) + c(-1,1) * x$xstep/2
      Y$yrange <- range(Y$yrow) + c(-1,1) * x$ystep/2
      return(Y)
    }
    # return pixel values (possibly as matrix)
    return(y)
  }

  if(!is.null(ip <- as.ppp(i, W=as.owin(x), fatal=FALSE, check=FALSE))) {
    # 'i' is a point pattern 
    # Look up the greyscale values for the points of the pattern
    values <- lookup.im(x, ip$x, ip$y, naok=TRUE)
    if(drop) 
      values <- values[!is.na(values)]
    if(length(values) == 0) 
      # ensure the zero-length vector is of the right type
      values <- 
      switch(x$type,
             factor={ factor(, levels=levels(x)) },
             integer = { integer(0) },
             logical = { logical(0) },
             real = { numeric(0) },
             complex = { complex(0) },
             character = { character(0) },
             { values }
             )
    return(values)
  }
  stop("The subset operation is undefined for this type of index")
}


"[<-.im" <- function(x, i, value) {
  X <- x
  W <- as.owin(X)
  if(missing(i)) {
    # set all pixels to 'value'
    v <- X$v
    v[!is.na(v)] <- value
    X$v <- v
    return(X)
  }
  if(verifyclass(i, "owin", fatal=FALSE)) {
    # 'i' is a window
    if(is.empty(i))
      return(X)
    xx <- as.vector(raster.x(W))
    yy <- as.vector(raster.y(W))
    ok <- inside.owin(xx, yy, i)
    X$v[ok] <- value
    return(X)
  }
  if(verifyclass(i, "im", fatal=FALSE) && i$type == "logical") {
    # convert logical vector to window where entries are TRUE
    i <- as.owin(eval.im(ifelse(i, 1, NA)))
    # continue as above
    xx <- as.vector(raster.x(W))
    yy <- as.vector(raster.y(W))
    ok <- inside.owin(xx, yy, i)
    X$v[ok] <- value
    return(X)
  }
  # try indexing as a matrix
  out <- try(X$v[i] <- value, silent=TRUE)
  if(!inherits(out, "try-error")) 
    return(X)

  if(!is.null(ip <- as.ppp(i, W=W, fatal=FALSE, check=TRUE))) {
    # 'i' is a point pattern
    # test whether all points are inside window FRAME
    ok <- inside.owin(ip$x, ip$y, as.rectangle(W))
    if(any(!ok)) {
      warning("Some points are outside the outer frame of the image")
      if(length(value) == ip$n)
        value <- value[ok]
      ip <- ip[ok]
    }
    # determine row & column positions for each point 
    loc <- nearest.pixel(ip$x, ip$y, X)
    # set values
    X$v[cbind(loc$row, loc$col)] <- value
    return(X)
  }
  stop("The subset operation is undefined for this type of index")
}

################################################################
########   other tools
################################################################

#
# This function is similar to nearest.raster.point except for
# the third argument 'im' and the different idiom for calculating
# row & column - which could be used in nearest.raster.point()

nearest.pixel <- function(x,y,im) {
  verifyclass(im, "im")
  nr <- im$dim[1]
  nc <- im$dim[2]
  cc <- round(1 + (x - im$xcol[1])/im$xstep)
  rr <- round(1 + (y - im$yrow[1])/im$ystep)
  cc <- pmax(1,pmin(cc, nc))
  rr <- pmax(1,pmin(rr, nr))
  return(list(row=rr, col=cc))
}

# Explores the 3 x 3 neighbourhood of nearest.pixel
# and finds the nearest pixel that is not NA

nearest.valid.pixel <- function(x,y,im) {
  rc <- nearest.pixel(x,y,im)
  rr <- rc$row
  cc <- rc$col
  # check whether any pixels are outside image domain
  outside <- is.na(im$v)
  miss <- outside[cbind(rr, cc)]
  if(!any(miss))
    return(rc)
  # for offending pixels, explore 3 x 3 neighbourhood
  nr <- im$dim[1]
  nc <- im$dim[2]
  xcol <- im$xcol
  yrow <- im$yrow
  for(i in which(miss)) {
    rows <- rr[i] + c(-1,0,1)
    cols <- cc[i] + c(-1,0,1)
    rows <- unique(pmax(1, pmin(rows, nr)))
    cols <- unique(pmax(1, pmin(cols, nc)))
    rcp <- expand.grid(row=rows, col=cols)
    ok <- !outside[as.matrix(rcp)]
    if(any(ok)) {
      # At least one of the neighbours is valid
      # Find the closest one
      rcp <- rcp[ok,]
      dsq <- with(rcp, (x[i] - xcol[col])^2 + (y[i] - yrow[row])^2)
      j <- which.min(dsq)
      rc$row[i] <- rcp$row[j]
      rc$col[i] <- rcp$col[j]
    }
  }
  return(rc)
}
  

# This function is a generalisation of inside.owin()
# to images other than binary-valued images.

lookup.im <- function(Z, x, y, naok=FALSE, strict=TRUE) {
  verifyclass(Z, "im")

  if(Z$type == "factor")
    Z <- repair.old.factor.image(Z)
  
  if(length(x) != length(y))
    stop("x and y must be numeric vectors of equal length")

  # initialise answer to NA 
  if(Z$type != "factor") {
    niets <- NA
    mode(niets) <- mode(Z$v)
  } else {
    niets <- factor(NA, levels=levels(Z))
  }
  value <- rep(niets, length(x))
               
  # test whether inside bounding rectangle
  xr <- Z$xrange
  yr <- Z$yrange
  eps <- sqrt(.Machine$double.eps)
  frameok <- (x >= xr[1] - eps) & (x <= xr[2] + eps) & 
             (y >= yr[1] - eps) & (y <= yr[2] + eps)
  
  if(!any(frameok)) {
    # all points OUTSIDE range - no further work needed
    if(!naok)
      warning("Internal error: all values NA")
    return(value)  # all NA
  }

  # consider only those points which are inside the frame
  xf <- x[frameok]
  yf <- y[frameok]
  # map locations to raster (row,col) coordinates
  if(strict)
    loc <- nearest.pixel(xf,yf,Z)
  else
    loc <- nearest.valid.pixel(xf,yf,Z)
  # look up image values
  vf <- Z$v[cbind(loc$row, loc$col)]
  
  # insert into answer
  value[frameok] <- vf

  if(!naok && any(is.na(value)))
    warning("Internal error: NA's generated")

  return(value)
}
  

rasterx.im <- function(x) {
  verifyclass(x, "im")
  v <- x$v
  xx <- x$xcol
  matrix(xx[col(v)], ncol=ncol(v), nrow=nrow(v))
}

rastery.im <- function(x) {
  verifyclass(x, "im")
  v <- x$v
  yy <- x$yrow
  matrix(yy[row(v)], ncol=ncol(v), nrow=nrow(v))
}

##############

# methods for other functions

xtfrm.im <- function(x) { as.numeric(as.matrix.im(x)) }

as.matrix.im <- function(x, ...) {
  return(x$v)
}

as.data.frame.im <- function(x, ...) {
  verifyclass(x, "im")
  v <- x$v
  xx <- x$xcol[col(v)]
  yy <- x$yrow[row(v)]
  ok <- !is.na(v)
  xx <- as.vector(xx[ok])
  yy <- as.vector(yy[ok])
  # extract pixel values without losing factor info
  vv <- v[ok]
  dim(vv) <- NULL
  # 
  data.frame(x=xx, y=yy, value=vv, ...)
}
  
mean.im <- function(x, ...) {
  verifyclass(x, "im")
  xvalues <- x[drop=TRUE]
  return(mean(xvalues))
}

sum.im <- function(x, ...) {
  verifyclass(x, "im")
  xvalues <- x[drop=TRUE]
  return(sum(xvalues, ...))
}

median.im <- function(x, ...) {
  verifyclass(x, "im")
  xvalues <- x[drop=TRUE]
  return(median(xvalues))
}

range.im <- function(x, ...) {
  verifyclass(x, "im")
  xvalues <- x[drop=TRUE]
  return(range(xvalues))
}

max.im <- function(x, ...) {
  verifyclass(x, "im")
  xvalues <- x[drop=TRUE]
  return(max(xvalues, ...))
}

min.im <- function(x, ...) {
  verifyclass(x, "im")
  xvalues <- x[drop=TRUE]
  return(min(xvalues, ...))
}

hist.im <- function(x, ..., probability=FALSE) {
  xname <- paste(deparse(substitute(x), 500), collapse="\n")
  verifyclass(x, "im")
  main <- paste("Histogram of", xname)
  # default plot arguments
  # extract pixel values
  values <- as.matrix(x)
  dim(values) <- NULL
  # barplot or histogram
  if(x$type %in% c("logical", "factor")) {
    # barplot
    tab <- table(values)
    probs <- tab/sum(tab)
    if(probability) {
      heights <- probs
      ylab <- "Probability"
    } else {
      heights <- tab
      ylab <- "Number of pixels"
    }
    mids <- do.call("barplot",
                   resolve.defaults(list(heights),
                                    list(...),
                                    list(xlab=paste("Pixel value"),
                                         ylab=ylab,
                                         main=main)))
    out <- list(counts=tab, probs=probs, heights=heights,
                mids=mids, xname=xname)
    class(out) <- "barplotdata"
  } else {
    # histogram
    values <- values[!is.na(values)]
    plotit <- resolve.defaults(list(...), list(plot=TRUE))$plot
    if(plotit) {
      ylab <- if(probability) "Probability density" else "Number of pixels"
      out <- do.call("hist.default",
                   resolve.defaults(list(values),
                                    list(...),
                                    list(probability=probability),
                                    list(xlab=paste("Pixel value"),
                                         ylab=ylab,
                                         main=main)))
      out$xname <- xname
    } else {
      # plot.default whinges if `probability' given when plot=FALSE
      out <- do.call("hist.default",
                   resolve.defaults(list(values),
                                    list(...)))
      # hack!
      out$xname <- xname
    }
  }
  return(invisible(out))
}

plot.barplotdata <- function(x, ...) {
  do.call("barplot",
          resolve.defaults(list(height=x$heights),
                           list(...),
                           list(main=paste("Histogram of ", x$xname))))
}

cut.im <- function(x, ...) {
  verifyclass(x, "im")
  vcut <- cut(as.numeric(as.matrix(x)), ...)
  return(im(vcut, xcol=x$xcol, yrow=x$yrow, unitname=unitname(x)))
}

quantile.im <- function(x, ...) {
  verifyclass(x, "im")
  q <- do.call("quantile",
               resolve.defaults(list(as.numeric(as.matrix(x))),
                                list(...),
                                list(na.rm=TRUE)))
  return(q)
}

integral.im <- function(x, ...) {
  verifyclass(x, "im")
  typ <- x$type
  if(!any(typ == c("integer", "real", "complex", "logical")))
    stop(paste("Don't know how to integrate an image of type", sQuote(typ)))
  a <- with(x, sum(v, na.rm=TRUE) * xstep * ystep)
  return(a)
}

conform.imagelist <- function(X, Zlist) {
  # determine points of X where all images in Zlist are defined
  ok <- rep(TRUE, length(X$x))
  for(i in seq(Zlist)) {
    Zi <- Zlist[[i]]
    ZiX <- Zi[X, drop=FALSE]
    ok <- ok & !is.na(ZiX)
  }
  return(ok)
}

split.im <- function(x, f, ..., drop=FALSE) {
  stopifnot(is.im(x))
  if(inherits(f, "tess")) 
    subsets <- tiles(f)
  else if(is.im(f)) {
    if(f$type != "factor")
      f <- eval.im(factor(f))
    subsets <- tiles(tess(image=f))
  } else stop("f should be a tessellation or a factor-valued image")
  if(!is.subset.owin(as.owin(x), as.owin(f)))
    stop("f does not cover the window of x")
  n <- length(subsets)
  out <- vector(mode="list", length=n)
  names(out) <- names(subsets)
  for(i in 1:n)
    out[[i]] <- x[subsets[[i]], drop=drop]
  if(drop)
    return(out)
  else 
    return(as.listof(out))
}

by.im <- function(data, INDICES, FUN, ...) {
  stopifnot(is.im(data))
  V <- split(data, INDICES)
  U <- lapply(V, FUN, ...)
  return(as.listof(U))
}

rebound.im <- function(x, rect) {
  stopifnot(is.im(x))
  stopifnot(is.owin(rect))
  rect <- as.rectangle(rect)
  stopifnot(is.subset.owin(as.rectangle(x), rect))
  # compute number of extra rows/columns
  dx <- x$xstep
  nleft  <- max(0, floor((x$xrange[1]-rect$xrange[1])/dx))
  nright <- max(0, floor((rect$xrange[2]-x$xrange[2])/dx))
  dy <- x$ystep
  nbot <- max(0, floor((x$yrange[1]-rect$yrange[1])/dy))
  ntop <- max(0, floor((rect$yrange[2]-x$yrange[2])/dy))
  # expand pixel data matrix
  nr <- x$dim[1]
  nc <- x$dim[2]
  nrnew <- nr + ntop + nbot
  ncnew <- nc + nleft + nright
  vnew <- cbind(matrix(NA, nr, nleft), x$v, matrix(NA, nr, nright))
  vnew <- rbind(matrix(NA, nbot, ncnew), vnew, matrix(NA, ntop, ncnew))
  # extend x, y coordinate vectors
  xcolnew <- c(if(nleft > 0) x$xcol[1] - (nleft:1) * dx else NULL,
               x$xcol,
               if(nright > 0) x$xcol[nc] + (1:nright) * dx else NULL)
  yrownew <- c(if(nbot > 0) x$yrow[1] - (nbot:1) * dy else NULL,
               x$yrow,
               if(ntop > 0) x$yrow[nr] + (1:ntop) * dy else NULL)
  # rebuild image object
  xnew <- list(v      = vnew,
               dim    = c(nrnew, ncnew),
               xrange = rect$xrange,
               yrange = rect$yrange,
               xstep  = dx,
               ystep  = dy,
               xcol   = xcolnew,
               yrow   = yrownew,
               type   = x$type,
               units   = unitname(x))
  class(xnew) <- "im"
  return(xnew)
}

sort.im <- function(x, ...) {
  verifyclass(x, "im")
  sort(as.vector(as.matrix(x)), ...)
}

dim.im <- function(x) { x$dim }
