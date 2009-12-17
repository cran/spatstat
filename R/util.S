#
#    util.S    miscellaneous utilities
#
#    $Revision: 1.50 $    $Date: 2009/12/15 23:57:21 $
#
#  (a) for matrices only:
#
#    matrowany(X) is equivalent to apply(X, 1, any)
#    matrowall(X) "   "  " "  "  " apply(X, 1, all)
#    matcolany(X) "   "  " "  "  " apply(X, 2, any)
#    matcolall(X) "   "  " "  "  " apply(X, 2, all)
#
#  (b) for 3D arrays only:
#    apply23sum(X)  "  "   "  " apply(X, c(2,3), sum)
#
#  (c) weighted histogram
#    whist()
#
#  (d) for matrices:
#    matrixsample()
#           subsamples or supersamples a matrix
#
#
matrowsum <- function(x) {
  x %*% rep(1, ncol(x))
}

matcolsum <- function(x) {
  rep(1, nrow(x)) %*% x
}
  
matrowany <- function(x) {
  (matrowsum(x) > 0)
}

matrowall <- function(x) {
  (matrowsum(x) == ncol(x))
}

matcolany <- function(x) {
  (matcolsum(x) > 0)
}

matcolall <- function(x) {
  (matcolsum(x) == nrow(x))
}

########
    # hm, this is SLOWER

apply23sum <- function(x) {
  dimx <- dim(x)
  if(length(dimx) != 3)
    stop("x is not a 3D array")
  result <- array(0, dimx[-1])

  nz <- dimx[3]
  for(k in 1:nz) {
    result[,k] <- matcolsum(x[,,k])
  }
  result
}
    
#######################
#
#    whist      weighted histogram
#

whist <- function(x, breaks, weights=NULL) {
    N <- length(breaks)
    if(length(x) == 0) 
      h <- rep(0, N+1)
    else {
      # classify data into histogram cells (breaks need not span range of data)
      cell <- findInterval(x, breaks, rightmost.closed=TRUE)
      cell <- factor(cell, levels=0:N)
      # compute weighted histogram
      if(is.null(weights))
        h <- table(cell)
      else 
        h <- unlist(lapply(split(weights, cell), sum, na.rm=TRUE))
    }
    h <- as.numeric(h)
    y <- h[2:N]
    attr(y, "low") <- h[1]
    attr(y, "high") <- h[N+1]
    return(y)
}

######################
#
#   matrixsample         subsample or supersample a matrix
#

matrixsample <- function(mat, newdim, phase=c(0,0)) {
  olddim <- dim(mat)
  oldlength <- prod(olddim)
  if(all(olddim >= newdim) && all(olddim %% newdim == 0)) {
    # new matrix is periodic subsample of old matrix
    ratio <- olddim/newdim
    phase <- pmin(pmax(phase, 0), ratio-1)
    ii <- seq(1+phase[1], olddim[1], by=ratio[1])
    jj <- seq(1+phase[2], olddim[2], by=ratio[2])
    return(mat[ii,jj])
  } else if(all(olddim <= newdim) && all(newdim %% olddim == 0)) {
    # new matrix is repetition of old matrix
    ratio <- newdim/olddim
    replicate.rows <- function(m, nrep, l=length(m)) {
      matrix(rep(m, rep(nrep, l)), nrow(m) * nrep, ncol(m))
    }
    mrow <- replicate.rows(mat, ratio[1], oldlength)
    mfinal <- t(replicate.rows(t(mrow), ratio[2], oldlength * ratio[1]))
    return(mfinal)
  } else {
    # general case
    newmat <- matrix(, newdim[1], newdim[2])
    newrow <- row(newmat)
    newcol <- col(newmat)
    oldrow <- phase[1] + ceiling((newrow * olddim[1])/newdim[1])
    oldcol <- phase[2] + ceiling((newcol * olddim[2])/newdim[2])
    oldrow[oldrow < 1] <- 1
    oldrow[oldrow > olddim[1]] <- olddim[1]
    oldcol[oldcol < 1] <- 1
    oldcol[oldcol > olddim[2]] <- olddim[2]
    newmat <- matrix(mat[cbind(as.vector(oldrow), as.vector(oldcol))],
                     newdim[1], newdim[2])
    return(newmat)
  }
}

pointgrid <- function(W, ngrid) {
  W <- as.owin(W)
  masque <- as.mask(W, dimyx=ngrid)
  xx <- raster.x(masque)
  yy <- raster.y(masque)
  xx <- xx[masque$m]
  yy <- yy[masque$m]
  return(ppp(xx, yy, W))
}

   
commasep <- function(x) {
  px <- paste(x)
  nx <- length(px)
  if(nx <= 1) return(px)
  commas <- c(rep(", ", length(px)-2),
              " and ", "")
  return(paste(paste(px, commas, sep=""), collapse=""))
}

paren <- function(x) { paste("(", x, ")", sep="") }

ordinal <- function(k) {
  last <- abs(k) %% 10
  lasttwo <- abs(k) %% 100
  isteen <- (lasttwo > 10 & lasttwo < 20)
  ending <- ifelse(isteen, "th",
                   ifelse(last == 1, "st",
                          ifelse(last == 2, "nd",
                                 ifelse(last == 3, "rd",
                                        "th"))))
  return(paste(k, ending, sep=""))
}

# equivalent to rev(cumsum(rev(x)))

revcumsum <- function(x) {
  n <- length(x)
  if(identical(storage.mode(x), "integer")) {
    z <- .C("irevcumsum",
            x=as.integer(x),
            as.integer(n),
            PACKAGE="spatstat")
    return(z$x)
  } else {
    z <- .C("drevcumsum",
            x=as.double(x),
            as.integer(n),
            PACKAGE="spatstat")
    return(z$x)
  }
}

prolongseq <- function(x, newrange) {
  stopifnot(length(newrange) == 2 && newrange[1] < newrange[2])
  stopifnot(length(x) >= 2)
  dx <- diff(x)
  if(any(dx <= 0))
    stop("x must be an increasing sequence")
  if(diff(range(dx)) > 0.01 * abs(mean(dx)))
    stop("x must be evenly spaced")
  dx <- mean(dx)

  # add or trim data to left
  if(x[1] > newrange[1]) {
    leftbit <- seq(x[1], newrange[1], by= -dx) 
    x <- c(rev(leftbit), x[-1])
  } else 
    x <- x[x >= newrange[1]]

  # add or trim data to right
  nx <- length(x)
  if(newrange[2] > x[nx]) {
    rightbit <- seq(x[nx], newrange[2], by= dx)
    x <- c(x[-nx], rightbit)
  } else 
    x <- x[x <= newrange[2]]

  return(x)
}

intersect.ranges <- function(a, b, fatal=TRUE) {
  lo <- max(a[1],b[1])
  hi <- min(a[2],b[2])
  if(lo >= hi) {
    if(fatal) stop("Intersection is empty")
    else return(NULL)
  }
  return(c(lo, hi))
}

.Spatstat.ProgressBar <- NULL

progressreport <- function(i, n, every=max(1, ceiling(n/100)),
                           nperline=min(charsperline,
                             every * ceiling(charsperline /(every+3))),
                           charsperline=60,
                           style=spatstat.options("progress")) {
  switch(style,
         txtbar={
           if(i == 1)
             .Spatstat.ProgressBar <<- txtProgressBar(1, n, 1, style=3)
           setTxtProgressBar(.Spatstat.ProgressBar, i)
           if(i == n)
             close(.Spatstat.ProgressBar)
         },
         tty={
           if(i == n) 
             cat(paste(n, ".\n", sep=""))
           else if(every == 1 || i <= 3)
             cat(paste(i, ",", if(i %% nperline == 0) "\n" else " ", sep=""))
           else {
             if(i %% every == 0) 
               cat(i)
             else
               cat(".")
             if(i %% nperline == 0)
               cat("\n")
           }
           flush.console()
         },
         stop(paste("Unrecognised option for style:", dQuote(style)))
         )
  return(invisible(NULL))
}
  
numalign <- function(i, nmax, zero="0") {
  stopifnot(i <= nmax)
  nplaces <- as.integer(ceiling(log10(nmax+1)))
  out <- blank <- paste(rep(zero, nplaces), collapse="")
  istring <- paste(i)
  ilen <- nchar(istring)
  substr(out, nplaces-ilen+1, nplaces) <- istring
  return(out)
}

ensure2vector <- function(x) {
  xname <- deparse(substitute(x))
  if(!is.numeric(x))
    stop(paste(xname, "is not numeric"))
  n <- length(x)
  if(n == 0 || n > 2)
    stop(paste(xname, "should be of length 1 or 2"))
  if(n == 1)
    return(rep(x,2))
  return(x)
}

check.nvector <- function(v, npoints, fatal=TRUE) {
  # vector of numeric values for each point
  vname <- sQuote(deparse(substitute(v)))
  whinge <- NULL
  if(!is.numeric(v))
    whinge <- paste(vname, "is not numeric")
  else if(length(v) != npoints)
    whinge <- paste("The length of", vname,
                    "should equal the number of data points")
  else if(any(is.na(v)))
    whinge <- paste("Some values of", vname, "are NA or NaN")
  #
  if(!is.null(whinge)) {
    if(fatal) stop(whinge) else return(FALSE)
  }
  return(TRUE)
}

check.nmatrix <- function(m, npoints, fatal=TRUE) {
  # matrix of values for each pair of points
  mname <- sQuote(deparse(substitute(m)))
  whinge <- NULL
  if(!is.matrix(m))
    whinge <- paste(mname, "should be a matrix")
  else if(nrow(m) != ncol(m))
    whinge <- paste(mname, "should be a square matrix")
  else if(nrow(m) != npoints)
    whinge <- paste("Dimensions of", mname,
               "do not match number of data points")
  #
  if(!is.null(whinge)) {
    if(fatal) stop(whinge) else return(FALSE)
  }
  return(TRUE)
}

check.named.vector <- function(x, nam, context="", namopt=character(0)) {
  xtitle <- deparse(substitute(x))
  check.named.thing(x, nam, namopt, sQuote(xtitle),
                    is.numeric(x), "vector", context)
  opt <- namopt %in% names(x)
  return(x[c(nam, namopt[opt])])
}

check.named.list <- function(x, nam, context="", namopt=character(0)) {
  xtitle <- deparse(substitute(x))
  check.named.thing(x, nam, namopt, sQuote(xtitle),
                    is.list(x), "list", context)  
  opt <- namopt %in% names(x)
  return(x[c(nam, namopt[opt])])
}

check.named.thing <- function(x, nam, namopt=character(0), xtitle=NULL,
                              valid=TRUE, type="object", context="",
                              fatal=TRUE) {
  if(is.null(xtitle))
    xtitle <- sQuote(deparse(substitute(x)))
  # check whether names(x) contains all obligatory names 'nam'
  # and possibly some of the optional names 'namopt'
  namesx <- names(x)
  omitted <- !(nam %in% namesx)
  foreign <- !(namesx %in% c(nam, namopt))
  if(valid && !any(omitted) && !any(foreign))
    return(character(0))
  # some condition violated
  if(nzchar(context))
    xtitle <- paste(context, xtitle)
  whinge <- paste(xtitle,
                  "must be a named", paste(type, ",", sep=""),
                  "with components", commasep(nam))
  if(length(namopt) > 0)
    whinge <- paste(whinge, paren(paste("and optionally", commasep(namopt))))
  if(any(omitted)) {
    grizzle <- paste(ngettext(sum(omitted), "parameter", "parameters"),
                     commasep(nam[omitted]),
                     "omitted")
    whinge <- paste(whinge, grizzle, sep="; ")
  }
  if(any(foreign)) {
    grizzle <- paste(ngettext(sum(foreign), "component", "components"),
                     commasep(namesx[foreign]),
                     "not recognised")
    whinge <- paste(whinge, grizzle, sep="; ")
  }
  if(fatal)
    stop(whinge, call.=FALSE)
  return(whinge)
}

evenly.spaced <- function(x, tol=1e-07) {
  # test whether x is evenly spaced and increasing
  dx <- diff(x)
  if(any(dx <= .Machine$double.eps))
    return(FALSE)
  # The following test for equal spacing is used in hist.default
  if(diff(range(dx)) > tol * mean(dx))
    return(FALSE)
  return(TRUE)
}

adjustthinrange <- function(ur,vstep,vr) {
  if(diff(ur) >= vstep) return(ur)
  ur <- mean(ur) + c(-1,1) * vstep/2
  if(ur[1] < vr[1]) ur <- vr[1] + c(0,1)*vstep
  if(ur[2] > vr[2]) ur <- vr[2] - c(1,0)*vstep
  return(ur)
}
