#
#   localpcf.R
#
#  $Revision: 1.9 $  $Date: 2011/05/18 08:03:13 $
#
#

localpcf <- function(X, ..., delta=NULL, rmax=NULL, nr=512, stoyan=0.15) {
  if(length(list(...)) > 0)
    warning("Additional arguments ignored")
  m <- localpcfmatrix(X, delta=delta, rmax=rmax, nr=nr, stoyan=stoyan)
  r <- attr(m, "r")
  delta <- attr(m, "delta")
  nX <- npoints(X)
  # border correction
  dbord <- bdist.points(X)
  m[r[row(m)] > dbord[col(m)]] <- NA
  #
  df <- data.frame(m, r=r, theo=rep(1, length(r)))
  icode <- unlist(lapply(seq_len(nX), numalign, nmax=nX))
  nama <- paste("est", icode, sep="")
  desc <- paste("estimate of %s for point", icode)
  labl <- paste("%s[", icode, "](r)", sep="")
  names(df) <- c(nama, "r", "theo")
  desc <- c(desc, "distance argument r", "theoretical Poisson %s")
  labl <- c(labl, "r", "%s[pois](r)")
  # create fv object
  g <- fv(df, "r", substitute(localg(r), NULL),
          "theo", , c(0, max(r)), labl, desc, fname="localg")
  # default is to display them all
  attr(g, "fmla") <- . ~ r
  fvnames(g, ".") <- names(df)[names(df) != "r"]
  unitname(g) <- unitname(X)
  attr(g, "delta") <- delta
  return(g)
}

localpcfmatrix <- function(X, i=seq_len(npoints(X)), ...,
                           delta=NULL, rmax=NULL,
                           nr=512, stoyan=0.15) {
  missi <- missing(i)
  nX <- npoints(X)
  W <- as.owin(X)
  lambda <- nX/area.owin(W)
  if(is.null(delta)) 
    delta <- stoyan/sqrt(lambda)
  if(is.null(rmax)) 
    rmax <- rmax.rule("K", W, lambda)
  # sort points in increasing order of x coordinate
  oX <- order(X$x)
  Xsort <- X[oX]
  idXsort <- (1:nX)[oX]
  if(missi) {
    Y <- X
    oY <- oX
    Ysort   <- Xsort
    idYsort <- idXsort
  } else {
    # i is some kind of index
    Y <- X[i]
    oY <- order(Y$x)
    Ysort <- Y[oY]
    idYsort <- ((1:nX)[i])[oY]
  }
  nY <- npoints(Y)
  force(nr)
  # call C
  zz <- .C("locpcfx",
           nn1 = as.integer(nY),
           x1  = as.double(Ysort$x),
           y1  = as.double(Ysort$y),
           id1 = as.integer(idYsort),
           nn2 = as.integer(nX),
           x2  = as.double(Xsort$x),
           y2  = as.double(Xsort$y),
           id2 = as.integer(idXsort),
           nnr = as.integer(nr),
           rmaxi=as.double(rmax),
           del=as.double(delta),
           pcf=as.double(double(nr * nY)),
           PACKAGE="spatstat")
  out <- matrix(zz$pcf, nr, nY)
  if(!missi) {
    # reorder columns to match original
    out[, oY] <- out
  }
  out <- out/(2 * pi * lambda)
  attr(out, "r") <- seq(from=0, to=rmax, length.out=nr)
  attr(out, "delta") <- delta
  class(out) <- c("localpcfmatrix", class(out))
  return(out)
}

print.localpcfmatrix <- function(x, ...) {
  cat("Matrix of local pair correlation estimates\n")
  nc <- ncol(x)
  nr <- nrow(x)
  cat(paste("pcf estimates for", nc, ngettext(nc, "point", "points"), "\n"))
  rval <- attr(x, "r")
  cat(paste("r values from 0 to", max(rval), "in", nrow(x), "steps\n"))
  return(invisible(NULL))
}

plot.localpcfmatrix <- function(x, ...) {
  xname <- deparse(substitute(x))
  rval <- attr(x, "r")
  do.call("matplot",
          resolve.defaults(list(rval, x),
                           list(...),
                           list(type="l", main=xname,
                                xlab="r", ylab="pair correlation")))
}

"[.localpcfmatrix" <-
  function(x, i, ...) {
    r     <- attr(x, "r")
    delta <- attr(x, "delta")
    class(x) <- "matrix"
    if(missing(i)) {
      x <- x[ , ...]
    } else {
      x <- x[i, ...]
      if(is.matrix(i))
        return(x)
      r <- r[i]
    }
    if(!is.matrix(x))
      x <- matrix(x, nrow=length(r))
    attr(x, "r") <- r
    attr(x, "delta") <- delta
    class(x) <- c("localpcfmatrix", class(x))
    return(x)
}

