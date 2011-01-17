#
#   measures.R
#
#  signed/vector valued measures with atomic and diffuse components
#
#  $Revision: 1.15 $  $Date: 2011/01/15 03:16:06 $
#
msr <- function(qscheme, discrete, continuous, check=TRUE) {
  if(!inherits(qscheme, "quad"))
    stop("qscheme should be a quadrature scheme")
  nquad <- n.quad(qscheme)
  U <- union.quad(qscheme)
  W <- w.quad(qscheme)
  Z <- is.data(qscheme)
  ndata <- sum(Z)
  # ensure conformable vectors/matrices
  if(is.vector(discrete) && is.vector(continuous)) {
    if(check) {
      check.nvector(discrete,   ndata, things="data points", naok=TRUE)
      check.nvector(continuous, nquad, things="quadrature points", naok=TRUE)
    }
    discretepad <- rep(0, nquad)
    discretepad[Z] <- discrete
  } else {
    discrete <- as.matrix(discrete)
    continuous <- as.matrix(continuous)
    if(check) {
      check.nmatrix(discrete, ndata, things="data points",
                    naok=TRUE, squarematrix=FALSE)
      check.nmatrix(continuous, nquad, things="quadrature points",
                    naok=TRUE, squarematrix=FALSE)
    }
    nd <- ncol(discrete)
    nc <- ncol(continuous)
    if(nd != nc) {
      if(nd == 1) {
        discrete <- matrix(rep(discrete, nc), ndata, nc)
      } else if(nc == 1) {
        continuous <- matrix(rep(continuous, nd), nquad, nd)
      } else stop(paste("Incompatible numbers of columns in",
                        sQuote("discrete"), paren(nd), "and",
                        sQuote("continuous"), paren(nc)))
    }
    discretepad <- matrix(0, nquad, max(nd, nc))
    discretepad[Z, ] <- discrete
  }

  #
  #
  # Discretised measure (value of measure for each quadrature tile)
  val <- discretepad + W * continuous
  #
  out <- list(loc = U,
              val = val,
              atoms = Z,
              discrete = discretepad,
              continuous = continuous,
              wt = W)
  class(out) <- "msr"
  return(out)
}

print.msr <- function(x, ...) {
  n <- npoints(x$loc)
  d <- ncol(as.matrix(x$val))
  descrip <- if(d == 1) "Scalar" else paste(d, "dimensional vector", sep="-")
  cat(paste(descrip, "-valued measure\n", sep=""))
  if(d > 1 && !is.null(cn <- colnames(x$val)))
    cat(paste("vector components:", commasep(sQuote(cn)), "\n"))
  cat(paste("Approximated by", n, "quadrature points\n"))
  print(as.owin(x$loc))
  cat(paste(sum(x$atoms), "atoms\n"))
  cat(paste("Total mass:\n"))
  if(d == 1) {
    cat(paste("discrete =", signif(sum(x$discrete), 5),
              "\tcontinuous =", signif(sum(x$wt * x$continuous), 5),
              "\ttotal =", signif(sum(x$val), 5), "\n"))
  } else {
    if(is.null(cn)) cn <- paste("component", 1:d)
    for(j in 1:d) {
      cat(paste(cn[j], ":\t",
                "discrete =", signif(sum(x$discrete[,j]), 5),
                "\tcontinuous =", signif(sum((x$wt * x$continuous)[,j]), 5),
                "\ttotal =", signif(sum(x$val[,j]), 5), "\n"))
    }
  }
  return(invisible(NULL))
}

plot.msr <- function(x, ...) {
  xname <- deparse(substitute(x))
  d <- ncol(as.matrix(x$val))  
  if(d == 1) {
    smo <- smooth.ppp(x$loc %mark% x$continuous, sigma=max(nndist(x$loc)), ...)
    xtra <- unique(c(names(formals(plot.default)),
                     names(formals(image.default))))
    do.call.matched("plot.im",
                    resolve.defaults(list(x=smo),
                                     list(...),
                                     list(main=xname)),
                    extrargs=xtra)
    xtra <- unique(c(names(formals(plot.owin)),
                     names(formals(points)),
                     names(formals(symbols))))
    do.call.matched("plot.ppp",
                    resolve.defaults(list(x=x$loc %mark% x$discrete),
                                     list(add=TRUE),
                                     list(...)),
                    extrargs=xtra)
  } else {
    # split into a list of real-valued measures
    lis <- list()
    for(j in 1:d) 
      lis[[j]] <- x[,j]
    lis <- as.listof(lis)
    if(!is.null(cn <- colnames(x$val)))
      names(lis) <- cn
    do.call("plot.listof", resolve.defaults(list(lis),
                                            list(...),
                                            list(main=xname)))
  }
  return(invisible(NULL))  
}

"[.msr" <- function(x, i, j, ...) {
  valu  <- as.matrix(x$val)
  disc  <- as.matrix(x$discrete)
  cont  <- as.matrix(x$continuous)
  wt    <- x$wt
  atoms <- x$atoms
  #
  if(!missing(j)) {
    valu <- valu[, j]
    disc <- disc[, j]
    cont <- cont[, j]
  }
  loc <- x$loc
  if(!missing(i)) {
    # use [.ppp to identify which points are retained
    locn  <- loc %mark% seq(npoints(loc))
    loci  <- locn[i]
    loc   <- unmark(loci)
    id    <- marks(loci)
    # extract
    valu  <- valu[id, ]
    disc  <- disc[id, ]
    cont  <- cont[id, ]
    wt    <- wt[id]
    atoms <- atoms[id]
  }
  out <- list(loc=loc,
              val=valu,
              atoms=atoms,
              discrete=disc,
              continuous=cont,
              wt=wt)
  class(out) <- "msr"
  return(out)    
}

dim.msr <- function(x) { dim(as.matrix(x$val)) }

dimnames.msr <- function(x) { list(NULL, colnames(x$val)) }

smooth.msr <- function(X, ...) {
  verifyclass(X, "msr")
  loc <- X$loc
  val <- X$val
  d <- ncol(as.matrix(val))
  if(d == 1) {
    result <- density(loc, weights=val, ...)
  } else {
    result <- list()
    for(j in 1:d) 
      result[[j]] <- density(loc, weights=val[,j], ...)
    result <- as.listof(result)
    names(result) <- colnames(X)
  }
  return(result)
}
