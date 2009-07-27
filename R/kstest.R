#
#  kstest.R
#
#  $Revision: 1.31 $  $Date: 2009/07/24 19:43:47 $
#
#

# --------- old -------------

ks.test.ppm <- function(...) {
  .Deprecated("kstest.ppm", package="spatstat")
  kstest.ppm(...)
}

# ---------------------------

kstest <- function(...) {
  UseMethod("kstest")
}

kstest.ppp <-
  function(X, covariate, ..., jitter=TRUE) {
    Xname <- deparse(substitute(X))
    covname <- deparse(substitute(covariate))
    do.call("kstestEngine",
            resolve.defaults(list(ppm(X), covariate, jitter=jitter),
                             list(...),
                             list(modelname="CSR",
                                  covname=covname, dataname=Xname)))
}

kstest.ppm <- function(model, covariate, ..., jitter=TRUE) {
  modelname <- deparse(substitute(model))
  covname <- deparse(substitute(covariate))
  verifyclass(model, "ppm")
  if(is.poisson.ppm(model) && is.stationary.ppm(model))
    modelname <- "CSR"
  do.call("kstestEngine",
          resolve.defaults(list(model, covariate, jitter=jitter),
                           list(...),
                           list(modelname=modelname,
                                covname=covname)))
}

kstestEngine <- function(model, covariate, ...,
                         jitter=TRUE, 
                         modelname, covname, dataname) {
  verifyclass(model, "ppm")
  csr <- is.poisson.ppm(model) && is.stationary.ppm(model)
  
  if(missing(modelname))
    modelname <- if(csr) "CSR" else deparse(substitute(model))
  if(missing(covname))
    covname <- deparse(substitute(covariate))
  if(missing(dataname))
    dataname <- paste(model$call[[2]])

  if(!is.poisson.ppm(model))
    stop("Only implemented for Poisson point process models")

  X <- data.ppm(model)
  W <- as.owin(model)
  
  if(is.im(covariate)) {
    type <- "im"
    # evaluate at data points by interpolation
    ZX <- interp.im(covariate, X$x, X$y)
    # covariate values for pixels inside window
    Z <- covariate[W, drop=FALSE]
    # corresponding mask
    W <- as.owin(Z)
  } else if(is.function(covariate)) {
    type <- "function"
    # evaluate exactly at data points
    ZX <- covariate(X$x, X$y)
    # window
    W <- as.mask(W)
    # covariate in window
    Z <- as.im(covariate, W=W)
    # collapse function body to single string
    if(length(covname) > 1)
      covname <- paste(covname, collapse="")
  } else
     stop("covariate should be an image or a function(x,y)")

  # values of covariate in window
  Zvalues <- as.vector(Z[W, drop=TRUE])
  # corresponding fitted intensity values
  lambda <- as.vector(predict(model, locations=W)[W, drop=TRUE])

  # apply jittering to avoid ties
  if(jitter) {
    dZ <- 0.3 * quantile(diff(sort(unique(c(ZX, Zvalues)))), 1/min(20, X$n))
    ZX <- ZX + rnorm(length(ZX), sd=dZ)
    Zvalues <- Zvalues + rnorm(length(Zvalues), sd=dZ)
  }
  
  # form weighted cdf of Z values in window
  FZ <- ewcdf(Zvalues, lambda/sum(lambda))
  # Ensure support of cdf includes the range of the data
  xxx <- knots(FZ)
  yyy <- FZ(xxx)
  if(min(xxx) > min(ZX)) {
    xxx <- c(min(ZX), xxx)
    yyy <- c(0, yyy)
  }
  if(max(xxx) < max(ZX)) {
    xxx <- c(xxx, max(ZX))
    yyy <- c(yyy, 1)
  }
  # make piecewise linear approximation of cdf
  FZ <- approxfun(xxx, yyy, rule=2)
  # now apply cdf
  U <- FZ(ZX)
  # Test uniformity of transformed values
  result <- ks.test(U, "punif", ...)

  # modify the 'htest' entries
  result$method <- paste("Spatial Kolmogorov-Smirnov test of",
                         if(csr) "CSR" else "inhomogeneous Poisson process")
  result$data.name <-
    paste("covariate", sQuote(paste(covname, collapse="")),
          "evaluated at points of", sQuote(dataname), "\n\t",
          "and transformed to uniform distribution under",
          if(csr) modelname else sQuote(modelname))
  
  # additional class 'kstest'
  class(result) <- c("kstest", class(result))
  attr(result, "prep") <-
    list(Zimage=Z,
         Zvalues=Zvalues,
         lambda=lambda,
         ZX=ZX, FZ=FZ, FZX=ecdf(ZX),
         U=U, type=type)
  attr(result, "info") <- list(modelname=modelname, covname=covname,
                               dataname=dataname, csr=csr)
  return(result)        
}

plot.kstest <- function(x, ..., lwd=par("lwd"), col=par("col"), lty=par("lty"),
                                lwd0=lwd, col0=col, lty0=lty) {
  prep <- attr(x, "prep")
  info <- attr(x, "info")
  FZ <- prep$FZ
  xxx <- get("x", environment(FZ))
  yyy <- get("y", environment(FZ))
  main <- c(x$method,
            paste("based on distribution of covariate", sQuote(info$covname)),
            paste("p-value=", signif(x$p.value, 4)))
  do.call("plot.default",
          resolve.defaults(
                           list(x=xxx, y=yyy, type="l"),
                           list(...),
                           list(lwd=lwd0, col=col0, lty=lty0),
                           list(xlab=info$covname, ylab="probability",
                                main=main)))
  FZX <- prep$FZX
  if(is.null(FZX))
    FZX <- ecdf(prep$ZX)
  plot(FZX, add=TRUE, do.points=FALSE, lwd=lwd, col=col, lty=lty)
  return(invisible(NULL))
}

