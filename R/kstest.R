#
#  kstest.R
#
#  $Revision: 1.35 $  $Date: 2010/03/08 03:50:57 $
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
    if(!is.marked(X, dfok=TRUE)) {
      # unmarked
      model <- ppm(X)
      modelname <- "CSR"
    } else if(is.multitype(X)) {
      # multitype
      model <- ppm(X, ~marks)
      modelname <- "CSRI"
    } else {
      # marked - general case
      X <- unmark(X)
      warning("marks ignored")
      model <- ppm(X)
      modelname <- "CSR"
    }
    do.call("kstestEngine",
            resolve.defaults(list(model, covariate, jitter=jitter),
                             list(...),
                             list(modelname=modelname,
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

  if(!is.marked(model)) {
    # ...................  unmarked .......................
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
  } else {
    # ...................  marked .......................
    if(!is.multitype(model))
      stop("Only implemented for multitype models (factor marks)")
    marx <- marks(X, dfok=FALSE)
    possmarks <- levels(marx)
    npoints <- X$n
    # single image: replicate 
    if(is.im(covariate))
      covariate <- lapply(possmarks, function(x,v){v}, v=covariate)
    #
    if(is.list(covariate) && all(unlist(lapply(covariate, is.im)))) {
      # list of images
      type <- "im"
      if(length(covariate) != length(possmarks))
        stop("Number of images does not match number of possible marks")
      # evaluate covariate at each data point by interpolation
      ZX <- numeric(npoints)
      for(k in seq(possmarks)) {
        ii <- (marx == possmarks[k])
        ZX[ii] <- interp.im(covariate[[k]], x=X$x[ii], y=X$y[ii])
      }
      # restrict covariate images to window 
      Z <- lapply(covariate, function(x,W){x[W, drop=FALSE]}, W=W)
      # extract pixel locations and pixel values
      Zframes <- lapply(Z, as.data.frame)
      # covariate values at each pixel inside window
      Zvalues <- unlist(lapply(Zframes, function(df) { df[ , 3] }))
      # pixel locations 
      locn <- lapply(Zframes, function(df) { df[ , 1:2] })
      # tack on mark values
      for(k in seq(possmarks))
        locn[[k]] <- cbind(locn[[k]], data.frame(marks=possmarks[k]))
      loc <- do.call("rbind", locn)
      # corresponding fitted intensity values
      lambda <- predict(model, locations=loc)
    } else if(is.function(covariate)) {
      type <- "function"
      # evaluate exactly at data points
      ZX <- covariate(X$x, X$y, marx)
      # same window
      W <- as.mask(W)
      # covariate in window
      Z <- list()
      g <- function(x,y,m,f) { f(x,y,m) }
      for(k in seq(possmarks))
        Z[[k]] <- as.im(g, m=possmarks[k], f=covariate, W=W)
      Zvalues <- unlist(lapply(Z, function(z) { as.data.frame(z)[,3] }))
      # corresponding fitted intensity values
      lambda <- predict(model, locations=W)
      lambda <- unlist(lapply(lambda, function(z) { as.data.frame(z)[,3] }))
      if(length(lambda) != length(Zvalues))
        stop("Internal error: length(lambda) != length(Zvalues)")
      # collapse function body to single string
      if(length(covname) > 1)
        covname <- paste(covname, collapse="")
    } else
    stop("covariate should be an image or a function(x,y)")
  }    
  # ..........................................................

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

