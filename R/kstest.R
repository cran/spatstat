#
#  kstest.R
#
#  $Revision: 1.46 $  $Date: 2010/06/08 10:02:40 $
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
    covname <- singlestring(deparse(substitute(covariate)))
    if(is.character(covariate)) covname <- covariate
    if(!is.marked(X, dfok=TRUE)) {
      # unmarked
      model <- ppm(X)
      modelname <- "CSR"
    } else if(is.multitype(X)) {
      # multitype
      mf <- summary(X)$marks$frequency
      if(all(mf > 0)) {
        model <- ppm(X, ~marks)
        modelname <- "CSRI"
      } else {
        warning("Ignoring marks, because some mark values have zero frequency")
        X <- unmark(X)
        model <- ppm(X)
        modelname <- "CSR"
      } 
    } else {
      # marked - general case
      X <- unmark(X)
      warning("marks ignored")
      model <- ppm(X)
      modelname <- "CSR"
    }
    do.call("spatialCDFtest",
            resolve.defaults(list(model, covariate, test="ks"),
                             list(jitter=jitter),
                             list(...),
                             list(modelname=modelname,
                                  covname=covname, dataname=Xname)))
}

kstest.ppm <- function(model, covariate, ..., jitter=TRUE) {
  modelname <- deparse(substitute(model))
  covname <- singlestring(deparse(substitute(covariate)))
  if(is.character(covariate)) covname <- covariate
  verifyclass(model, "ppm")
  if(is.poisson.ppm(model) && is.stationary.ppm(model))
    modelname <- "CSR"
  do.call("spatialCDFtest",
          resolve.defaults(list(model, covariate, test="ks"),
                           list(jitter=jitter),
                           list(...),
                           list(modelname=modelname,
                                covname=covname)))
}

spatialCDFtest <- function(model, covariate, test, ...,
                           jitter=TRUE, 
                           modelname=NULL, covname=NULL, dataname=NULL) {
  if(!is.poisson.ppm(model))
    stop("Only implemented for Poisson point process models")
  # conduct test based on comparison of CDF's of covariate values
  test <- pickoption("test", test, c(ks="ks"))
  # compute the essential data
  fra <- spatialCDFframe(model, covariate,
                         jitter, modelname, covname, dataname)
  values <- fra$values
  info   <- fra$info
  # Test uniformity of transformed values
  U <- values$U
  switch(test,
         ks={ result <- ks.test(U, "punif", ...) },
         stop("Unrecognised test option"))

  # modify the 'htest' entries
  csr <- info$csr
  result$method <- paste("Spatial Kolmogorov-Smirnov test of",
                         if(csr) "CSR" else "inhomogeneous Poisson process")
  result$data.name <-
    paste("covariate", sQuote(singlestring(info$covname)),
          "evaluated at points of", sQuote(info$dataname), "\n\t",
          "and transformed to uniform distribution under",
          if(csr) info$modelname else sQuote(info$modelname))
  
  # additional class 'kstest'
  class(result) <- c("kstest", class(result))
  attr(result, "frame") <- fra
  return(result)        
}

spatialCDFframe <- function(model, covariate, ...) {
  # evaluate CDF of covariate values at data points and at pixels
  stuff <- evalCovar(model, covariate, ...)
  # extract 
  values <- stuff$values
  info   <- stuff$info
  Zimage  <- values$Zimage
  Zvalues <- values$Zvalues
  lambda  <- values$lambda
  ZX      <- values$ZX
  type    <- values$type
  # compute empirical cdf of Z values at points of X
  FZX <- ecdf(ZX)
  # form weighted cdf of Z values in window
  FZ <- ewcdf(Zvalues, lambda/sum(lambda))
  # Ensure support of cdf includes the range of the data
  xxx <- knots(FZ)
  yyy <- FZ(xxx)
  minZX <- min(ZX, na.rm=TRUE)
  minxxx <- min(xxx, na.rm=TRUE)
  if(minxxx > minZX) {
    xxx <- c(minZX, xxx)
    yyy <- c(0, yyy)
  }
  maxZX <- max(ZX, na.rm=TRUE)
  maxxxx <- max(xxx, na.rm=TRUE)
  if(maxxxx < maxZX) {
    xxx <- c(xxx, maxZX)
    yyy <- c(yyy, 1)
  }
  # make piecewise linear approximation of cdf
  FZ <- approxfun(xxx, yyy, rule=2)
  # now apply cdf
  U <- FZ(ZX)

  # pack up
  stuff$values$FZ  <- FZ
  stuff$values$FZX <- FZX
  stuff$values$U   <- U
  class(stuff) <- "spatialCDFframe"
  return(stuff)
}

evalCovar <- function(model, covariate, 
                      jitter=TRUE, 
                      modelname=NULL, covname=NULL,
                      dataname=NULL) {
  # evaluate covariate values at data points and at pixels
  verifyclass(model, "ppm")
  csr <- is.poisson.ppm(model) && is.stationary.ppm(model)

  # determine names
  if(is.null(modelname))
    modelname <- if(csr) "CSR" else deparse(substitute(model))
  if(is.null(covname)) {
    covname <- singlestring(deparse(substitute(covariate)))
    if(is.character(covariate)) covname <- covariate
  }
  if(is.null(dataname))
    dataname <- paste(model$call[[2]])
  info <-  list(modelname=modelname, covname=covname,
                dataname=dataname, csr=csr)

  
  # evaluate covariate 
  X <- data.ppm(model)
  W <- as.owin(model)

  if(is.character(covariate)) {
    # One of the characters 'x' or 'y'
    # Turn it into a function.
    ns <- length(covariate)
    if(ns == 0) stop("covariate is empty")
    if(ns > 1) stop("more than one covariate specified")
    covname <- covariate
    covariate <- switch(covariate,
                     x=function(x,y,m){x},
                     y=function(x,y,m){y},
                     stop(paste("Unrecognised covariate", dQuote(covariate))))
  } 
  
  if(!is.marked(model)) {
    # ...................  unmarked .......................
    if(is.im(covariate)) {
      type <- "im"
      # evaluate at data points by interpolation
      ZX <- interp.im(covariate, X$x, X$y)
      # fix boundary glitches
      if(any(uhoh <- is.na(ZX)))
        ZX[uhoh] <- safelookup(covariate, X[uhoh])
      # covariate values for pixels inside window
      Z <- covariate[W, drop=FALSE]
      # corresponding mask
      W <- as.owin(Z)
    } else if(is.function(covariate)) {
      type <- "function"
      # evaluate exactly at data points
      ZX <- covariate(X$x, X$y)
      if(!all(is.finite(ZX)))
        warning("covariate function returned NA or Inf values")
      # window
      W <- as.mask(W)
      # covariate in window
      Z <- as.im(covariate, W=W)
      # collapse function body to single string
      covname <- singlestring(covname)
    } else stop(paste("The covariate should be",
                      "an image, a function(x,y)",
                      "or one of the characters",
                      sQuote("x"), "or", sQuote("y")))
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
        covariate.k <- covariate[[k]]
        values <- interp.im(covariate.k, x=X$x[ii], y=X$y[ii])
        # fix boundary glitches
        if(any(uhoh <- is.na(values)))
          values[uhoh] <- safelookup(covariate.k, X[ii][uhoh])
        ZX[ii] <- values
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
      covname <- singlestring(covname)
    } else stop(paste("For a multitype point process model,",
                      "the covariate should be an image, a list of images,",
                      "a function(x,y,m)", 
                      "or one of the characters",
                      sQuote("x"), "or", sQuote("y")))
  }    
  # ..........................................................

  # apply jittering to avoid ties
  if(jitter) {
    nX <- length(ZX)
    dZ <- 0.3 * quantile(diff(sort(unique(c(ZX, Zvalues)))), 1/min(20, nX))
    ZX <- ZX + rnorm(nX, sd=dZ)
    Zvalues <- Zvalues + rnorm(length(Zvalues), sd=dZ)
  }

  # wrap up 
  values <- list(Zimage=Z,
                 Zvalues=Zvalues,
                 lambda=lambda,
                 ZX=ZX, type=type)
  return(list(values=values, info=info))
}

plot.kstest <- function(x, ..., lwd=par("lwd"), col=par("col"), lty=par("lty"),
                                lwd0=lwd, col0=col, lty0=lty) {
  fram <- attr(x, "frame")
  if(!is.null(fram)) {
    values <- fram$values
    info <- fram$info
  } else {
    # old style
    values <- attr(x, "prep")
    info <- attr(x, "info")
  }
  FZ <- values$FZ
  xxx <- get("x", environment(FZ))
  yyy <- get("y", environment(FZ))
  covname <- info$covname
  covdescrip <- switch(covname,
                       x="x coordinate",
                       y="y coordinate",
                       paste("covariate", dQuote(covname)))
  main <- c(x$method,
            paste("based on distribution of", covdescrip),
            paste("p-value=", signif(x$p.value, 4)))
  do.call("plot.default",
          resolve.defaults(
                           list(x=xxx, y=yyy, type="l"),
                           list(...),
                           list(lwd=lwd0, col=col0, lty=lty0),
                           list(xlab=info$covname, ylab="probability",
                                main=main)))
  FZX <- values$FZX
  if(is.null(FZX))
    FZX <- ecdf(values$ZX)
  plot(FZX, add=TRUE, do.points=FALSE, lwd=lwd, col=col, lty=lty)
  return(invisible(NULL))
}

