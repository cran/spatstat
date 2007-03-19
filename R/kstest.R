#
#  kstest.R
#
#  $Revision: 1.11 $  $Date: 2007/03/19 06:28:09 $
#
#
ks.test.ppm <- function(model, covariate, ...) {
  modelname <- deparse(substitute(model))
  covname <- deparse(substitute(covariate))
  stopifnot(is.ppm(model))
  if(!is.poisson.ppm(model))
    stop("Only implemented for Poisson point process models")
  X <- data.ppm(model)
  W <- as.owin(model)
  if(is.im(covariate)) {
    type <- "im"
    # evaluate at data points by interpolation
    ZX <- interp.im(covariate, X$x, X$y)
    # window
    W <- as.owin(Z)
    # covariate in window
    Z <- covariate[W, drop=FALSE]
  } else if(is.function(covariate)) {
    type <- "function"
    # evaluate exactly at data points
    ZX <- covariate(X$x, X$y)
    # window
    W <- as.mask(W)
    # covariate in window
    Z <- as.im(covariate, W=W)
  } else
     stop("covariate should be an image or a function(x,y)")

  # values of covariate in window
  Zvalues <- as.vector(Z[W, drop=TRUE])
  # corresponding fitted intensity values
  lambda <- as.vector(predict(model, locations=W)[W, drop=TRUE])
  # form weighted cdf of Z values in window
  FZ <- ewcdf(Zvalues, lambda/sum(lambda))
  # smooth out jumps
  xxx <- knots(FZ)
  yyy <- FZ(xxx)
  xxx <- c(min(ZX), xxx, max(ZX))
  yyy <- c(0, yyy, 1)
  if(any(dup <- duplicated(xxx))) {
    xxx <- xxx[!dup]
    yyy <- yyy[!dup]
  }
  FZ <- approxfun(xxx, yyy, rule=2)
  # now apply cdf
  U <- FZ(ZX)
  # Test uniformity of transformed values
  result <- ks.test(U, "punif", ...)
  result$data.name <-
    paste("predicted cdf of covariate", sQuote(paste(covname, collapse="")),
          "evaluated at data points of", sQuote(modelname))
  attr(result, "prep") <-
    list(Zvalues=Zvalues, lambda=lambda, ZX=ZX, U=U, type=type)
  return(result)        
}

