
ks.test.ppm <- function(model, covariate, ...) {
  modelname <- deparse(substitute(model))
  covname <- deparse(substitute(covariate))
  stopifnot(is.ppm(model))
  if(!is.poisson.ppm(model))
    stop("Only implemented for Poisson point process models")
  X <- data.ppm(model)
  W <- as.owin(model)
  if(is.im(covariate)) {
    Z <- covariate
    Z <- Z[W, drop=FALSE]
    lambda <- predict(model, locations=as.owin(Z))
    Zvalues <- as.matrix(Z)
    lambda <- as.vector(as.matrix(lambda))
    FZ <- ewcdf(Zvalues, lambda/sum(lambda))
    U <- FZ(Z[X])
  } else if(is.function(covariate)) {
    W <- as.mask(W)
    lambda <- predict(model, locations=W)
    lambda <- as.vector(as.matrix(lambda))
    xx <- as.vector(raster.x(W)[W$m])
    yy <- as.vector(raster.y(W)[W$m])
    Zvalues <- covariate(xx, yy)
    FZ <- ewcdf(Zvalues, lambda/sum(lambda))
    # smooth out jumps
    xxx <- knots(FZ)
    yyy <- FZ(xxx)
    FZ <- approxfun(xxx, yyy, rule=2)
    # evaluate at data points
    ZX <- covariate(X$x, X$y)
    U <- FZ(ZX)
  } else
     stop("covariate should be an image or a function(x,y)")
  result <- ks.test(U, "punif", ...)
  result$data.name <-
    paste("predicted cdf of covariate", sQuote(paste(covname, collapse="")),
          "evaluated at data points of", sQuote(modelname))
  return(result)        
}
