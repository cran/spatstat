#
# bermantest.R
#
# Test statistics from Berman (1986)
#
#  $Revision: 1.4 $  $Date: 2009/07/24 21:12:57 $
#
#

# ---------------------------

bermantest <- function(...) {
  UseMethod("bermantest")
}

bermantest.ppp <-
  function(X, covariate,
           which=c("Z1", "Z2"),
           alternative=c("two.sided", "less", "greater"),
           ...) {
    Xname <- deparse(substitute(X))
    covname <- deparse(substitute(covariate))
    which <- match.arg(which)
    alternative <- match.arg(alternative)

    do.call("bermantestEngine",
            resolve.defaults(list(ppm(X), covariate, which, alternative),
                             list(...),
                             list(modelname="CSR",
                                  covname=covname, dataname=Xname)))
}

bermantest.ppm <- function(model, covariate,
                           which=c("Z1", "Z2"),
                           alternative=c("two.sided", "less", "greater"),
                           ...) {
  modelname <- deparse(substitute(model))
  covname <- deparse(substitute(covariate))
  verifyclass(model, "ppm")
  which <- match.arg(which)
  alternative <- match.arg(alternative)
  if(is.poisson.ppm(model) && is.stationary.ppm(model))
    modelname <- "CSR"
  do.call("bermantestEngine",
          resolve.defaults(list(model, covariate, which, alternative),
                           list(...),
                           list(modelname=modelname,
                                covname=covname)))
}

bermantestEngine <- function(model, covariate,
                             which=c("Z1", "Z2"),
                             alternative=c("two.sided", "less", "greater"),
                             ...,
                             modelname, covname, dataname) {

  
  if(missing(modelname))
    modelname <- if(csr) "CSR" else deparse(substitute(model))
  if(missing(covname))
    covname <- deparse(substitute(covariate))
  if(missing(dataname))
    dataname <- paste(model$call[[2]])

  which <- match.arg(which)
  alternative <- match.arg(alternative)

  verifyclass(model, "ppm")
  if(!is.poisson.ppm(model))
    stop("Only implemented for Poisson point process models")
  csr <- is.stationary.ppm(model)

  # data point pattern
  X <- data.ppm(model)
  npoints <- X$n
  
  # ........... first do Kolmogorov-Smirnov ...............
  ksout <- kstestEngine(model, covariate, ...,
                        alternative=alternative,
                        modelname=modelname,
                        covname=covname,
                        dataname=dataname)
  ksprep <- attr(ksout, "prep")
  # covariate image
  Z <- ksprep$Zimage
  # values of covariate at data points
  ZX <- ksprep$ZX
  # transformed to Unif[0,1] under H0
  U  <- ksprep$U
  # domain
  W <- as.owin(Z)
  # intensity of model
  lambda <- predict(model, locations=W)

  switch(which,
         Z1={
           #......... Berman Z1 statistic .....................
           method <- paste("Berman Z1 test of",
                           if(csr) "CSR" else "inhomogeneous Poisson process")
           # sum of covariate values at data points
           Sn <- sum(ZX)
           # predicted mean and variance
           ESn   <- summary(eval.im(lambda * Z  ))$integral
           varSn <- summary(eval.im(lambda * Z^2))$integral
           # standardise
           statistic <- (Sn - ESn)/sqrt(varSn)
           names(statistic) <- "Z1"
           p.value <- switch(alternative,
                            two.sided=2 * pnorm(-abs(statistic)),
                            less=pnorm(statistic),
                            greater=pnorm(statistic, lower.tail=FALSE))
           altblurb <- switch(alternative,
                              two.sided="two-sided",
                              less="mean value of covariate at random points is less than predicted under model",
                              greater="mean value of covariate at random points is greater than predicted under model")
           valuename <- paste("covariate",
                              sQuote(paste(covname, collapse="")),
                              "evaluated at points of",
                              sQuote(dataname))
         },
         Z2={
           #......... Berman Z2 statistic .....................
           method <- paste("Berman Z2 test of",
                           if(csr) "CSR" else "inhomogeneous Poisson process")
           statistic <- sqrt(12/npoints) * (sum(U) - npoints/2)
           names(statistic) <- "Z2"
           p.value <- switch(alternative,
                            two.sided=2 * pnorm(-abs(statistic)),
                            less=pnorm(statistic),
                            greater=pnorm(statistic, lower.tail=FALSE))
           altblurb <- switch(alternative,
                              two.sided="two-sided",
                              less="covariate values at random points have lower quantiles than predicted under model",
                              greater="covariate values at random points have higher quantiles than predicted under model")
           valuename <- paste("covariate",
                              sQuote(paste(covname, collapse="")),
                              "evaluated at points of",
                              sQuote(dataname), "\n\t",
                              "and transformed to uniform distribution under",
                              if(csr) modelname else sQuote(modelname))
         })
           
  out <- list(statistic=statistic,
              p.value=p.value,
              alternative=altblurb,
              method=method,
              data.name=valuename,
              ks=ksout)
  class(out) <- "htest"
  return(out)
}


