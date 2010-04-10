#
# bermantest.R
#
# Test statistics from Berman (1986)
#
#  $Revision: 1.6 $  $Date: 2010/04/08 10:18:40 $
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
  
  # ........... first assemble data ...............
  fram <- spatialCDFframe(model, covariate, ...,
                        modelname=modelname,
                        covname=covname,
                        dataname=dataname)
  fprep <- fram$prep
  # covariate image
  Z <- fprep$Zimage
  # values of covariate at data points
  ZX <- fprep$ZX
  # transformed to Unif[0,1] under H0
  U  <- fprep$U
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
           slamZ  <- summary(eval.im(lambda * Z))
           slamZ2 <- summary(eval.im(lambda * Z^2))
           ESn   <- slamZ$integral
           varSn <- slamZ2$integral 
           # working, for plot method
           working <- list(meanZX=mean(ZX),
                           meanZ=slamZ$mean)
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
           working <- list(meanU=mean(U))
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
              which=which,
              working=working,
              data.name=valuename,
              fram=fram)
  class(out) <- c("htest", "bermantest")
  return(out)
}

plot.bermantest <-
  function(x, ..., lwd=par("lwd"), col=par("col"), lty=par("lty"),
           lwd0=lwd, col0=col, lty0=lty)
{
  fram <- x$fram
  if(!is.null(fram)) {
    prep <- fram$prep
    info <- fram$info
  } else {
    # old style
    ks <- x$ks
    prep <- attr(ks, "prep")
    info <- attr(ks, "info")
  }
  work <- x$working
  switch(x$which,
         Z1={
           # plot cdf's of Z
           FZ <- prep$FZ
           xxx <- get("x", environment(FZ))
           yyy <- get("y", environment(FZ))
           main <- c(x$method,
                     paste("based on distribution of covariate",
                           sQuote(info$covname)),
                     paste("Z1 statistic =", signif(x$statistic, 4)),
                     paste("p-value=", signif(x$p.value, 4)))
           do.call("plot.default",
                   resolve.defaults(
                                    list(x=xxx, y=yyy, type="l"),
                                    list(...),
                                    list(lwd=lwd0, col=col0, lty=lty0),
                                    list(xlab=info$covname,
                                         ylab="probability",
                                         main=main)))
           FZX <- prep$FZX
           if(is.null(FZX))
             FZX <- ecdf(prep$ZX)
           plot(FZX, add=TRUE, do.points=FALSE, lwd=lwd, col=col, lty=lty)
           abline(v=work$meanZ, lwd=lwd0,col=col0, lty=lty0)
           abline(v=work$meanZX, lwd=lwd,col=col, lty=lty)
         },
         Z2={
           # plot cdf of U
           U <- prep$U
           cdfU <- ecdf(U)
           main <- c(x$method,
                     paste("based on distribution of covariate",
                           sQuote(info$covname)),
                     paste("Z2 statistic =", signif(x$statistic, 4)),
                     paste("p-value=", signif(x$p.value, 4)))
           do.call("plot.ecdf",
                   resolve.defaults(
                                    list(cdfU),
                                    list(...),
                                    list(do.points=FALSE, asp=1),
                                    list(lwd=lwd, col=col, lty=lty),
                                    list(xlab="U", ylab="relative frequency"),
                                    list(main=main)))
           abline(0,1,lwd=lwd0,col=col0,lty=lty0)
           abline(v=0.5, lwd=lwd0,col=col0,lty=lty0)
           abline(v=work$meanU, lwd=lwd,col=col,lty=lty)
         })
  return(invisible(NULL))
}



