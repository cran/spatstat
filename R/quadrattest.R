#
#  quadrattest.R
#
#  $Revision: 1.9 $  $Date: 2008/07/25 19:48:53 $
#


quadrat.test <- function(X, ...) {
  UseMethod("quadrat.test")
}

quadrat.test.ppp <- function(X, nx=5, ny=nx, xbreaks=NULL, ybreaks=NULL, ...)
{
  Xname <- deparse(substitute(X))
  do.call("quadrat.testEngine",
          resolve.defaults(list(X, nx=nx, ny=ny,
                                xbreaks=xbreaks, ybreaks=ybreaks),
                           list(...), 
                           list(Xname=Xname, fitname="CSR")))
}

quadrat.test.ppm <- function(X, nx=5, ny=nx, xbreaks=NULL, ybreaks=NULL, ...)
{
  fitname <- deparse(substitute(X))
  dataname <- paste("data from", fitname)
  if(!is.poisson.ppm(X))
    stop("Test is only defined for Poisson point process models")
  do.call("quadrat.testEngine",
          resolve.defaults(list(data.ppm(X), nx=nx, ny=ny,
                                xbreaks=xbreaks, ybreaks=ybreaks, 
                                fit=X),
                           list(...),
                           list(Xname=dataname, fitname=fitname)))
}

quadrat.testEngine <- function(X, nx, ny, xbreaks, ybreaks,
                                ..., fit=NULL, Xname=NULL, fitname=NULL) {
  if(length(list(...)) > 0) {
    nama <- names(list(...))
    warning(paste(ngettext(length(nama), "argument", "arguments"),
                  paste(dQuote(nama), collapse=", "),
                  "ignored"))
  }
  Xcount <- quadratcount(X, nx, ny, xbreaks, ybreaks)
  xbreaks <- attr(Xcount, "xbreaks")
  ybreaks <- attr(Xcount, "ybreaks")
  W <- X$window
  # determine expected values under model
  if(is.null(fit)) {
    testname <- "Chi-squared test of CSR using quadrat counts"
    if(W$type == "rectangle") {
      areas <- outer(diff(xbreaks), diff(ybreaks), "*")
      fitmeans <- X$n * areas/sum(areas)
      df <- length(fitmeans) - 1
    } else {
      W <- as.mask(W)
      xx <- as.vector(raster.x(W))[W$m]
      yy <- as.vector(raster.y(W))[W$m]
      areas <- quadrat.countEngine(xx, yy, xbreaks, ybreaks)
      fitmeans <- X$n * areas/sum(areas)
      df <- length(fitmeans) - 1
    }
  } else {
    testname <- paste("Chi-squared test of fitted model",
                      sQuote(fitname),
                      "using quadrat counts")
    Q <- quad.ppm(fit, drop=TRUE)
    xx <- x.quad(Q)
    yy <- y.quad(Q)
    ww <- w.quad(Q)
    lambda <- fitted(fit, drop=TRUE)
    masses <- lambda * ww
    fitmeans <- quadrat.countEngine(xx, yy, xbreaks, ybreaks, weights=masses)
    df <- length(fitmeans) - length(coef(fit))
  }
  OBS <- as.vector(Xcount)
  EXP <- as.vector(fitmeans)
  if(df < 1)
    stop(paste("Not enough quadrats: degrees of freedom df =", df))
  if(any(EXP < 5))
    warning(paste("Some expected counts are small;",
                  "chi^2 approximation may be inaccurate"),
            call.=FALSE)
  X2 <- sum((OBS - EXP)^2/EXP)
  PVAL <- pchisq(X2, df, lower=FALSE)
  names(X2) <- "X-squared"
  names(df) <- "df"
  result <- structure(list(statistic = X2,
                           parameter = df,
                           p.value = PVAL,
                           method = testname,
                           data.name = Xname,
                           observed = OBS,
                           expected = EXP,
                           residuals = (OBS - EXP)/sqrt(EXP)),
                      class = "htest")
  class(result) <- c("quadrattest", class(result))
  attr(result, "quadrats") <- Xcount
  return(result)
}


plot.quadrattest <- function(x, ...) {
  xname <- deparse(substitute(x))
  quads <- attr(x, "quadrats")
  if(inherits(quads, "quadratcount")) {
    # plot observed counts
    do.call("plot.quadratcount",
            resolve.defaults(list(quads, dx=-0.3, dy=0.3),
                             list(...),
                             list(main=xname)))
    # plot expected counts
    do.call("plot.quadratcount",
            resolve.defaults(list(quads, dx=0.3, dy=0.3),
                             list(add=TRUE),
                             list(entries=signif(x$expected,2)),
                             list(...)))
    # plot Pearson residuals
    do.call("plot.quadratcount",
            resolve.defaults(list(quads, dx=0, dy=-0.3),
                             list(add=TRUE),
                             list(entries=signif(x$residuals,2)),
                             list(...)))
    return(invisible(NULL))
  }
  components <- attr(x, "components")
  do.call("plot",
          resolve.defaults(list(components),
                           list(...),
                           list(main=xname)))
}


