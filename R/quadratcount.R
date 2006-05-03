#
#  quadratcount.R
#
#  $Revision: 1.5 $  $Date: 2006/04/12 09:54:08 $
#

quadratcount <- function(X, nx=5, ny=nx, xbreaks=NULL, ybreaks=NULL)  {
  verifyclass(X, "ppp")

  xr <- X$window$xrange
  yr <- X$window$yrange
  b <- quadrat.breaks(xr, yr, nx, ny, xbreaks, ybreaks)
  Xcount <- quadrat.count.engine(X$x, X$y, b$xbreaks, b$ybreaks)
  return(Xcount)
}

quadrat.test <- function(X, nx=5, ny=nx, xbreaks=NULL, ybreaks=NULL,
                         fit) {
  Xname <- deparse(substitute(X))
  if(!missing(fit)) {
    fitname <- deparse(substitute(fit))
    verifyclass(fit, "ppm")
    if(!is.poisson.ppm(fit))
      stop("Test is only defined for Poisson point process models")
    if(inherits(X, "ppm"))
      stop("I am confused: both X and fit are ppm objects")
  }
  else if(inherits(X, "ppm")) {
    fit <- X
    X <- data.ppm(fit)
    fitname <- Xname
    Xname <- paste("data from", fitname)
  } else 
    fitname <- NULL
  Xcount <- quadratcount(X, nx, ny, xbreaks, ybreaks)
  xbreaks <- attr(Xcount, "xbreaks")
  ybreaks <- attr(Xcount, "ybreaks")
  W <- X$window
  # determine expected values under model
  if(missing(fit)) {
    testname <- "Chi-squared test of CSR using quadrat counts"
    if(W$type == "rectangle") {
      areas <- outer(diff(xbreaks), diff(ybreaks), "*")
      fitmeans <- X$n * areas/sum(areas)
      df <- length(fitmeans) - 1
    } else {
      W <- as.mask(W)
      xx <- as.vector(raster.x(W))[W$m]
      yy <- as.vector(raster.y(W))[W$m]
      areas <- quadrat.count.engine(xx, yy, xbreaks, ybreaks)
      fitmeans <- X$n * areas/sum(areas)
      df <- length(fitmeans) - 1
    }
  } else {
    testname <- paste("Chi-squared test of fitted model",
                      sQuote(fitname),
                      "using quadrat counts")
    Q <- quad.ppm(fit)
    xx <- x.quad(Q)
    yy <- y.quad(Q)
    ww <- w.quad(Q)
    lambda <- fitted(fit)
    masses <- lambda * ww
    fitmeans <- quadrat.count.engine(xx, yy, xbreaks, ybreaks, weights=masses)
    df <- length(fitmeans) - length(coef(fit))
  }
  OBS <- as.vector(Xcount)
  EXP <- as.vector(fitmeans)
  if(df < 1)
    stop(paste("Not enough quadrats: degrees of freedom df =", df))
  if(any(EXP < 5))
    warning("Some expected counts are small; chi^2 approximation may be inaccurate")
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
  return(result)
}

quadrat.breaks <- function(xr, yr, nx=5, ny=nx, xbreaks=NULL, ybreaks=NULL) {
  if(is.null(xbreaks))
    xbreaks <- seq(xr[1], xr[2], length=nx+1)
  else if(min(xbreaks) > xr[1] || max(xbreaks) < xr[2])
    stop("xbreaks do not span the range of x coordinates in the window")
  if(is.null(ybreaks))
    ybreaks <- seq(yr[1], yr[2], length=ny+1)
  else if(min(ybreaks) > yr[1] || max(ybreaks) < yr[2])
    stop("ybreaks do not span the range of y coordinates in the window")
  return(list(xbreaks=xbreaks, ybreaks=ybreaks))
}

quadrat.count.engine <- function(x, y, xbreaks, ybreaks, weights) {
  if(min(x) < min(xbreaks) || max(x) > max(xbreaks))
    stop("xbreaks do not span the actual range of x coordinates in data")
  if(min(y) < min(ybreaks) || max(y) > max(ybreaks))
    stop("ybreaks do not span the actual range of y coordinates in data")
  xg <- cut(x, breaks=xbreaks, include.lowest=TRUE)
  yg <- cut(y, breaks=ybreaks, include.lowest=TRUE)
  if(missing(weights)) 
    sumz <- table(list(x=xg, y=yg))
  else {
    sumz <- tapply(weights, list(x=xg, y=yg), sum)
    if(any(nbg <- is.na(sumz)))
      sumz[nbg] <- 0
  }
  attr(sumz, "xbreaks") <- xbreaks
  attr(sumz, "ybreaks") <- ybreaks
  return(sumz)
}
