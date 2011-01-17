#
#  quadrattest.R
#
#  $Revision: 1.26 $  $Date: 2010/12/13 07:44:31 $
#


quadrat.test <- function(X, ...) {
  UseMethod("quadrat.test")
}

quadrat.test.splitppp <- function(X, ...)
{
  as.listof(lapply(X, quadrat.test.ppp, ...))
}

quadrat.test.ppp <- function(X, nx=5, ny=nx, ...,
                             xbreaks=NULL, ybreaks=NULL, tess=NULL)
{
  Xname <- deparse(substitute(X))
  do.call("quadrat.testEngine",
          resolve.defaults(list(X, nx=nx, ny=ny,
                                xbreaks=xbreaks, ybreaks=ybreaks,
                                tess=tess),
                           list(...), 
                           list(Xname=Xname, fitname="CSR")))
}

quadrat.test.ppm <- function(X, nx=5, ny=nx, ...,
                             xbreaks=NULL, ybreaks=NULL, tess=NULL)
{
  fitname <- deparse(substitute(X))
  dataname <- paste("data from", fitname)
  if(!is.poisson.ppm(X))
    stop("Test is only defined for Poisson point process models")
  if(is.marked(X))
    stop("Sorry, not yet implemented for marked point process models")
  do.call("quadrat.testEngine",
          resolve.defaults(list(data.ppm(X), nx=nx, ny=ny,
                                xbreaks=xbreaks, ybreaks=ybreaks,
                                tess=tess,
                                fit=X),
                           list(...),
                           list(Xname=dataname, fitname=fitname)))
}

quadrat.test.quadratcount <- function(X, ...) {
  if((nxtra <- length(list(...))) > 0) {
    warning(paste(ngettext(nxtra, "argument", "arguments"),
                  "ignored"))
  }
  quadrat.testEngine(Xcount=X)
}


quadrat.testEngine <- function(X, nx, ny, ..., Xcount=NULL,
                               xbreaks=NULL, ybreaks=NULL, tess=NULL,
                               fit=NULL, Xname=NULL, fitname=NULL) {
  if(length(list(...)) > 0) {
    nama <- names(list(...))
    warning(paste(ngettext(length(nama), "argument", "arguments"),
                  paste(dQuote(nama), collapse=", "),
                  "ignored"))
  }
  if(is.null(Xcount))
    Xcount <- quadratcount(X, nx=nx, ny=ny, xbreaks=xbreaks, ybreaks=ybreaks,
                           tess=tess)
  tess <- attr(Xcount, "tess")
  # determine expected values under model
  if(is.null(fit)) {
    testname <- "Chi-squared test of CSR using quadrat counts"
    if(tess$type == "rect") 
      areas <- outer(diff(tess$xgrid), diff(tess$ygrid), "*")
    else 
      areas <- unlist(lapply(tiles(tess), area.owin))
    fitmeans <- sum(Xcount) * areas/sum(areas)
    df <- length(fitmeans) - 1
  } else {
    if(!is.ppm(fit))
      stop("fit should be a ppm object")
    if(!is.poisson.ppm(fit))
      stop("Quadrat test only supported for Poisson point process models")
    if(is.marked(fit))
      stop("Sorry, not yet implemented for marked point process models")
    testname <- paste("Chi-squared test of fitted Poisson model",
                      sQuote(fitname),
                      "using quadrat counts")
    Q <- quad.ppm(fit, drop=TRUE)
    ww <- w.quad(Q)
    lambda <- fitted(fit, drop=TRUE)
    masses <- lambda * ww
    # sum weights of quadrature points in each tile 
    if(tess$type == "rect") {
      xx <- x.quad(Q)
      yy <- y.quad(Q)
      xbreaks <- tess$xgrid
      ybreaks <- tess$ygrid
      fitmeans <- rectquadrat.countEngine(xx, yy, xbreaks, ybreaks,
                                          weights=masses)
      fitmeans <- as.vector(t(fitmeans))
    } else {
      U <- as.ppp(Q)
      V <- marks(cut(U, tess), dfok=FALSE)
      fitmeans <- tapply(masses, list(tile=V), sum)
      fitmeans[is.na(fitmeans)] <- 0
    }
    df <- length(fitmeans) - length(coef(fit))
  }
  OBS <- as.vector(t(as.table(Xcount)))
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
  attr(result, "quadratcount") <- Xcount
  return(result)
}

print.quadrattest <- function(x, ...) {
  NextMethod("print")
  pm <- function(te, ..., brief=TRUE) {
    print(te, ..., brief=brief)
  }
  cat("Quadrats: ")
  pm(as.tess(x), ...)
  return(invisible(NULL))
}

plot.quadrattest <- function(x, ...) {
  xname <- deparse(substitute(x))
  Xcount <- attr(x, "quadratcount")
  # plot tessellation
  tess  <- attr(Xcount, "tess")
  do.call("plot.tess",
          resolve.defaults(list(tess),
                           list(...),
                           list(main=xname)))
  # compute locations for text
  til <- tiles(tess)
  ok <- unlist(lapply(til, function(x) { !is.null(x) && area.owin(x) > 0 }))
  incircles <- lapply(til[ok], incircle)
  x0 <- unlist(lapply(incircles, function(z) { z$x }))
  y0 <- unlist(lapply(incircles, function(z) { z$y }))
  ra <- unlist(lapply(incircles, function(z) { z$r }))
  #
  dotext <- function(dx, dy, values, x0, y0, ra, ...) {
    do.call.matched("text.default",
                    resolve.defaults(list(x=x0 + dx * ra, y = y0 + dy * ra),
                                     list(labels=paste(as.vector(values))),
                                     list(...)))
  }
  # plot observed counts
  cos30 <- sqrt(2)/2
  sin30 <- 1/2
  f <- 0.4
  dotext(-f * cos30, f * sin30, as.vector(t(as.table(Xcount)))[ok],
         x0, y0, ra, adj=c(1,0), ...)
  # plot expected counts
  dotext(f * cos30, f * sin30, round(x$expected,1)[ok], x0, y0, ra, adj=c(0,0),
         ...)
  # plot Pearson residuals
  dotext(0, -f,  signif(x$residuals,2)[ok],x0, y0, ra, ...)

  return(invisible(NULL))
}


