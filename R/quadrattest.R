#
#   quadrattest.R
#
#   $Revision: 1.31 $  $Date: 2012/08/21 02:36:52 $
#

quadrat.test <- function(X, ...) {
   UseMethod("quadrat.test")
}

quadrat.test.splitppp <- function(X, ...)
{
   as.listof(lapply(X, quadrat.test.ppp, ...))
}

quadrat.test.ppp <-
  function(X, nx=5, ny=nx,
           alternative = c("two.sided", "regular", "clustered"),
           method = c("Chisq", "MonteCarlo"),
           conditional=TRUE,
           ..., 
           xbreaks=NULL, ybreaks=NULL,
           tess=NULL, nsim=1999)
{
   Xname <- short.deparse(substitute(X))
   method <- match.arg(method)
   alternative <- match.arg(alternative)
   do.call("quadrat.testEngine",
          resolve.defaults(list(X, nx=nx, ny=ny,
                                alternative=alternative,
                                method=method,
                                conditional=conditional,
                                xbreaks=xbreaks, ybreaks=ybreaks,
                                tess=tess,
                                nsim=nsim),
                           list(...), 
                           list(Xname=Xname, fitname="CSR")))
}

quadrat.test.ppm <-
  function(X, nx=5, ny=nx,
           alternative = c("two.sided", "regular", "clustered"),      
           method=c("Chisq", "MonteCarlo"),
           conditional=TRUE, ...,
           xbreaks=NULL, ybreaks=NULL,
           tess=NULL, nsim=1999)
{
   fitname <- short.deparse(substitute(X))
   dataname <- paste("data from", fitname)
   method <- match.arg(method)
   alternative <- match.arg(alternative)
   if(!is.poisson.ppm(X))
    stop("Test is only defined for Poisson point process models")
   if(is.marked(X))
    stop("Sorry, not yet implemented for marked point process models")
   do.call("quadrat.testEngine",
          resolve.defaults(list(data.ppm(X), nx=nx, ny=ny,
                                alternative=alternative,
                                method=method,
                                conditional=conditional,
                                xbreaks=xbreaks, ybreaks=ybreaks,
                                tess=tess,
                                nsim=nsim, 
                                fit=X),
                           list(...),
                           list(Xname=dataname, fitname=fitname)))
}

quadrat.test.quadratcount <-
  function(X,
           alternative = c("two.sided", "regular", "clustered"),
           method=c("Chisq", "MonteCarlo"),
           conditional=TRUE,
           ...,
           nsim=1999) {
   trap.extra.arguments(...)
   method <- match.arg(method)
   alternative <- match.arg(alternative)
   quadrat.testEngine(Xcount=X,
                      alternative=alternative,
                      method=method, conditional=conditional, nsim=nsim)
}

quadrat.testEngine <- local({
  
  quadrat.testEngine <- function(X, nx, ny,
                                 alternative = c("two.sided",
                                                "regular", "clustered"),
                                 method=c("Chisq", "MonteCarlo"),
                                 conditional=TRUE, ...,
                                 nsim=1999,
                                 Xcount=NULL,
                                 xbreaks=NULL, ybreaks=NULL, tess=NULL,
                                 fit=NULL, Xname=NULL, fitname=NULL) {
    trap.extra.arguments(...)
    method <- match.arg(method)
    alternative <- match.arg(alternative)
    if(method == "MonteCarlo") {
      check.1.real(nsim)
      explain.ifnot(nsim > 0)
    }
    if(is.null(Xcount))
      Xcount <- quadratcount(X, nx=nx, ny=ny, xbreaks=xbreaks, ybreaks=ybreaks,
                             tess=tess)
    tess <- attr(Xcount, "tess")
    testname <- switch(method,
                       Chisq = "Chi-squared test",
                       MonteCarlo = paste(
                         if(conditional) "Conditional" else "Unconditional",
                         "Monte Carlo test")
                       )
    # determine expected values under model
    if(is.null(fit)) {
      nullname <- "CSR"
      if(tess$type == "rect") 
        areas <- outer(diff(tess$xgrid), diff(tess$ygrid), "*")
      else 
        areas <- unlist(lapply(tiles(tess), area.owin))
      fitmeans <- sum(Xcount) * areas/sum(areas)
      df <- switch(method,
                   Chisq      = length(fitmeans) - 1,
                   MonteCarlo = NULL)
    } else {
      if(!is.ppm(fit))
        stop("fit should be a ppm object")
      if(!is.poisson.ppm(fit))
        stop("Quadrat test only supported for Poisson point process models")
      if(is.marked(fit))
        stop("Sorry, not yet implemented for marked point process models")
      nullname <- paste("fitted Poisson model", sQuote(fitname))
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
      df <- switch(method,
                   Chisq = {
                     length(fitmeans) - length(coef(fit))
                   },
                   MonteCarlo = {
                     NA
                   })
    }
    OBS <- as.vector(t(as.table(Xcount)))
    EXP <- as.vector(fitmeans)
    if(method == "Chisq") {
      # standard warnings
      if(df < 1)
        stop(paste("Not enough quadrats: degrees of freedom df =", df))
      if(any(EXP < 5))
        warning(paste("Some expected counts are small;",
                      "chi^2 approximation may be inaccurate"),
                call.=FALSE)
    }
    X2 <- sum((OBS - EXP)^2/EXP)
    names(X2) <- "X-squared"
    # p-value
    PVAL <- switch(method,
                   MonteCarlo = mcpval(X2,OBS,EXP,nsim,conditional,alternative),
                   Chisq      = {
                     pup <- pchisq(X2, df, lower.tail=FALSE)
                     plo <- pchisq(X2, df, lower.tail=TRUE)
                     switch(alternative,
                            regular   = plo,
                            clustered = pup,
                            two.sided = 2 * min(pup, plo))
                   })
    if(!is.null(df))
      names(df) <- "df"
    testname <- paste(testname, "of", nullname, "using quadrat counts")
    result <- structure(list(statistic = X2,
                             parameter = df,
                             p.value = PVAL,
                             method = testname,
                             data.name = Xname,
                             alternative = alternative,
                             observed = OBS,
                             expected = EXP,
                             residuals = (OBS - EXP)/sqrt(EXP)),
                        class = "htest")
    class(result) <- c("quadrattest", class(result))
    attr(result, "quadratcount") <- Xcount
    return(result)
  }

  mcpval <- function(stat,obsd,expctd,nsim,conditional,alternative) {
    nsim <- as.integer(nsim)
    if(conditional) {
      npts <- sum(obsd)
      p <- expctd/sum(expctd)
      X <- rmultinom(n=nsim,size=npts,prob=p)
    } else {
      ne <- length(expctd)
      X  <- matrix(rpois(nsim*ne,expctd),nrow=ne)
    }
    simstats <- apply((X-expctd)^2/expctd,2,sum)
    if(any(duplicated(simstats)))
      simstats <- jitter(simstats)
    phi <- (1 + sum(simstats >= stat))/(1+nsim)
    plo <- (1 + sum(simstats <= stat))/(1+nsim)
    pval <- switch(alternative,
                   clustered = phi,
                   regular   = plo,
                   two.sided = 2 * min(phi,plo))
    return(pval)
  }
  
  quadrat.testEngine
})

print.quadrattest <- function(x, ...) {
   NextMethod("print")
   cat("Quadrats: ")
   do.call("print.tess",
           resolve.defaults(list(as.tess(x)),
                            list(...),
                            list(brief=TRUE)))
   return(invisible(NULL))
}

plot.quadrattest <- local({
  
  plot.quadrattest <- function(x, ...) {
    xname <- short.deparse(substitute(x))
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
    # plot observed counts
    cos30 <- sqrt(2)/2
    sin30 <- 1/2
    f <- 0.4
    dotext(-f * cos30, f * sin30, as.vector(t(as.table(Xcount)))[ok],
           x0, y0, ra, adj=c(1,0), ...)
    # plot expected counts
    dotext(f * cos30, f * sin30, round(x$expected,1)[ok],
           x0, y0, ra, adj=c(0,0),
           ...)
    # plot Pearson residuals
    dotext(0, -f,  signif(x$residuals,2)[ok],x0, y0, ra, ...)
    return(invisible(NULL))
  }
  dotext <- function(dx, dy, values, x0, y0, ra, ...) {
    do.call.matched("text.default",
                    resolve.defaults(list(x=x0 + dx * ra, y = y0 + dy * ra),
                                     list(labels=paste(as.vector(values))),
                                     list(...)))
  }
  plot.quadrattest
})

