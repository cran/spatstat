#
# kppm.R
#
# kluster point process models
#
# $Revision: 1.4 $ $Date: 2008/07/18 04:28:53 $
#

kppm <- function(X, trend = ~1, clusters="Thomas", covariates=NULL, ...) {
  Xname <- deparse(substitute(X))
  clusters <- pickoption("cluster type", clusters,
                         c(Thomas="Thomas",
                           MatClust="MatClust"))
  if(is.marked(X))
    stop("Sorry, cannot handle marked point patterns")
  po <- ppm(X, trend=trend, covariates=covariates)
  stationary <- is.null(trend) || identical.formulae(trend, ~1)
  if(stationary) {
    lambda <- summary(po)$trend$value
    Ki <- Kest(X)
  } else {
    lambda <- predict(po, ngrid=spatstat.options("npixel"))
    Ki <- Kinhom(X, lambda)
  }
  switch(clusters,
         Thomas={
           startpar <- c(kappa=1, sigma2 = 4 * mean(nndist(X))^2)
           mcfit <- 
               thomas.estK(Ki, startpar=startpar, lambda=lambda, ...)
         },
         MatClust={
           startpar <- c(kappa=1, R = 2 * mean(nndist(X)))
           mcfit <- 
               matclust.estK(Ki, startpar=startpar, lambda=lambda, ...)
         })
  kappa <- mcfit$par[["kappa"]]
  mu <- if(stationary) lambda/kappa else eval.im(lambda/kappa)

  out <- list(stationary=stationary,
              clusters=clusters,
              Xname=Xname,
              X=X,
              po=po,
              lambda=lambda,
              mu=mu,
              Ki=Ki,
              mcfit=mcfit)
  class(out) <- c("kppm", class(out))
  return(out)
}

print.kppm <- function(x, ...) {
  ps <- summary(x$po)

  if(x$stationary)
    cat("Stationary cluster point process model\n")
  else 
    cat("Inhomogeneous cluster point process model\n")

  if(nchar(x$Xname) < 20)
    cat(paste("Fitted to point pattern dataset", sQuote(x$Xname), "\n"))

  if(!x$stationary) {
    cat("Trend formula:")
    print(ps$trend$formula)
    cat(paste("\n", ps$trend$label, ":", sep = ""))
    tv <- ps$trend$value
    if (is.list(tv)) {
      cat("\n")
      for (i in seq(tv)) print(tv[[i]])
    }
    else if (is.numeric(tv) && length(tv) == 1 && is.null(names(tv)))
      cat("\t", paste(tv), "\n")
    else {
      cat("\n")
      print(tv)
    }
  }
  mcfit <- x$mcfit
  cat(paste("Cluster model:", mcfit$info$modelname, "\n"))
  mp <- mcfit$modelpar
  if (!is.null(mp)) {
    cat("Fitted parameters:\n")
    print(mp[!is.na(mp)])
  }
  invisible(NULL)
}

plot.kppm <- function(x, ...) {
  modelname <- deparse(substitute(x))
  if(x$stationary)
    plot(x$mcfit, main=c(modelname, "K-function"))
  else {
    if(interactive()) opa <- par(ask=TRUE)
    plot(x$po, main=c(modelname, "Intensity"), how="image")
    plot(x$mcfit, main=c(modelname, "Inhomogeneous K-function"))
    if(interactive()) par(opa)
  }
  return(invisible(NULL))
}

predict.kppm <- function(object, ...) {
  predict(object$po, ...)
}

simulate.kppm <- function(object, nsim=1, ...) {
  out <- list()
  win <- object$X$window
  switch(object$clusters,
         Thomas={
           mp <- as.list(object$mcfit$modelpar)
           kappa <- mp$kappa
           sigma <- mp$sigma
           mu    <- object$mu
           for(i in 1:nsim)
             out[[i]] <- rThomas(kappa,sigma,mu,win)
         },
         MatClust={
           mp <- as.list(object$mcfit$modelpar)
           kappa <- mp$kappa
           r <-     mp$R
           mu    <- object$mu
           for(i in 1:nsim)
             out[[i]] <- rMatClust(kappa,r,mu,win)
         })
  class(out) <- c("listof", class(out))
  return(out)
}

formula.kppm <- function(x, ...) {
  formula(x$po, ...)
}

terms.kppm <- function(x, ...) {
  terms(x$po, ...)
}

update.kppm <- function(object, trend=~1, ..., clusters=NULL) {
  trend <- update(formula(object), trend)
  if(is.null(clusters))
    clusters <- object$clusters
  out <- kppm(object$X, trend, clusters, ...)
  out$Xname <- object$Xname
  return(out)
}
