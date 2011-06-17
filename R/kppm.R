#
# kppm.R
#
# kluster/kox point process models
#
# $Revision: 1.38 $ $Date: 2011/06/14 07:39:50 $
#

kppm <- function(X, trend = ~1, clusters="Thomas", covariates=NULL, ...,
                 statistic="K", statargs=list()) {
  Xname <- deparse(substitute(X))
  clusters <- pickoption("cluster type", clusters,
                         c(Thomas="Thomas",
                           MatClust="MatClust",
                           LGCP="LGCP"))
  statistic <- pickoption("summary statistic", statistic,
                          c(K="K", g="pcf", pcf="pcf"))
  if(is.marked(X))
    stop("Sorry, cannot handle marked point patterns")
  po <- ppm(X, trend=trend, covariates=covariates)
  stationary <- is.stationary(po)
  if(stationary) {
    lambda <- summary(po)$trend$value
    StatFun <- if(statistic == "K") "Kest" else "pcf"
    StatName <- if(statistic == "K") "K-function" else "pair correlation function"
    Stat <- do.call(StatFun, append(list(X), statargs))
  } else {
    lambda <- predict(po, ngrid=rev(spatstat.options("npixel")))
    StatFun <- if(statistic == "K") "Kinhom" else "pcfinhom"
    StatName <- if(statistic == "K") "inhomogeneous K-function" else "inhomogeneous pair correlation function"
    Stat <- do.call(StatFun, append(list(X, lambda), statargs))
  }
  switch(clusters,
         Thomas={
           startpar0 <- c(kappa=1, sigma2 = 4 * mean(nndist(X))^2)
           FitFun <- if(statistic == "K") "thomas.estK" else "thomas.estpcf"
           mcfit <-
             do.call(FitFun,
                     resolve.defaults(
                                      list(X=Stat, lambda=lambda),
                                      list(...),
                                      list(startpar=startpar0, dataname=Xname)))
           # kappa = mean number of points per cluster
           kappa <- mcfit$par[["kappa"]]
           # mu = intensity of parents
           mu <- if(stationary) lambda/kappa else eval.im(lambda/kappa)
           isPCP <- TRUE
         },
         MatClust={
           startpar0 <- c(kappa=1, R = 2 * mean(nndist(X)))
           FitFun <- if(statistic == "K") "matclust.estK" else "matclust.estpcf"
           mcfit <-
             do.call(FitFun,
                     resolve.defaults(
                                      list(X=Stat, lambda=lambda),
                                      list(...),
                                      list(startpar=startpar0, dataname=Xname)))
           # kappa = mean number of points per cluster
           kappa <- mcfit$par[["kappa"]]
           # mu = intensity of parents
           mu <- if(stationary) lambda/kappa else eval.im(lambda/kappa)
           isPCP <- TRUE
         },
         LGCP={
           startpar0 <- c(sigma2=1, alpha=2 * mean(nndist(X)))
           FitFun <- if(statistic == "K") "lgcp.estK" else "lgcp.estpcf"
           mcfit <-
             do.call(FitFun,
                     resolve.defaults(
                                      list(X=Stat, lambda=lambda),
                                      list(...),
                                      list(startpar=startpar0, dataname=Xname)))
           sigma2 <- mcfit$par[["sigma2"]]
           # mu = mean of log intensity 
           mu <- if(stationary) log(lambda) - sigma2/2 else eval.im(log(lambda) - sigma2/2)
           isPCP <- FALSE
         })


  out <- list(stationary=stationary,
              clusters=clusters,
              isPCP=isPCP,
              Xname=Xname,
              statistic=statistic,
              X=X,
              po=po,
              lambda=lambda,
              mu=mu,
              Stat=Stat,
              StatName=StatName,
              mcfit=mcfit)
  class(out) <- c("kppm", class(out))
  return(out)
}

print.kppm <- function(x, ...) {

  isPCP <- x$isPCP
  # handle outdated objects - which were all cluster processes
  if(is.null(isPCP)) isPCP <- TRUE
  
  cat(paste(if(x$stationary) "Stationary" else "Inhomogeneous",
            if(isPCP) "cluster" else "Cox",
            "point process model\n"))

  if(nchar(x$Xname) < 20)
    cat(paste("Fitted to point pattern dataset", sQuote(x$Xname), "\n"))

  cat(paste("Fitted using the", x$StatName, "\n"))

  # ............... trend .........................

  print(x$po, what="trend")

  # ..................... clusters ................
  
  mcfit <- x$mcfit
  cat(paste(if(isPCP) "Cluster" else "Cox",
            "model:", mcfit$info$modelname, "\n"))
  cm <- mcfit$covmodel
  if(!is.null(cm)) {
    # Covariance model - LGCP only
    cat(paste("\tCovariance model:", cm$model, "\n"))
    margs <- cm$margs
    if(!is.null(margs)) {
      nama <- names(margs)
      tags <- ifelse(nzchar(nama), paste(nama, "="), "")
      tagvalue <- paste(tags, margs)
      cat(paste("\tCovariance parameters:",
                paste(tagvalue, collapse=", "),
                "\n"))
    }
  }
  pa <- mcfit$par
  if (!is.null(pa)) {
    cat(paste("Fitted",
              if(isPCP) "cluster" else "correlation",
              "parameters:\n"))
    print(pa)
  }

  if(!is.null(mu <- x$mu)) {
    if(isPCP) {
      cat("Mean cluster size: ")
      if(!is.im(mu)) cat(mu, "points\n") else cat("[pixel image]\n")
    } else {
      cat("Fitted mean of log of random intensity: ")
      if(!is.im(mu)) cat(mu, "\n") else cat("[pixel image]\n")
    }
  }
  invisible(NULL)
}

plot.kppm <- function(x, ..., what=c("intensity", "statistic")) {
  modelname <- deparse(substitute(x))
  plotem <- function(x, ..., main=dmain, dmain) { plot(x, ..., main=main) }
  what <- pickoption("plot type", what,
                    c(statistic="statistic",
                      intensity="intensity"),
                    multi=TRUE)
  if(x$stationary && ("statistic" %in% what))
    plotem(x$mcfit, ..., dmain=c(modelname, x$StatName))
  else {
    pauseit <- interactive() && (length(what) > 1)
    if(pauseit) opa <- par(ask=TRUE)
    for(style in what)
      switch(style,
             intensity={
               plotem(x$po, ...,
                      dmain=c(modelname, "Intensity"),
                      how="image")
             },
             statistic={
               plotem(x$mcfit, ...,
                      dmain=c(modelname, x$StatName))
             })
    if(pauseit) par(opa)
  }
  return(invisible(NULL))
}

predict.kppm <- function(object, ...) {
  predict(object$po, ...)
}

simulate.kppm <- function(object, nsim=1, seed=NULL, ...,
                          window=NULL, covariates=NULL) {
  # .... copied from simulate.lm ....
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if (is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  # ..................................
  # determine parameters
  if(is.null(window)) {
    win <- object$X$window
  } else {
    stopifnot(is.owin(window))
    win <- window
  }
  mp <- as.list(object$mcfit$modelpar)

  # parameter 'mu'
  # = parent intensity of cluster process
  # = mean log intensity of log-Gaussian Cox process
  
  if(is.null(covariates) && (object$stationary || is.null(window))) {
    # use existing 'mu' (scalar or image)
    mu <- object$mu
  } else {
    # recompute 'mu' using new data
    switch(object$clusters,
           Thomas=,
           MatClust={
             # Poisson cluster process
             kappa <- mp$kappa
             lambda <- predict(object, window=win, covariates=covariates)
             mu <- eval.im(lambda/kappa)
           },
           LGCP={
             # log-Gaussian Cox process
             sigma2 <- mp$sigma2
             lambda <- predict(object, window=win, covariates=covariates)
             mu <- eval.im(log(lambda) - sigma2/2)
           }
           )
  }
  
  out <- list()
  switch(object$clusters,
         Thomas={
           kappa <- mp$kappa
           sigma <- mp$sigma
           for(i in 1:nsim)
             out[[i]] <- rThomas(kappa,sigma,mu,win)
         },
         MatClust={
           kappa <- mp$kappa
           r <-     mp$R
           for(i in 1:nsim)
             out[[i]] <- rMatClust(kappa,r,mu,win)
         },
         LGCP={
           sigma2 <- mp$sigma2
           alpha  <- mp$alpha
           cm <- object$mcfit$covmodel
           model <- cm$model
           margs <- cm$margs
           param <- c(0, sigma2, 0, alpha, unlist(margs))
           for(i in 1:nsim)
             out[[i]] <- rLGCP(model=model, mu=mu, param=param, ..., win=win)
         })
  out <- as.listof(out)
  names(out) <- paste("Simulation", 1:nsim)
  attr(out, "seed") <- RNGstate
  return(out)
}

formula.kppm <- function(x, ...) {
  formula(x$po, ...)
}

terms.kppm <- function(x, ...) {
  terms(x$po, ...)
}

update.kppm <- function(object, trend=~1, ..., clusters=NULL) {
  if(!missing(trend))
    trend <- update(formula(object), trend)
  if(is.null(clusters))
    clusters <- object$clusters
  out <- do.call(kppm,
                 resolve.defaults(list(trend=trend, clusters=clusters),
                                  list(...),
                                  list(X=object$X)))
  out$Xname <- object$Xname
  return(out)
}

unitname.kppm <- function(x) {
  return(unitname(x$X))
}

"unitname<-.kppm" <- function(x, value) {
  unitname(x$X) <- value
  unitname(x$mcfit) <- value
  return(x)
}

coef.kppm <- function(object, ...) {
  return(coef(object$po))
}


Kmodel.kppm <- function(model, ...) {
  Kpcf.kppm(model, what="K")
}

pcfmodel.kppm <- function(model, ...) {
  Kpcf.kppm(model, what="pcf")
}

Kpcf.kppm <- function(model, what=c("K", "pcf")) {
  what <- match.arg(what)
  # Extract function definition from internal table
  clusters <- model$clusters
  tableentry <- .Spatstat.ClusterModelInfo[[clusters]]
  if(is.null(tableentry))
    stop("No information available for", sQuote(clusters), "cluster model")
  fun <- tableentry[[what]]
  if(is.null(fun))
    stop("No expression available for", what, "for", sQuote(clusters),
         "cluster model")
  # Extract model parameters
  par <- model$mcfit$par
  # Extract auxiliary definitions (if applicable)
  funaux <- tableentry$funaux
  # Extract covariance model (if applicable)
  cm <- model$mcfit$covmodel
  model <- cm$model
  margs <- cm$margs
  #
  f <- function(r) as.numeric(fun(par=par, rvals=r,
                                  funaux=funaux, model=model, margs=margs))
  return(f)
}

is.stationary.kppm <- function(x) {
  return(x$stationary)
}

is.poisson.kppm <- function(x) {
  switch(x$clusters,
         Thomas=,
         MatClust={
           # Poisson cluster process
           mu <- x$mu
           return(!is.null(mu) && (max(mu) == 0))
         },
         LGCP = {
           # log-Gaussian Cox process
           sigma2 <- x$mcfit$par[["sigma2"]]
           return(sigma2 == 0)
         },
         return(FALSE))
}
