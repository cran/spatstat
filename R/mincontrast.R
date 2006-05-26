#
#  mincontrast.R
#
#  Functions for estimation by minimum contrast
#

##################  base ################################

mincontrast <- function(observed, theoretical, startpar,
                        ...,
                        ctrl=list(q = 1/4, p = 2, rmin=NULL, rmax=NULL),
                        fvlab=list(label=NULL, desc="minimum contrast fit"),
                        explain=list(dataname=NULL, modelname=NULL, fname=NULL)) {
  verifyclass(observed, "fv")
  stopifnot(is.function(theoretical))
  if(!any("par" %in% names(formals(theoretical))))
    stop(paste("Theoretical function does not include an argument called",
               sQuote("par")))

  # enforce defaults
  ctrl <- resolve.defaults(ctrl, list(q = 1/4, p = 2, rmin=NULL, rmax=NULL))
  fvlab <- resolve.defaults(fvlab,
                            list(label=NULL, desc="minimum contrast fit"))
  explain <- resolve.defaults(explain,
                              list(dataname=NULL, modelname=NULL, fname=NULL))
  
  # determine range of r values
  rmin <- ctrl$rmin
  rmax <- ctrl$rmax
  if(!is.null(rmin) && !is.null(rmax))
    stopifnot(rmin < rmax && rmin >= 0)
  else {
    alim <- attr(observed, "alim")
    if(is.null(rmin)) rmin <- alim[1]
    if(is.null(rmax)) rmax <- alim[2]
  }
  # extract vector of r values
  argu <- attr(observed, "argu")
  rvals <- observed[[argu]]
  # extract vector of observed values of statistic
  valu <- attr(observed, "valu")
  obs <- observed[[valu]]
  # restrict to [rmin, rmax]
  if(max(rvals) < rmax)
    stop(paste("rmax=", signif(rmax,4),
               "exceeds the range of available data",
               "= [", signif(min(rvals),4), ",", signif(max(rvals),4), "]"))
  sub <- (rvals >= rmin) & (rvals <= rmax)
  rvals <- rvals[sub]
  obs <- obs[sub]
  # for efficiency
  obsq <- obs^(ctrl$q)
  # define objective function
  objective <- function(par, obsq, theoretical, rvals, qq, pp, rmin, rmax, ...) {
    theo <- theoretical(par=par, rvals, ...)
    if(!is.vector(theo) || !is.numeric(theo))
      stop("theoretical function did not return a numeric vector")
    if(length(theo) != length(obs))
      stop("theoretical function did not return the correct number of values")
    discrep <- (abs(theo^qq - obsq))^pp
    return(sum(discrep))
  }
  # go
  minimum <- optim(startpar, fn=objective,
                   obsq=obsq, theoretical=theoretical,
                   rvals=rvals,
                   qq=ctrl$q, pp=ctrl$p, rmin=rmin, rmax=rmax, ...)
  # evaluate the fitted theoretical curve
  fittheo <- theoretical(minimum$par, rvals, ...)
  # pack it up as an `fv' object
  label <- fvlab$label
  desc  <- fvlab$desc
  if(is.null(label))
    label <- paste("fit(", argu, ")", collapse="")
  fitfv <- bind.fv(observed[sub, ],
                   data.frame(fit=fittheo),
                   label, desc)
  result <- list(par=minimum$par,
                 fit=fitfv,
                 opt=minimum,
                 ctrl=list(p=ctrl$p,q=ctrl$q,rmin=rmin,rmax=rmax),
                 info=explain
                 )
  
  class(result) <- c("minconfit", class(result))
  return(result)
}

print.minconfit <- function(x, ...) {
  # explanatory
  cat(paste("Minimum contrast fit ",
            "(",
            "object of class ",
            dQuote("minconfit"),
            ")",
            "\n", sep=""))
  mo <- x$info$modelname
  fu <- x$info$fname
  da <- x$info$dataname
  if(!is.null(mo))
    cat(paste("Model:", mo, "\n"))
  if(!is.null(fu) && !is.null(da))
    cat(paste("Fitted by matching theoretical", fu, "function to", da))
  else {
    if(!is.null(fu))
      cat(paste(" based on", fu))
    if(!is.null(da))
      cat(paste(" fitted to", da))
  }
  cat("\n")
  # Values
  cat("Parameters fitted by minimum contrast ($par):\n")
  print(x$par, ...)
  mp <- x$modelpar
  if(!is.null(mp)) {
    cat(paste("Derived parameters of",
              if(!is.null(mo)) mo else "model",
              "($modelpar):\n"))
    print(mp)
  }
  # Diagnostics
  cgce <- x$opt$convergence
  switch(paste(cgce),
         "0"={
           cat(paste("Converged successfully after",
                     x$opt$counts[["function"]],
                     "iterations.\n"))
         },
         "1"={
           cat("Warning: iteration limit maxit was reached.\n")
         },
         "10"={
           cat("Warning: Nelder-Mead simplex was degenerate\n")
         },
         "51"={
           cat(paste("Warning from L-BGFS-B method:\n",
                     x$opt$message, "\n"))
         },
         "52"={
           cat(paste("Error message from L-BGFS-B method:\n",
                     x$opt$message, "\n"))
         },
         cat(paste("Unrecognised error code", cgce, "\n"))
         )
  # Algorithm parameters
  ct <- x$ctrl
  cat(paste("Domain of integration:",
            "[",
            signif(ct$rmin,4),
            ",",
            signif(ct$rmax,4),
            "]\n"))
  cat(paste("Exponents:",
            "p=", paste(signif(ct$p, 3), ",",  sep=""),
            "q=", signif(ct$q,3), "\n"))
  invisible(NULL)
}
              

plot.minconfit <- function(x, ...) {
  plot(x$fit, ...)
}


############### applications (specific models) ##################


thomas.estK <- function(X, startpar=c(kappa=1,sigma2=1),
                        lambda=NULL, q=1/4, p=2, rmin=NULL, rmax=NULL){

  dataname <- deparse(substitute(X))

  if(inherits(X, "fv")) {
    K <- X
  } else if(inherits(X, "ppp")) {
    K <- Kest(X)
    dataname <- paste("Kest(", dataname, ")", sep="")
    lambda <- summary(X)$intensity
  } else 
    stop("Unrecognised format for argument X")

  startpar <- check.named.vector(startpar, c("kappa","sigma2"))

  theoret <- function(par,rvals){
    if(any(par <= 0))
      return(rep(Inf, length(rvals)))
    pi*rvals^2+(1-exp(-rvals^2/(4*par[2])))/par[1]
  }
  result <- mincontrast(K, theoret, startpar,
                        ctrl=list(q=q, p=p,rmin=rmin, rmax=rmax),
                        fvlab=list(label="Kfit(r)",
                          desc="minimum contrast fit of Thomas process"),
                        explain=list(dataname=dataname,
                          fname="K", modelname="Thomas process"))
  # imbue with meaning
  par <- result$par
  names(par) <- c("kappa", "sigma2")
  result$par <- par
  # infer model parameters
  mu <- if(is.numeric(lambda) && length(lambda) == 1)
           lambda/result$par[["kappa"]] else NA
  result$modelpar <- c(kappa=par[["kappa"]],
                       sigma=sqrt(par[["sigma2"]]),
                       mu=mu)
  result$internal <- list(model="Thomas")
  return(result)
}


lgcp.estK <- function(X, startpar=list(sigma2=1,alpha=1),
                      lambda=NULL, q=1/4, p=2, rmin=NULL, rmax=NULL) {
  
  dataname <- deparse(substitute(K))
  if(inherits(X, "fv")) {
    K <- X
  } else if(inherits(X, "ppp")) {
    K <- Kest(X)
    dataname <- paste("Kest(", dataname, ")", sep="")
    lambda <- summary(X)$intensity
  } else 
    stop("Unrecognised format for argument X")

  startpar <- check.named.vector(startpar, c("sigma2","alpha"))

  integrand<-function(r,par){
    2*pi*r*exp(par[1]*exp(-r/par[2]))
  }
  theoret <- function(par, rvals) {
    if(any(par <= 0))
      return(rep(Inf, length(rvals)))
    th <- numeric(length(rvals))
    th[1] <- if(rvals[1] == 0) 0 else 
             integrate(integrand,lower=0,upper=rvals[1],par=par)$value
    for (i in 2:length(rvals))
      th[i]=th[i-1]+integrate(integrand,lower=rvals[i-1],upper=rvals[i],par=par)$value
    return(th)
  }
  result <- mincontrast(K, theoret, startpar,
                        ctrl=list(q=q, p=p, rmin=rmin, rmax=rmax),
                        fvlab=list(label="Kfit(r)",
                          desc="minimum contrast fit of LGCP"),
                        explain=list(dataname=dataname,
                          fname="K",
                          modelname="log-Gaussian Cox process"))
  # imbue with meaning
  par <- result$par
  names(par) <- c("sigma2", "alpha")
  result$par <- par
  # infer model parameters
  mu <- if(is.numeric(lambda) && length(lambda) == 1 && lambda > 0)
           log(lambda) - par[["sigma2"]]/2 else NA
  result$modelpar <- c(sigma2=par[["sigma2"]],
                       alpha=par[["alpha"]],
                       mu=mu)
  result$internal <- list(model="lcgp")
  return(result)
}


matclust.estK <- function(X, startpar=c(kappa=1,R=1),
                        lambda=NULL, q=1/4, p=2, rmin=NULL, rmax=NULL){

  dataname <- deparse(substitute(X))

  if(inherits(X, "fv")) {
    K <- X
  } else if(inherits(X, "ppp")) {
    K <- Kest(X)
    dataname <- paste("Kest(", dataname, ")", sep="")
    lambda <- summary(X)$intensity
  } else 
    stop("Unrecognised format for argument X")

  startpar <- check.named.vector(startpar, c("kappa","R"))

  hfun <- function(zz) {
    ok <- (zz < 1)
    h <- numeric(length(zz))
    h[!ok] <- 1
    z <- zz[ok]
    h[ok] <- 2 + (1/pi) * (
                           (8 * z^2 - 4) * acos(z)
                           - 2 * asin(z)
                           + 4 * z * sqrt((1 - z^2)^3)
                           - 6 * z * sqrt(1 - z^2)
                           )
    return(h)
  }
  theoret <- function(par,rvals){
    if(any(par <= 0))
      return(rep(Inf, length(rvals)))
    pi * rvals^2 + (1/par[1]) * hfun(rvals/(2 * par[2]))
  }
  result <- mincontrast(K, theoret, startpar,
                        ctrl=list(q=q, p=p,rmin=rmin, rmax=rmax),
                        fvlab=list(label="Kfit(r)",
                          desc="minimum contrast fit of Matern Cluster process"),
                        explain=list(dataname=dataname,
                          fname="K", modelname="Matern Cluster process"))
  # imbue with meaning
  par <- result$par
  names(par) <- c("kappa", "R")
  result$par <- par
  # infer model parameters
  mu <- if(is.numeric(lambda) && length(lambda) == 1)
           lambda/result$par[["kappa"]] else NA
  result$modelpar <- c(kappa=par[["kappa"]],
                       R=par[["R"]],
                       mu=mu)
  result$internal <- list(model="MatClust")
  return(result)
}


check.named.vector <- function(x, nam) {
  xname <- deparse(substitute(x))
  if(!is.numeric(x) || length(x) != 2 || !all(nam %in% names(x))) {
    whinge <- paste(xname,
                    "must be a named vector, with components",
                    commasep(nam))
    stop(whinge)
  }
  return(x[nam])
}

