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
  argu <- fvnames(observed, ".x")
  rvals <- observed[[argu]]
  # extract vector of observed values of statistic
  valu <- fvnames(observed, ".y")
  obs <- observed[[valu]]
  # restrict to [rmin, rmax]
  if(max(rvals) < rmax)
    stop(paste("rmax=", signif(rmax,4),
               "exceeds the range of available data",
               "= [", signif(min(rvals),4), ",", signif(max(rvals),4), "]"))
  sub <- (rvals >= rmin) & (rvals <= rmax)
  rvals <- rvals[sub]
  obs <- obs[sub]
  # sanity clause
  if(!all(ok <- is.finite(obs))) {
    whinge <- paste("Some values of the empirical function",
                    sQuote(explain$fname),
                    "were infinite or NA.")
    iMAX <- max(which(ok))
    iMIN <- min(which(!ok)) + 1
    if(iMAX > iMIN && all(ok[iMIN:iMAX])) {
      rmin <- rvals[iMIN]
      rmax <- rvals[iMAX]
      obs   <- obs[iMIN:iMAX]
      rvals <- rvals[iMIN:iMAX]
      sub[sub] <- ok
      warning(paste(whinge,
                    "Range of r values was reset to",
                    prange(c(rmin, rmax))),
              call.=FALSE)
    } else stop(paste(whinge, "Please choose a narrower range [rmin, rmax]"),
                call.=FALSE)
  }
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
  cm <- x$covmodel
  if(!is.null(mo))
    cat(paste("Model:", mo, "\n"))
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
  xname <- deparse(substitute(x))
  do.call("plot.fv",
          resolve.defaults(list(x$fit),
                           list(...),
                           list(main=xname)))
}

unitname.minconfit <- function(x) {
  unitname(x$fit)
}

"unitname<-.minconfit" <- function(x, value) {
  unitname(x$fit) <- value
  return(x)
}


############### applications (specific models) ##################


# table of explicitly-known K functions and pcf

.Spatstat.ClusterModelInfo <-
  list(
       Thomas=list(
         # Thomas process: par = (kappa, sigma2)
         K = function(par,rvals, ...){
           if(any(par <= 0))
             return(rep(Inf, length(rvals)))
           pi*rvals^2+(1-exp(-rvals^2/(4*par[2])))/par[1]
         },
         pcf= function(par,rvals, ...){
           if(any(par <= 0))
             return(rep(Inf, length(rvals)))
           1 + exp(-rvals^2/(4 * par[2]))/(4 * pi * par[1] * par[2])
         }
         ),
       MatClust=list(
         # Matern Cluster process: par = (kappa, R)
         K = function(par,rvals, ..., funaux){
           if(any(par <= 0))
             return(rep(Inf, length(rvals)))
           hfun <- funaux$Hfun
           pi * rvals^2 + (1/par[1]) * hfun(rvals/(2 * par[2]))
         },
         pcf= function(par,rvals, ..., funaux){
             if(any(par <= 0))
               return(rep(Inf, length(rvals)))
             doh <- funaux$DOH
             y <- (1/(4*pi*rvals * par[1] * par[2])) * doh(rvals/(2 * par[2]))
             ifelse(rvals < .Machine$double.eps, 1, 1+y)
           },
         funaux=list(
           Hfun=function(zz) {
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
           },
           DOH=function(zz) {
             ok <- (zz < 1)
             h <- numeric(length(zz))
             h[!ok] <- 0
             z <- zz[ok]
             h[ok] <- (16/pi) * (z * acos(z) - (z^2) * sqrt(1 - z^2))
             return(h)
           })

         ),
       LGCP=list(
         # Log Gaussian Cox process: par = (sigma2, alpha)
         # calls Covariance() from RandomFields package
         K = function(par, rvals, ..., model, margs) {
           if(any(par <= 0))
             return(rep(Inf, length(rvals)))
           if(model == "exponential") {
             # For efficiency
             integrand <- function(r,par,...) 2*pi*r*exp(par[1]*exp(-r/par[2]))
           } else {
             integrand <- function(r,par,model,margs)
               2*pi *r *exp(Covariance(r,model=model,
                                       param=c(0.0,par[1],0.0,par[2],margs)))
           }
           th <- numeric(length(rvals))
           th[1] <- if(rvals[1] == 0) 0 else 
                    integrate(integrand,lower=0,upper=rvals[1],
                              par=par,model=model,margs=margs)$value
           for (i in 2:length(rvals)) {
             delta <- integrate(integrand,
                                lower=rvals[i-1],upper=rvals[i],
                                par=par,model=model,margs=margs)
             th[i]=th[i-1]+delta$value
           }
           return(th)
         },
         pcf= function(par, rvals, ..., model, margs) {
           if(any(par <= 0))
             return(rep(Inf, length(rvals)))
           if(model == "exponential") {
             # For efficiency and to avoid need for RandomFields package
             gtheo <- exp(par[1]*exp(-rvals/par[2]))
           } else {
             gtheo <- exp(Covariance(rvals,model=model,
                                     param=c(0.0,par[1],0.0,par[2],margs)))
           }
           return(gtheo)
         },
         parhandler=function(model, ...) {
           if(!is.character(model))
             stop("Covariance function model should be specified by name")
           margs <- c(...)
           if(model != "exponential") {
             # check validity
             ok <- try(Covariance(0, model=model,param=c(0,1,0,1,margs)))
             if(inherits(ok, "try-error"))
               stop("Error in evaluating Covariance")
           }
           return(list(model=model, margs=margs))
         }
         )
  )

getdataname <- function(defaultvalue, ..., dataname=NULL) {
  if(!is.null(dataname)) dataname else defaultvalue
}
  
thomas.estK <- function(X, startpar=c(kappa=1,sigma2=1),
                        lambda=NULL, q=1/4, p=2, rmin=NULL, rmax=NULL, ...) {

  dataname <- getdataname(deparse(substitute(X), width=20, nlines=1), ...)

  if(inherits(X, "fv")) {
    K <- X
    if(!(attr(K, "fname") %in% c("K", "K[inhom]")))
      warning("Argument X does not appear to be a K-function")
  } else if(inherits(X, "ppp")) {
    K <- Kest(X)
    dataname <- paste("Kest(", dataname, ")", sep="")
    if(is.null(lambda))
      lambda <- summary(X)$intensity
  } else 
    stop("Unrecognised format for argument X")

  startpar <- check.named.vector(startpar, c("kappa","sigma2"))

  theoret <- .Spatstat.ClusterModelInfo$Thomas$K
  
  result <- mincontrast(K, theoret, startpar,
                        ctrl=list(q=q, p=p,rmin=rmin, rmax=rmax),
                        fvlab=list(label="%s[fit](r)",
                          desc="minimum contrast fit of Thomas process"),
                        explain=list(dataname=dataname,
                          fname=attr(K, "fname"),
                          modelname="Thomas process"), ...)
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

lgcp.estK <- function(X, startpar=c(sigma2=1,alpha=1),
                      covmodel=list(model="exponential"), 
                      lambda=NULL, q=1/4, p=2, rmin=NULL, rmax=NULL, ...) {

  dataname <- getdataname(deparse(substitute(X), width=20, nlines=1), ...)
  
  if(inherits(X, "fv")) {
    K <- X
    if(!(attr(K, "fname") %in% c("K", "K[inhom]")))
      warning("Argument X does not appear to be a K-function")
  } else if(inherits(X, "ppp")) {
    K <- Kest(X)
    dataname <- paste("Kest(", dataname, ")", sep="")
    if(is.null(lambda))
      lambda <- summary(X)$intensity
  } else 
    stop("Unrecognised format for argument X")

  startpar <- check.named.vector(startpar, c("sigma2","alpha"))

  # digest parameters of Covariance model and test validity
  ph <- .Spatstat.ClusterModelInfo$LGCP$parhandler
  cmodel <- do.call(ph, covmodel)
  
  theoret <- .Spatstat.ClusterModelInfo$LGCP$K

  result <- mincontrast(K, theoret, startpar,
                        ctrl=list(q=q, p=p, rmin=rmin, rmax=rmax),
                        fvlab=list(label="%s[fit](r)",
                          desc="minimum contrast fit of LGCP"),
                        explain=list(dataname=dataname,
                          fname=attr(K, "fname"),
                          modelname="log-Gaussian Cox process"),
                        ...,
                        model=cmodel$model,
                        margs=cmodel$margs)
  # imbue with meaning
  par <- result$par
  names(par) <- c("sigma2", "alpha")
  result$par <- par
  result$covmodel <- cmodel
  # infer model parameters
  mu <- if(is.numeric(lambda) && length(lambda) == 1 && lambda > 0)
           log(lambda) - par[["sigma2"]]/2 else NA
  result$modelpar <- c(sigma2=par[["sigma2"]],
                       alpha=par[["alpha"]],
                       mu=mu)
  result$internal <- list(model="lgcp")
  return(result)
}



matclust.estK <- function(X, startpar=c(kappa=1,R=1),
                          lambda=NULL, q=1/4, p=2, rmin=NULL, rmax=NULL, ...) {

  dataname <- getdataname(deparse(substitute(X), width=20, nlines=1), ...)

  if(inherits(X, "fv")) {
    K <- X
    if(!(attr(K, "fname") %in% c("K", "K[inhom]")))
      warning("Argument X does not appear to be a K-function")
  } else if(inherits(X, "ppp")) {
    K <- Kest(X)
    dataname <- paste("Kest(", dataname, ")", sep="")
    if(is.null(lambda))
      lambda <- summary(X)$intensity
  } else 
    stop("Unrecognised format for argument X")

  startpar <- check.named.vector(startpar, c("kappa","R"))

  theoret <- .Spatstat.ClusterModelInfo$MatClust$K
  funaux <- .Spatstat.ClusterModelInfo$MatClust$funaux
  
  result <- mincontrast(K, theoret, startpar,
                        ctrl=list(q=q, p=p,rmin=rmin, rmax=rmax),
                        fvlab=list(label="%s[fit](r)",
                          desc="minimum contrast fit of Matern Cluster process"),
                        explain=list(dataname=dataname,
                          fname=attr(K, "fname"),
                          modelname="Matern Cluster process"),
                        ...,
                        funaux=funaux)
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


## versions using pcf (suggested by Jan Wild)

thomas.estpcf <- function(X, startpar=c(kappa=1,sigma2=1),
                          lambda=NULL, q=1/4, p=2, rmin=NULL, rmax=NULL, ...,
                          pcfargs=list()){

  dataname <- getdataname(deparse(substitute(X), width=20, nlines=1), ...)

  if(inherits(X, "fv")) {
    g <- X
    if(!(attr(g, "fname") %in% c("g", "g[inhom]")))
      warning("Argument X does not appear to be a pair correlation function")
  } else if(inherits(X, "ppp")) {
    g <- do.call("pcf.ppp", append(list(X), pcfargs))
    dataname <- paste("pcf(", dataname, ")", sep="")
    if(is.null(lambda))
      lambda <- summary(X)$intensity
  } else 
    stop("Unrecognised format for argument X")

  startpar <- check.named.vector(startpar, c("kappa","sigma2"))

  theoret <- .Spatstat.ClusterModelInfo$Thomas$pcf
  
  # avoid using g(0) as it may be infinite
  argu <- fvnames(g, ".x")
  rvals <- g[[argu]]
  if(rvals[1] == 0 && (is.null(rmin) || rmin == 0)) {
    rmin <- rvals[2]
  }
  result <- mincontrast(g, theoret, startpar,
                        ctrl=list(q=q, p=p,rmin=rmin, rmax=rmax),
                        fvlab=list(
                          label="%s[fit](r)",
                          desc="minimum contrast fit of Thomas process"),
                        explain=list(
                          dataname=dataname,
                          fname=attr(g, "fname"),
                          modelname="Thomas process"), ...)
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

matclust.estpcf <- function(X, startpar=c(kappa=1,R=1),
                            lambda=NULL, q=1/4, p=2, rmin=NULL, rmax=NULL, ...,
                            pcfargs=list()){

  dataname <- getdataname(deparse(substitute(X), width=20, nlines=1), ...)

  if(inherits(X, "fv")) {
    g <- X
    if(!(attr(g, "fname") %in% c("g", "g[inhom]")))
      warning("Argument X does not appear to be a pair correlation function")
  } else if(inherits(X, "ppp")) {
    g <- do.call("pcf.ppp", append(list(X), pcfargs))
    dataname <- paste("pcf(", dataname, ")", sep="")
    if(is.null(lambda))
      lambda <- summary(X)$intensity
  } else 
    stop("Unrecognised format for argument X")

  startpar <- check.named.vector(startpar, c("kappa","R"))

  theoret <- .Spatstat.ClusterModelInfo$MatClust$pcf
  funaux <- .Spatstat.ClusterModelInfo$MatClust$funaux
  
  # avoid using g(0) as it may be infinite
  argu <- fvnames(g, ".x")
  rvals <- g[[argu]]
  if(rvals[1] == 0 && (is.null(rmin) || rmin == 0)) {
    rmin <- rvals[2]
  }
  result <- mincontrast(g, theoret, startpar,
                        ctrl=list(q=q, p=p,rmin=rmin, rmax=rmax),
                        fvlab=list(label="%s[fit](r)",
                          desc="minimum contrast fit of Matern Cluster process"),
                        explain=list(dataname=dataname,
                          fname=attr(g, "fname"),
                          modelname="Matern Cluster process"),
                        ...,
                        funaux=funaux)
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

lgcp.estpcf <- function(X, startpar=c(sigma2=1,alpha=1),
                      covmodel=list(model="exponential"), 
                        lambda=NULL, q=1/4, p=2, rmin=NULL, rmax=NULL, ...,
                        pcfargs=list()) {
  
  dataname <- getdataname(deparse(substitute(X), width=20, nlines=1), ...)
  
  if(inherits(X, "fv")) {
    g <- X
    if(!(attr(g, "fname") %in% c("g", "g[inhom]")))
      warning("Argument X does not appear to be a pair correlation function")
  } else if(inherits(X, "ppp")) {
    g <- do.call("pcf.ppp", append(list(X), pcfargs))
    dataname <- paste("pcf(", dataname, ")", sep="")
    if(is.null(lambda))
      lambda <- summary(X)$intensity
  } else 
    stop("Unrecognised format for argument X")

  startpar <- check.named.vector(startpar, c("sigma2","alpha"))

  # digest parameters of Covariance model and test validity
  ph <- .Spatstat.ClusterModelInfo$LGCP$parhandler
  cmodel <- do.call(ph, covmodel)
  
  theoret <- .Spatstat.ClusterModelInfo$LGCP$pcf
  
  result <- mincontrast(g, theoret, startpar,
                        ctrl=list(q=q, p=p, rmin=rmin, rmax=rmax),
                        fvlab=list(label="%s[fit](r)",
                          desc="minimum contrast fit of LGCP"),
                        explain=list(dataname=dataname,
                          fname=attr(g, "fname"),
                          modelname="log-Gaussian Cox process"),
                        ...,
                        model=cmodel$model,
                        margs=cmodel$margs)
  # imbue with meaning
  par <- result$par
  names(par) <- c("sigma2", "alpha")
  result$par <- par
  result$covmodel <- cmodel
  # infer model parameters
  mu <- if(is.numeric(lambda) && length(lambda) == 1 && lambda > 0)
           log(lambda) - par[["sigma2"]]/2 else NA
  result$modelpar <- c(sigma2=par[["sigma2"]],
                       alpha=par[["alpha"]],
                       mu=mu)
  result$internal <- list(model="lgcp")
  return(result)
}


