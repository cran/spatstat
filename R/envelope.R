#
#   envelope.R
#
#   computes simulation envelopes (finally)
#
#   $Revision: 1.10 $  $Date: 2005/06/08 17:44:08 $
#

envelope <- function(Y, fun=Kest, nsim=99, nrank=1, verbose=TRUE,
                     ..., simulate=NULL,
                     start=NULL, control=list(nrep=1e5, expand=1.5)) {

  # determine type of simulation
  if(inherits(Y, "ppm")) {
    model <- Y
    X <- data.ppm(model)
    csr <- FALSE
  } else if(inherits(Y, "ppp")) {
    X <- Y
    model <- ppm(Y, ~1, Poisson())
    csr <- TRUE
  } else 
    stop("\`object\' should be a point process model or a point pattern")

  metrop <- is.null(simulate)
  if(!metrop) {
    csr <- FALSE
    model <- NULL
    if(!is.expression(simulate))
      stop("\`simulate\' should be an expression")
  }

  # validate other arguments
  if(is.character(fun))
    fun <- get(fun)
  if(!is.function(fun)) 
    stop("unrecognised format for function \`fun\'")

  if((nrank %% 1) != 0)
    stop("nrank must be an integer")
  stopifnot(nrank > 0 && nrank < nsim/2)

  # evaluate function for data pattern X
  funX <- fun(X, ...)
  if(!inherits(funX, "fv"))
    stop("\`fun\' must return an object of class \"fv\"")
  argname <- attr(funX, "argu")
  valname <- attr(funX, "valu")
  rvals <- funX[[argname]]
  fX    <- funX[[valname]]

  # accelerate the simulations for the uniform Poisson process
  if(csr && is.null(start) && is.null(control)) {
    metrop <- FALSE
    simulate <- expression(runifpoispp(lambda, win))
    lambda <- summary(X)$intensity
    win <- X$window
  }
  
  # initialise simulations
  if(metrop) {
    rmodel <- rmhmodel(model, verbose=FALSE)
    if(is.null(start))
      start <- list(n.start=X$n)
    rstart <- rmhstart(start)
    rcontr <- rmhcontrol(control)
  }
  
  simvals <- matrix(, nrow=nsim, ncol=length(rvals))

  # simulate
  if(verbose)
    cat(paste("Generating", nsim, "simulations",
              if(csr) "of CSR" else if(metrop) "of model" else NULL,
              "...\n"))

  for(i in 1:nsim) {
    if(metrop)
      Xsim <- rmh(rmodel, rstart, rcontr, verbose=FALSE)
    else {
      Xsim <- eval(simulate)
      if(!inherits(Xsim, "ppp"))
        stop("Evaluating the expression \`simulate\' did not yield a point pattern")
    }
    funXsim <- fun(Xsim, r=rvals, ...)
    simvals[i,] <- funXsim[[valname]]
    if(verbose)
      cat(paste(i, if(i < nsim) ", " else ".",
                if(i %% 10 == 0) "\n", sep=""))
  }
  if(verbose)
    cat("\nDone.\n")

  # compute order statistics
  orderstat <- function(x, n) { sort(x)[n] }
  lo <- apply(simvals, 2, orderstat, n=nrank)
  hi <- apply(simvals, 2, orderstat, n=nsim-nrank+1)
  if(csr) {
    results <- data.frame(r=rvals,
                          obs=fX,
                          theo=funX[["theo"]],
                          lo=lo,
                          hi=hi)
  } else {
    m <- apply(simvals, 2, mean, na.rm=TRUE)
    results <- data.frame(r=rvals,
                          obs=fX,
                          mmean=m,
                          lo=lo,
                          hi=hi)
  }
      
  result <- fv(results,
               argu="r",
               ylab=attr(funX, "ylab"),
               valu="obs",
               fmla=
               if(csr) deparse(cbind(obs, theo, lo, hi) ~ r)
               else deparse(cbind(obs, mmean, lo, hi) ~ r),
               alim=attr(funX, "alim"),
               labl=c("r", "obs(r)", if(csr) "theo(r)" else "mean(r)",
                 "lo(r)", "hi(r)"),
               desc=c("distance argument r",
                 "function value for data pattern",
                 if(csr) "theoretical value for CSR"
                 else "mean of simulations",
                 "lower envelope of simulations",
                 "upper envelope of simulations"))
  return(result)
}

  
  
