#
#   envelope.R
#
#   computes simulation envelopes (finally)
#
#   $Revision: 1.15 $  $Date: 2005/09/08 10:13:10 $
#

envelope <- function(Y, fun=Kest, nsim=99, nrank=1, verbose=TRUE,
                     ..., simulate=NULL,
                     start=NULL, control=list(nrep=1e5, expand=1.5)) {

  # determine type of simulation
  if(!is.null(simulate)) {
    csr <- FALSE
    model <- NULL
    if(!is.expression(simulate))
      stop("\`simulate\' should be an expression")
    if(inherits(Y, "ppp"))
      X <- Y
    else if(inherits(Y, "ppm"))
      X <- data.ppm(Y)
    else stop("\`object\' should be a point process model or a point pattern")
  } else if(inherits(Y, "ppm")) {
    csr <- FALSE
    model <- Y
    X <- data.ppm(model)
  } else if(inherits(Y, "ppp")) {
    csr <- TRUE
    X <- Y
    rmhstuff <- !is.null(start) || !missing(control)
    sY <- summary(Y)
    Yintens <- sY$intensity
    Ywin <- Y$window
    Ymarx <- Y$marks
    if(!is.marked(Y)) {
      # unmarked point pattern
      if(rmhstuff)
        model <- ppm(Y, ~1, Poisson())
      else {
        model <- NULL
        simulate <- expression(rpoispp(Yintens, win=Ywin))
      }
    } else {
      if(rmhstuff) {
        if(sY$is.multitype)
          model <- ppm(Y, ~marks, Poisson())
        else
          stop("Sorry, ppm is not implemented for marked point patterns where the marks are not a factor")
      } else {
        model <- NULL
        simulate <- expression({A <- rpoispp(Yintens, win=Ywin);
                                A %mark% sample(Ymarx, A$n, replace=TRUE)})
      }
    }
  } else 
  stop("\`object\' should be a point process model or a point pattern")

  metrop <- !is.null(model)

  # Name of function, for error messages
  fname <- if(is.name(substitute(fun))) deparse(substitute(fun))
           else if(is.character(fun)) fun else "fun"
  fname <- sQuote(fname)
  
  # validate other arguments
  if(is.character(fun))
    fun <- get(fun)
  if(!is.function(fun)) 
    stop(paste("unrecognised format for function", fname))
  if(!any(c("r", "...") %in% names(formals(fun))))
    stop(paste(fname, "should have an argument \`r\' or \`...\'"))
  
  if((nrank %% 1) != 0)
    stop("nrank must be an integer")
  stopifnot(nrank > 0 && nrank < nsim/2)

  # evaluate function for data pattern X
  funX <- fun(X, ...)
  if(!inherits(funX, "fv"))
    stop(paste(fname, "must return an object of class", sQuote("fv")))
  argname <- attr(funX, "argu")
  valname <- attr(funX, "valu")
  rvals <- funX[[argname]]
  fX    <- funX[[valname]]

  # initialise simulations
  if(metrop) {
    rmodel <- rmhmodel(model, verbose=FALSE)
    if(is.null(start))
      start <- list(n.start=X$n)
    rstart <- rmhstart(start)
    rcontr <- rmhcontrol(control)
  }

  nrvals <- length(rvals)
  simvals <- matrix(, nrow=nsim, ncol=nrvals)

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
        stop(paste("Evaluating the expression", sQuote("simulate"),
                   "did not yield a point pattern"))
    }
    # apply function
    funXsim <- fun(Xsim, r=rvals, ...)
    # sanity checks
    if(!inherits(funXsim, "fv"))
      stop(paste("When applied to a simulated pattern, the function",
                 fname, "did not return an object of class",
                 sQuote("fv")))
    argname.sim <- attr(funXsim, "argu")
    valname.sim <- attr(funXsim, "valu")
    if(argname.sim != argname)
      stop(paste("The objects returned by", fname,
                 "when applied to a simulated pattern",
                 "and to the data pattern",
                 "are incompatible. They have different argument names",
                 sQuote(argname.sim), "and", sQuote(argname), 
                 "respectively"))
    if(valname.sim != valname)
      stop(paste("When", fname, "is applied to a simulated pattern",
                 "it provides an estimate named", sQuote(valname.sim), 
                 "whereas the estimate for the data pattern is named",
                 sQuote(valname),
                 ". Try using the argument", sQuote("correction"),
                 "to make them compatible"))

    # extract the values for simulation i
    simvals.i <- funXsim[[valname]]
    if(length(simvals.i) != nrvals)
      stop("Vectors of function values have incompatible lengths")
      
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

  
  
