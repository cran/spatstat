#
#   envelope.R
#
#   computes simulation envelopes 
#
#   $Revision: 1.25 $  $Date: 2006/04/13 12:36:10 $
#

envelope <- function(Y, fun=Kest, nsim=99, nrank=1, verbose=TRUE,
                     ..., simulate=NULL, clipdata=TRUE,
                     start=NULL, control=list(nrep=1e5, expand=1.5),
                     transform=NULL, global=FALSE, ginterval=NULL) {
  cl <- match.call()
  Yname <- deparse(substitute(Y))
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
    sY <- summary(Y, checkdup=FALSE)
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

  # set up simulation parameters
  metrop <- !is.null(model)
  if(metrop) {
    rmodel <- rmhmodel(model, verbose=FALSE)
    if(is.null(start))
      start <- list(n.start=X$n)
    rstart <- rmhstart(start)
    rcontr <- rmhcontrol(control)
  }

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

  rgiven <- "r" %in% names(list(...))

  if(tran <- !is.null(transform)) {
    transform.funX <- eval(substitute(substitute(
                              eval(transform),
                              list("."=as.name("funX")))))
    transform.funXsim <- eval(substitute(substitute(
                              eval(transform),
                              list("."=as.name("funXsim")))))
  }
  if(!is.null(ginterval)) 
    stopifnot(is.numeric(ginterval) && length(ginterval) == 2)

  # clip the data?
  if(clipdata) {
    # Generate one realisation
    if(metrop)
      Xsim <- rmh(rmodel, rstart, rcontr, verbose=FALSE)
    else {
      Xsim <- eval(simulate)
      if(!inherits(Xsim, "ppp"))
        stop(paste("Evaluating the expression", sQuote("simulate"),
                   "did not yield a point pattern"))
    }
    # Extract window
    clipwin <- Xsim$window
    if(!is.subset.owin(clipwin, X$window))
      warning("Window containing simulated patterns is not a subset of data window")
  }
  
  
  # evaluate function for data pattern X 
  if(!clipdata)
    funX <- fun(X, ...)
  else
    funX <- fun(X[,clipwin], ...)
    
  if(!inherits(funX, "fv"))
    stop(paste(fname, "must return an object of class", sQuote("fv")))

  argname <- attr(funX, "argu")
  valname <- attr(funX, "valu")

  if(tran) {
    # extract only the recommended value
    if(csr) 
      funX <- funX[, c(argname, valname, "theo")]
    else
      funX <- funX[, c(argname, valname)]
    # apply the transformation to it
    funX <- do.call("eval.fv", list(transform.funX))
  }

  rvals <- funX[[argname]]
  fX    <- funX[[valname]]


  # default domain over which to maximise
  alim <- attr(funX, "alim")
  if(global && is.null(ginterval))
    ginterval <- alim
  
  ######### simulate #######################
  if(verbose)
    cat(paste("Generating", nsim, "simulations",
              if(csr) "of CSR" else if(metrop) "of model" else NULL,
              "...\n"))
  nrvals <- length(rvals)
  simvals <- matrix(, nrow=nsim, ncol=nrvals)

  # start simulation loop 
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
    if(rgiven) 
      funXsim <- fun(Xsim, ...)
    else
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

    if(tran) {
      # extract only the recommended value
      if(csr) 
        funXsim <- funXsim[, c(argname, valname, "theo")]
      else
        funXsim <- funXsim[, c(argname, valname)]
      # apply the transformation to it
      funXsim <- do.call("eval.fv", list(transform.funXsim))
    }
      
    # extract the values for simulation i
    simvals.i <- funXsim[[valname]]
    if(length(simvals.i) != nrvals)
      stop("Vectors of function values have incompatible lengths")
      
    simvals[i,] <- funXsim[[valname]]
    if(verbose)
      cat(paste(i, if(i < nsim) ", " else ".",
                if(i %% 10 == 0) "\n", sep=""))
  }
  ##  end simulation loop
  
  if(verbose)
    cat("\nDone.\n")

  # compute envelopes
  orderstat <- function(x, n) { sort(x)[n] }
  if(!global) {
    # POINTWISE ENVELOPES
    lo <- apply(simvals, 2, orderstat, n=nrank)
    hi <- apply(simvals, 2, orderstat, n=nsim-nrank+1)
    lo.name <- paste("lower pointwise envelope of simulations")
    hi.name <- paste("upper pointwise envelope of simulations")

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
  } else {
    # SIMULTANEOUS ENVELOPES
    domain <- (rvals >= ginterval[1]) & (rvals <= ginterval[2])
    funX <- funX[domain, ]
    simvals <- simvals[ , domain]
    if(csr)
      theory <- funX[["theo"]]
    else
      theory <- m <- apply(simvals, 2, mean, na.rm=TRUE)
    # compute max absolute deviations
    deviations <- sweep(simvals, 2, theory)
    suprema <- apply(abs(deviations),  1, max, na.rm=TRUE)
    # ranked deviations
    dmax <- orderstat(suprema, n=nsim-nrank+1)
    # simultaneous bands
    lo <- theory - dmax
    hi <- theory + dmax
    lo.name <- "lower critical boundary"
    hi.name <- "upper critical boundary"

    if(csr)
      results <- data.frame(r=rvals[domain],
                            obs=fX[domain],
                            theo=theory,
                            lo=lo,
                            hi=hi)
    else
      results <- data.frame(r=rvals[domain],
                            obs=fX[domain],
                            mmean=m,
                            lo=lo,
                            hi=hi)
  }
  
  result <- fv(results,
               argu="r",
               ylab=attr(funX, "ylab"),
               valu="obs",
               fmla= deparse(. ~ r),
               alim=attr(funX, "alim"),
               labl=c("r", "obs(r)", if(csr) "theo(r)" else "mean(r)",
                 "lo(r)", "hi(r)"),
               desc=c("distance argument r",
                 "function value for data pattern",
                 if(csr) "theoretical value for CSR"
                 else "mean of simulations",
                 lo.name, hi.name))
  attr(result, "dotnames") <- c("obs",
                                if(csr) "theo" else "mmean",
                                "hi", "lo")
  class(result) <- c("envelope", class(result))
  attr(result, "einfo") <- list(call=cl,
                                Yname=Yname,
                                csr=csr,
                                nrank=nrank,
                                nsim=nsim,
                                global=global)
  return(result)
}


print.envelope <- function(x, ...) {
  e <- attr(x, "einfo")
  g <- e$global
  nr <- e$nrank
  nsim <- e$nsim
  csr <- e$csr
  fname <- deparse(attr(x, "ylab"))
  type <- if(g) "Simultaneous" else "Pointwise"
  cat(paste(type, "critical envelopes for", fname, "\n"))
  cat(paste("Obtained from", nsim,
            "simulations of", if(csr) "CSR" else "fitted model", "\n"))
  alpha <- if(g) { nr/(nsim+1) } else { 2 * nr/(nsim+1) }
  cat(paste("Significance level of",
            if(!g) "pointwise",
            "Monte Carlo test:",
            paste(if(g) nr else 2 * nr,
                  "/", nsim+1, sep=""),
            "=", alpha, "\n"))
  cat(paste("Data:", e$Yname, "\n\n"))
  NextMethod("print")
}
                  
summary.envelope <- function(object, ...) {
  e <- attr(object, "einfo")
  g <- e$global
  nr <- e$nrank
  nsim <- e$nsim
  csr <- e$csr
  fname <- deparse(attr(object, "ylab"))
  type <- if(g) "Simultaneous" else "Pointwise"
  cat(paste(type, "critical envelopes for", fname, "\n"))
  cat(paste("Obtained from", nsim,
            "simulations of", if(csr) "CSR" else "fitted model", "\n"))
  ordinal <- function(k) {
    last <- k %% 10
    ending <- if(last == 1) "st" else
              if(last == 2) "nd" else
              if(last== 3) "rd" else "th"
    return(paste(k, ending, sep=""))
  }
  lo.ord <- if(nr == 1) "minimum" else paste(ordinal(nr), "smallest")
  hi.ord <- if(nr == 1) "maximum" else paste(ordinal(nsim-nr+1), "largest")
  if(g) 
    cat(paste("Envelopes computed as",
              if(csr) "theoretical curve" else "mean of simulations",
              "plus/minus", hi.ord,
              "simulated value of maximum absolute deviation\n"))
  else {
    cat(paste("Upper envelope: pointwise", hi.ord,"of simulated curves\n"))
    cat(paste("Lower envelope: pointwise", lo.ord,"of simulated curves\n"))
  }
  alpha <- if(g) { nr/(nsim+1) } else { 2 * nr/(nsim+1) }
  cat(paste("Significance level of Monte Carlo test:",
            paste(if(g) nr else 2 * nr,
                  "/", nsim+1, sep=""),
            "=", alpha, "\n"))
  cat(paste("Data:", e$Yname, "\n\n"))
  return(invisible(NULL))
}
  
  
  
  
