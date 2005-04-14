#
#   envelope.R
#
#   computes simulation envelopes (finally)
#
#   $Revision: 1.4 $  $Date: 2005/04/19 01:04:51 $
#

envelope <- function(Y, fun=Kest, nsim=99, nrank=1, verbose=TRUE, ...) {
  if(inherits(Y, "ppm")) {
    model <- Y
    csr <- FALSE
  } else if(inherits(Y, "ppp")) {
    model <- ppm(Y, ~1, Poisson())
    csr <- TRUE
  } else
    stop("\`object\' should be a point process model or a point pattern")

  if(!is.function(fun)) 
    stop("unrecognised format for function \`fun\'")

  # evaluate function for data pattern
  X <- data.ppm(model)
  funX <- fun(X, ...)
  if(!inherits(funX, "fv"))
    stop("\`fun\' must return an object of class \"fv\"")
  argname <- attr(funX, "argu")
  valname <- attr(funX, "valu")
  rvals <- funX[[argname]]
  fX    <- funX[[valname]]

  # simulate
  if(verbose)
    cat(paste("Simulating", nsim, "realisations of",
              if(csr) "CSR" else "model", "...\n"))
  rmodel <- rmhmodel(model, verbose=FALSE)
  rstart <- rmhstart(n.start=data.ppm(model)$n)
  rcontr <- rmhcontrol(nrep=1e5, expand=1.5)

  simvals <- matrix(, nrow=nsim, ncol=length(rvals))

  for(i in 1:nsim) {
    Xsim <- rmh(rmodel, rstart, rcontr, verbose=FALSE)
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
                 else "mean of simulations from model",
                 "lower envelope of simulations from model",
                 "upper envelope of simulations from model"))
  return(result)
}

  
  
