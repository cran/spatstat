#
#      alltypes.R
#
#   $Revision: 1.18 $   $Date: 2008/07/29 03:16:17 $
#
#
alltypes <- function(X, fun="K", ...,
                     dataname=NULL,verb=FALSE,envelope=FALSE) {
#
# Function 'alltypes' --- calculates a summary function for
# each type, or each pair of types, in a multitype point pattern
#
  verifyclass(X,"ppp")
  if(is.null(dataname))
    dataname <- deparse(substitute(X))

# --------------------------------------------------------------------  
# First inspect marks

  if(!is.marked(X)) {
    nmarks <- 0
    marklabels <- ""
  } else {
    mks <- marks(X)
    if(!is.factor(mks))
      stop("the marks must be a factor")
    ma <- levels(mks)
    nmarks <- length(ma)
    marklabels <- paste(ma)
  }

# ---------------------------------------------------------------------
# determine function name

  f.is.name <- is.name(substitute(fun))
  fname <-
    if(f.is.name)
      paste(as.name(substitute(fun)))
    else if(is.character(fun))
      fun
    else sQuote("fun") 

# ---------------------------------------------------------------------
# determine function to be called
  
  estimator <- 
    if(is.function(fun))
      fun
    else if(is.character(fun) && fun %in% c("F", "G", "J", "K", "pcf")) {
    # conventional abbreviations
      if(nmarks > 0) 
        switch(fun,
               F=Fest,
               G=Gcross,
               J=Jcross,
               K=Kcross,
               pcf=pcfcross)
      else
        switch(fun,
               F=Fest,
               G=Gest,
               J=Jest,
               K=Kest,
               pcf=pcf)
    } else if(is.character(fun))
      get(fun, mode="function")
    else 
      stop(paste(sQuote("fun"), "should be a function or a character string"))
  
# ------------------------------------------------------------------  
# determine how the function shall be called.
#
  indices.expected <- sum(c("i", "j") %in% names(formals(estimator)))

  apply.to.split   <- (indices.expected == 0 && nmarks > 1)
  if(apply.to.split)
    ppsplit <- split(X)
  
# --------------------------------------------------------------------  
# determine array dimensions and margin labels
  witch <-
    if(nmarks == 0)
      matrix(1, nrow=1, ncol=1, dimnames=list("",""))
    else if (nmarks == 1) 
      matrix(1, nrow=1, ncol=1, dimnames=list(marklabels, marklabels))
    else if(indices.expected != 2)
      matrix(1:nmarks, nrow=nmarks, ncol=1,
             dimnames=list(marklabels, ""))
    else 
      matrix(1:(nmarks^2),ncol=nmarks,nrow=nmarks, byrow=TRUE,
             dimnames <- list(marklabels, marklabels))

  # ------------ start computing -------------------------------  
  # if computing envelopes, first generate simulated patterns
  # using undocumented feature of envelope()
  if(envelope) {
    L <- envelope(X, fun=NULL, ..., verbose=verb,
                  internal=list(patterns=TRUE))
    intern <- attr(L, "internal")
  }

  # compute function array and build up 'fasp' object
  fns  <- list()
  deform <- list()
  k   <- 0

  for(i in 1:nrow(witch)) {
    Y <- if(apply.to.split) ppsplit[[i]] else X
    for(j in 1:ncol(witch)) {
      if(verb) cat("i =",i,"j =",j,"\n")
      k <- k+1
      fns[[k]] <- currentfv <- 
        if(!envelope) 
          switch(1+indices.expected,
                 estimator(Y, ...),
                 estimator(Y, i=ma[i], ...),
                 estimator(Y, i=ma[i], j=ma[j], ...))
        else
          switch(1+indices.expected,
                 envelope(Y, estimator, ...,
                          verbose=FALSE, simulate=L, internal=intern),
                 envelope(Y, estimator, i=ma[i], ...,
                          verbose=FALSE, simulate=L, internal=intern),
                 envelope(Y, estimator, i=ma[i], j=ma[j], ...,
                          verbose=FALSE, simulate=L, internal=intern)
                 )
      deform[[k]] <- attr(currentfv, "fmla")
    }
  }

  # wrap up into 'fasp' object
  title <- paste(if(nmarks > 1) "array of " else NULL,
                 if(envelope) "envelopes of " else NULL,
                 fname,
                 if(nmarks <= 1) " function " else " functions ",
                 "for ", dataname, ".", sep="")
  
  rslt <- fasp(fns, which=witch,
               formulae=deform,
               dataname=dataname,
               title=title)
  return(rslt)
}
