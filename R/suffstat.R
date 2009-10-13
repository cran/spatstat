#
#   suffstat.R
#
# calculate sufficient statistic
#
#  $Revision: 1.12 $  $Date: 2009/10/07 02:36:51 $
#
#

suffstat <- function(model, X=data.ppm(model)) {
  cl <- sys.call()
  callstring <- paste(deparse(cl), collapse=" ")

  verifyclass(model, "ppm")
  if(!missing(X))
    verifyclass(X, "ppp")
  else
    X <- NULL

  ss <- model$interaction$family$suffstat

  func <- if(!is.null(ss))
            ss
          else if(summary(model, quick=TRUE)$poisson)
            suffstat.poisson
          else
            suffstat.generic

  return(func(model, X, callstring))
}

suffstat.generic <- function(model, X=NULL, callstring="suffstat.generic") {
  # This should work for an arbitrary ppm
  # since it uses the fundamental relation between
  # conditional intensity and likelihood.
  # But it is computationally intensive.

  verifyclass(model, "ppm")
  coefnames <- names(coef(model))

  if(is.null(X)) {
    X <- data.ppm(model)
    modelX <- model
  } else {
    verifyclass(X, "ppp")
    # refit the model to determine which points are used in pseudolikelihood
    modelX <- update(model, X, method="mpl")
  }
  
  # find data points which do not contribute to pseudolikelihood
  mplsubset <- getglmdata(modelX)$.mpl.SUBSET
  mpldata   <- is.data(quad.ppm(modelX))
  contribute <- mplsubset[mpldata]

  if(!any(contribute)) 
    # result is zero vector
    return(0 * coef(model))

  # Add points one-by-one
  # If there are points which don't contribute, condition on them
  use <- which(contribute)   
  dontuse <- which(!contribute)
  for(i in seq(length(use))) {
    prior <- if(i == 1) c() else use[1:(i-1)]
    prior <- c(dontuse, prior)
    Xprior <- X[prior]
    Xcurrent <- X[use[i]]
    mom <- partialModelMatrix(Xprior, Xcurrent, model, "suffstat")
    lastrow <- length(prior) + 1
    momrow <- mom[lastrow, ]
    names(momrow) <- coefnames
    if(i == 1)
      result <- momrow
    else
      result <- momrow + result
  }
  attr(result, "mplsubset") <- NULL
  return(result)
}

killinteraction <- function(model) {
  verifyclass(model, "ppm")
  ispoisson <- summary(model, quick=TRUE)$poisson
  if(ispoisson)
    return(model)
  # surgery required
  newmodel <- model
  newmodel$interaction <- NULL
  if(!is.null(Vnames <- model$internal$Vnames)) {
    matches <- names(model$coef) %in% Vnames
    newmodel$coef <- model$coef[!matches]
    newmodel$internal$Vnames <- NULL
  }
  # the other 'internal' stuff may still be wrong (or `preserved')
  return(newmodel)
}

suffstat.poisson <- function(model, X, callstring="suffstat.poisson") {
  verifyclass(model, "ppm")
  if(is.null(X))
    X <- data.ppm(model)
  else 
    verifyclass(X, "ppp")
  
  su <- summary(model, quick=TRUE)
  if(!(su$poisson))
    stop("Model is not a Poisson process")

  Empty <- X[numeric(0)]
  mom <- partialModelMatrix(X, Empty, model, "suffstat")

  nmom <- ncol(mom)
  ncoef <- length(coef(model))
  if(nmom != ncoef)
    stop("Internal error: number of columns of model matrix does not match number of coefficients in fitted model")
  
  if(nmom > 1 && any(colnames(mom) != names(coef(model))))
    warning("Internal error: mismatch between column names of model matrix and names of coefficient vector in fitted model")
     
  o1sum   <- apply(mom, 2, sum)
  return(o1sum)
}

