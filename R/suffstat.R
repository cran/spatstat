#
#   suffstat.R
#
# calculate sufficient statistic
#
#  $Revision: 1.1 $  $Date: 2005/04/28 00:47:30 $
#
#

suffstat <- function(model, X) {
  verifyclass(model, "ppm")
  verifyclass(X, "ppp")

  ss <- model$interaction$family$suffstat

  func <- if(!is.null(ss))
            ss
          else if(summary(model, quick=TRUE)$poisson)
            suffstat.poisson
          else
            suffstat.generic

  return(func(model, X))
}

suffstat.generic <- function(model, X) {
  # This should work for an arbitrary ppm
  # since it uses the fundamental relation between
  # conditional intensity and likelihood.
  # But it is computationally intensive.
  
  verifyclass(model, "ppm")
  verifyclass(X, "ppp")

  shazzam <- function(Q, X, P, model) {
    trend <- model$trend
    inter <- model$interaction
    covar <- model$covariates
    prep  <- mpl.prepare(Q, X, P, trend, inter, covar,
                       correction=model$correction,
                       rbord=model$rbord)
    fmla    <- prep$fmla
    glmdata <- prep$glmdata
    mof <- model.frame(fmla, glmdata)
    mom <- model.matrix(fmla, mof)

    if(any(colnames(mom) != names(coef(model))))
      warning("Internal error: mismatch between column names of model matrix and names of coefficient vector in fitted model")

    dummy <- !is.data(Q)
    mom <- mom[dummy, ]
    return(mom)
  }
  
  n <- X$n

  # for i = 2, ..., n compute T(X[i], X[1:(i-1)])
  # where conditional intensity lambda(u,X) = exp(theta T(u,X))

  for(i in 2:n) {
    Xprior <- X[1:(i-1)]
    Q <- quad(Xprior, X[i])
    mom <- shazzam(Q, Xprior, X[1:i], model)
    if(i == 2)
      result <- mom
    else
      result <- mom + result
  }

  # to compute T(X[1], Empty) we know there are no interaction contributions
  # so evaluate T for the corresponding Poisson model.
  # This requires complicated surgery...
  
  poismodel <- killinteraction(model)
  Empty <- X[rep(FALSE, n)]
  Q <- quad(Empty, X[1])
  mom1  <- shazzam(Q, Empty, X[1], poismodel)

  # add to result
  witch <- match(names(mom1), names(result))
  if(any(is.na(witch)))
    stop("Internal error: names in zero term do not match names in other terms")
  result[witch] <- result[witch] + mom1

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
  # the other 'internal' stuff may still be wrong
  return(newmodel)
}

suffstat.poisson <- function(model, X) {
  verifyclass(model, "ppm")
  verifyclass(X, "ppp")
  
  su <- summary(model, quick=TRUE)
  if(!(su$poisson))
    stop("Model is not a Poisson process")

  trend <- model$trend
  covar <- model$covariates
 
  Empty <- X[rep(FALSE, X$n)]
  Q     <- quad(X, Empty)
  prep  <- mpl.prepare(Q, X, X, trend, Poisson(), covar,
                       correction=model$correction,
                       rbord=model$rbord)
  fmla    <- prep$fmla
  glmdata <- prep$glmdata

  mof <- model.frame(fmla, glmdata)
  mom <- model.matrix(fmla, mof)

  nmom <- ncol(mom)
  ncoef <- length(coef(model))
  if(nmom != ncoef)
    stop("Internal error: number of columns of model matrix does not match number of coefficients in fitted model")
  
  if(nmom > 1 && any(colnames(mom) != names(coef(model))))
    warning("Internal error: mismatch between column names of model matrix and names of coefficient vector in fitted model")
     
  o1sum   <- apply(mom, 2, sum)
  return(o1sum)
}

