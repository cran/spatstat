#
#  residppm.R
#
# computes residuals for fitted point process model
#
#
# $Revision: 1.13 $ $Date: 2011/01/15 03:14:47 $
#

residuals.ppm <- function(object, type="raw", ..., check=TRUE, drop=FALSE,
                 fittedvalues = fitted.ppm(object, check=check, drop=drop)) {
  
  verifyclass(object, "ppm")

  if(check && missing(fittedvalues) && damaged.ppm(object)) 
    stop("object format corrupted; try update(object, use.internal=TRUE)")

  type <- pickoption("type", type,
                     c(inverse="inverse",
                       raw="raw",
                       pearson="pearson",
                       Pearson="pearson",
                       score="score"))
  typenames <- c(inverse="inverse-lambda residuals",
                 raw="raw residuals",
                 pearson="Pearson residuals",
                 score="score residuals")
  typename <- typenames[[type]]

  # Extract quadrature points and weights
  Q <- quad.ppm(object, drop=drop)
  U <- union.quad(Q) # quadrature points
  Z <- is.data(Q) # indicator data/dummy
#  W <- w.quad(Q) # quadrature weights

  # Compute fitted conditional intensity at quadrature points
  lambda <- fittedvalues

  # indicator is 1 if lambda > 0
  # (adjusted for numerical behaviour of predict.glm)
  indicator <- (lambda > .Machine$double.eps)

  if(type == "score") {
    # need the covariates
    X <- model.matrix(object)
    if(drop) {
      gs <- getglmsubset(object)
      ok <- !is.na(gs) && gs
      X <- X[ok,]
    }
  }
      
  # Evaluate residual measure components

  discrete <- switch(type,
                     raw     = rep(1, sum(Z)), 
                     inverse = 1/lambda[Z],
                     pearson = 1/sqrt(lambda[Z]),
                     score   = X[Z, ]
                     )

  continuous <- switch(type,
                       raw     = -lambda,
                       inverse = -indicator,
                       pearson = -indicator * sqrt(lambda),
                       score   = -lambda * X)

  # Residual measure (return value)
  res <- msr(Q, discrete, continuous)

  # name the residuals
  attr(res, "type") <- type
  attr(res, "typename") <- typename

  return(res)
}

