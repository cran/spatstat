#
#  residppm.R
#
# computes residuals for fitted point process model
#
#
# $Revision: 1.8 $ $Date: 2010/11/01 03:31:48 $
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
  
  # Extract quadrature points and weights
  Q <- quad.ppm(object, drop=drop)
  U <- union.quad(Q) # quadrature points
  Z <- is.data(Q) # indicator data/dummy
  W <- w.quad(Q) # quadrature weights

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
                     raw     = as.integer(Z),
                     inverse = ifelse(Z, 1/lambda, 0),
                     pearson = ifelse(Z, 1/sqrt(lambda), 0),
                     score   = Z * X
                     )

  continuous <- switch(type,
                       raw     = -lambda,
                       inverse = -indicator,
                       pearson = -indicator * sqrt(lambda),
                       score   = -lambda * X)

  # Discretised residual measure (return value)
  res <- discrete + W * continuous

  # name the residuals
  attr(res, "type") <- type
  attr(res, "typename") <- paste(if(type == "pearson") "Pearson" else type,
                                 "residuals")

  # also give the components of the exact residual measure
  attr(res, "discrete") <- discrete
  attr(res, "continuous") <-  continuous
  attr(res, "atoms") <- as.logical(Z)
  
  return(res)
}

