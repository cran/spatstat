#
#  residppm.R
#
# computes residuals for fitted point process model
#
#
# $Revision: 1.4 $ $Date: 2006/06/02 08:36:24 $
#

residuals.ppm <- function(object, type="raw", ..., check=TRUE, 
                          fittedvalues = fitted.ppm(object, check=check)) {
  
  verifyclass(object, "ppm")

  if(check && missing(fittedvalues) && damaged.ppm(object)) 
    stop("object format corrupted; try update(object, use.internal=TRUE)")

  typetable <- c("inverse", "raw", "pearson", "Pearson")
  typemap <-   c("inverse", "raw", "pearson", "pearson")
  
  if(is.na(m <- pmatch(type, typetable)))
    stop(paste("Unrecognised choice of \'type\':", type))
  else
    type <- typemap[m]

  # Extract quadrature points and weights
  Q <- quad.ppm(object)
  U <- union.quad(Q) # quadrature points
  Z <- is.data(Q) # indicator data/dummy
  W <- w.quad(Q) # quadrature weights

  # Compute fitted conditional intensity at quadrature points
  lambda <- fittedvalues

  # indicator is 1 if lambda > 0
  # (adjusted for numerical behaviour of predict.glm)
  indicator <- (lambda > .Machine$double.eps)

  # Evaluate residual measure components
  discrete <- ifelse(Z,
                     switch(type,
                            raw     = 1,
                            inverse = 1/lambda,
                            pearson = 1/sqrt(lambda)
                            ),
                     0)
  continuous <- switch(type,
                       raw     = -lambda,
                       inverse = -indicator,
                       pearson = -indicator * sqrt(lambda))

  # Discretised residual measure (return value)
  res <- discrete + W * continuous

  # name the residuals
  attr(res, "type") <- type
  attr(res, "typename") <- paste(typetable[m], "residuals")

  # also give the components of the exact residual measure
  attr(res, "discrete") <- discrete
  attr(res, "continuous") <-  continuous
  attr(res, "atoms") <- as.logical(Z)
  
  return(res)
}

