# 
#  fitted.ppm.R
#
# method for 'fitted' for ppm objects
#
#   $Revision: 1.2 $   $Date: 2006/01/09 07:09:02 $
# 

fitted.ppm <- function(object, ..., type="lambda", dataonly=FALSE) {
  verifyclass(object, "ppm")
  
  uniform <- is.poisson.ppm(object) && no.trend.ppm(object)

  typelist <- c("lambda", "cif",    "trend")
  typevalu <- c("lambda", "lambda", "trend")
  if(is.na(m <- pmatch(type, typelist)))
    stop(paste("Unrecognised choice of \`type\':", type))
  type <- typevalu[m]

  if(uniform) {
    fitcoef <- coef.ppm(object)
    lambda <- exp(fitcoef[[1]])
    Q <- quad.ppm(object)
    lambda <- rep(lambda, n.quad(Q))
  } else {
    glmdata <- object$internal$glmdata
    glmfit  <- getglmfit(object)
    Vnames <- object$internal$Vnames
    interacting <- !is.null(Vnames)
    
    # Modification of `glmdata' may be required
    if(interacting) 
      switch(type,
           trend={
             # zero the interaction statistics
             glmdata[ , Vnames] <- 0
           },
           lambda={
             # Find any dummy points with zero conditional intensity
             forbid <- matrowany(as.matrix(glmdata[, Vnames]) == -Inf)
             # exclude from predict.glm
             glmdata <- glmdata[!forbid, ]
           })

    # Compute predicted [conditional] intensity values
    lambda <- predict(glmfit, newdata=glmdata, type="response")
    # Note: the `newdata' argument is necessary in order to obtain
    # predictions at all quadrature points. If it is omitted then
    # we would only get predictions at the quadrature points j
    # where glmdata$SUBSET[j]=TRUE.

    if(interacting && type=="lambda") {
     # reinsert zeroes
      lam <- numeric(length(forbid))
      lam[forbid] <- 0
      lam[!forbid] <- lambda
      lambda <- lam
    }

  }
  if(dataonly)
    lambda <- lambda[is.data(quad.ppm(object))]
  
  return(lambda)
}


