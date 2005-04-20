#
#	ppmclass.R
#
#	Class 'ppm' representing fitted point process models.
#
#
#	$Revision: 2.4 $	$Date: 2005/04/16 01:40:57 $
#
#       An object of class 'ppm' contains the following:
#
#            $method           model-fitting method (currently "mpl")
#
#            $coef             vector of fitted regular parameters
#                              as given by coef(glm(....))
#
#            $theta            vector of fitted regular parameters
#                              as given by dummy.coef(glm(....))
#
#            $trend            the trend formula
#                              or NULL 
#
#            $interaction      the interaction family 
#                              (an object of class 'interact') or NULL
#
#            $Q                the quadrature scheme used
#
#            $maxlogpl         the maximised value of log pseudolikelihood
#
#            $internal         list of internal calculation results
#
#            $correction       name of edge correction method used
#            $rbord            erosion distance for border correction (or NULL)
#
#            $the.call         the originating call to ppm()
#
#            $the.version      version of mpl() which yielded the fit
#
#
#------------------------------------------------------------------------

is.ppm <- function(x) { inherits(x, "ppm") }

print.ppm <- function(x, ...) {
	verifyclass(x, "ppm")

        s <- summary.ppm(x)
        
        notrend <-    s$no.trend
	stationary <- s$stationary
	poisson <-    s$poisson
        markeddata <- s$marked
        multitype  <- s$multitype
        
        markedpoisson <- poisson && markeddata

        # names of interaction variables if any
        Vnames <- s$Vnames
        # their fitted coefficients
        theta <- s$theta

        # ----------- Print model type -------------------
        
	cat(s$name)
        cat("\n")
        
        if(markeddata) mrk <- s$entries$marks
        if(multitype) {
            cat("Possible marks: \n")
            cat(paste(levels(mrk)))
          }

        # ----- trend --------------------------

        cat(paste("\n", s$trend$name, ":\n", sep=""))

	if(!notrend) {
		cat("Trend formula: ")
		print(s$trend$formula)
        }
        
        cat(paste("\n", s$trend$label, ":\n", sep=""))

        tv <- s$trend$value
        if(!is.list(tv))
          print(tv)
        else 
          for(i in seq(tv))
            print(tv[[i]])
        
        # ---- Interaction ----------------------------

	if(!poisson) {
          cat("\nInteraction:\n")
          print(s$entries$interaction)
        
          cat(paste(s$interaction$header, ":\n", sep=""))
          print(s$interaction$printable)
        }

	invisible(NULL)
}

quad.ppm <- function(object) {
  verifyclass(object, "ppm")
  object$Q
}

data.ppm <- function(object) { 
  verifyclass(object, "ppm")
  object$Q$data
}

dummy.ppm <- function(object) { 
  verifyclass(object, "ppm")
  object$Q$dummy
}
  
# method for 'coef'
coef.ppm <- function(object, ...) {
  verifyclass(object, "ppm")
  object$coef
}

# method for 'fitted'
fitted.ppm <- function(object, ..., type="lambda") {
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

    if(type == "trend") {
      # first zero the interaction statistics
      Vnames <- object$internal$Vnames
      glmdata[ , Vnames] <- 0
    }
    lambda <- predict(glmfit, newdata=glmdata, type="response")
    # Note: the `newdata' argument is necessary in order to obtain predictions
    # at all quadrature points. If it is omitted then we would only get
    # predictions at the quadrature points j where glmdata$SUBSET[j]=TRUE.
  }
  
  return(lambda)
}


getglmfit <- function(object) {
  verifyclass(object, "ppm")
  glmfit <- object$internal$glmfit
  if(is.null(glmfit))
    return(NULL)
  if(object$method == "ho")
    glmfit$coefficients <- object$coef
  return(glmfit)
}

    
# ??? method for 'effects' ???

