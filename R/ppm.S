#
#	ppm.S
#
#	Class 'ppm' representing fitted point process models.
#
#
#	$Revision: 1.20 $	$Date: 2003/05/02 11:24:07 $
#
#       An object of class 'ppm' contains the following:
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
#            $the.call         the originating call to mpl()
#
#            $the.version      version of mpl() which yielded the fit
#
#
#------------------------------------------------------------------------

print.ppm <- function(x, ...) {
	verifyclass(x, "ppm")

        notrend <- no.trend.ppm(x)
	stationary <- is.stationary.ppm(x)
	poisson <- is.poisson.ppm(x)
        markeddata <- is.marked(x$Q$data)
        
        markedpoisson <- poisson && markeddata
        # i.e. if the data are marked and mpl() is called with
        # a Poisson interaction, the fitted model is a marked Poisson process.

        # names of interaction variables if any
        Vnames <- x$internal$Vnames

        # ----------- Print model type -------------------
        
	cat(paste(
		if(stationary) "Stationary " else "Nonstationary ",
                if(markedpoisson) "marked " else "", 
                if(poisson) "Poisson process" else x$interaction$name,
                "\n", sep=""))

        if(markeddata) {
          mrk <- x$Q$data$marks
          if(is.factor(mrk)) {
            cat("Possible marks: \n")
            print(factor(levels(mrk)))
          }
        }

        if(exists("is.R") && is.R()) 
          theta <- x$coef # result of coef(glm(...))
        else
          theta <- x$theta # result of dummy.coef(glm(....))
          
        # ----- trivial case: uniform Poisson --------
        
	if(poisson && notrend) {
          lambda <- exp(theta[[1]])
          cat(paste("Uniform intensity", lambda))
          if(markeddata) {
            cat("(for each mark level).\n")
            m <- markspace.integral(x$Q$data)
            cat(paste("Total point intensity", lambda * m))
          } 
          cat("\n")
          return(invisible(NULL))
	}

        # ----- trend --------------------------

        cat(paste("\n\n-----------",
                  if(poisson) "Intensity" else "First order term",
                  "-------------\n\n"))

	if(!notrend) {
		cat("Trend formula:\n")
		print(x$trend)
        }

        # process is at least one of: marked, nonstationary, non-poisson
        
        if(stationary) {
          if(!markeddata) {
            # it can't be poisson
            cat("First order term:\n")
            beta <- exp(theta[[1]])
            print(unlist(list(beta=beta)))
          } else {
            # marked and stationary
            cat( if(poisson)
                "Fitted intensities:\n"
                 else
                "Fitted first order terms:\n")
            lev <- factor(levels(mrk))
            betas <- predict(x, newdata=data.frame(marks=lev), type="trend")
            names(betas) <- paste("beta_", as.character(lev), sep="")
            print(betas)
          }
        } else {
          # not stationary 
          cat("Fitted coefficients for trend formula:\n")
          # extract trend terms without trying to understand them much
          if(is.null(Vnames)) 
            trendbits <- theta
          else {
            agree <- outer(names(theta), Vnames, "==")
            whichbits <- apply(!agree, 1, all)
            trendbits <- theta[whichbits]
          }
          # decide whether there are 'labels within labels'
          unlabelled <- unlist(lapply(trendbits,
                                    function(x) { is.null(names(x)) } ))
          if(all(unlabelled))
            print(unlist(trendbits))
          else 
            for(i in seq(trendbits))
              if(unlabelled[i])
                print(unlist(trendbits[i]))
              else
                print(trendbits[[i]])
        }
        
        # ---- Interaction ----------------------------

	if(!poisson) {
          cat("\n\n-----------Interaction------------------\n\n")
          print(x$interaction)
        
          if(!is.null(x$interaction$interpret)) {
            sensible <- (x$interaction$interpret)(x$coef, x$interaction)
            cat(paste("Fitted ", sensible$inames, ":\n", sep=""))
            print(sensible$printable)
          } else {
            # panic - fallback
            cat("Fitted interaction terms:\n")
            print(exp(unlist(theta[x$internal$Vnames])))
          }
        }

        cat("\n\n----------- gory details -----\n")
        
        cat("\nFitted regular parameters (theta): \n")
        print(theta)

        cat("\nFitted exp(theta): \n")
        print(exp(unlist(theta)))

	invisible(NULL)
}

no.trend.ppm <- function(x) {
  # the result is TRUE only when the trend is ~1
  verifyclass(x, "ppm")
  length(termsinformula(x$trend)) == 0
}

is.stationary.ppm <- function(x) {
  verifyclass(x, "ppm")
  turms <- termsinformula(x$trend)
  (length(turms) == 0) ||
     ((length(turms) == 1) && turms[1] == "marks")
}

is.poisson.ppm <- function(x) {
  verifyclass(x, "ppm")
  is.null(x$interaction) || is.null(x$interaction$pot)
}

is.marked.ppm <- function(X, ...) {
  verifyclass(X, "ppm")
  is.marked.ppp(X$Q$data)
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
fitted.ppm <- function(object, ...) {
  verifyclass(object, "ppm")
  
  uniform <- is.poisson.ppm(object) && no.trend.ppm(object)

  if(uniform) {
    fitcoef <- coef.ppm(object)
    lambda <- exp(fitcoef[[1]])
    Q <- quad.ppm(object)
    lambda <- rep(lambda, n.quad(Q))
  } else {
    lambda <- predict(object$internal$glmfit,
                 newdata=object$internal$glmdata,
                 type="response")
    # Note: the `newdata' argument is necessary in order to obtain predictions
    # at all quadrature points. If it is omitted then we would only get
    # predictions at the quadrature points j where glmdata$SUBSET[j]=TRUE.
  }
  
  return(lambda)
}

# ??? method for 'effects' ???



