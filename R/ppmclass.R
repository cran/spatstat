#
#	ppmclass.R
#
#	Class 'ppm' representing fitted point process models.
#
#
#	$Revision: 2.18 $	$Date: 2008/03/05 20:42:42 $
#
#       An object of class 'ppm' contains the following:
#
#            $method           model-fitting method (currently "mpl")
#
#            $coef             vector of fitted regular parameters
#                              as given by coef(glm(....))
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

        # ----------- Print model type -------------------
        
	cat(s$name)
        cat("\n")
        
        if(markeddata) mrk <- s$entries$marks
        if(multitype) {
            cat("Possible marks: \n")
            cat(paste(levels(mrk)))
            cat("\n")
          }

        # ----- trend --------------------------

#       cat(paste("\n", s$trend$name, ":\n", sep=""))

	if(!notrend) {
		cat("\nTrend formula: ")
		print(s$trend$formula)
        }
        
        cat(paste("\n", s$trend$label, ":", sep=""))

        tv <- s$trend$value
        if(is.list(tv)) {
          cat("\n")
          for(i in seq(tv))
            print(tv[[i]])
        } else if(is.numeric(tv) && length(tv) == 1 && is.null(names(tv)))
          # append to end of current line
          cat("\t", paste(tv), "\n")
        else {
          cat("\n")
          print(tv)
        }

        cat("\n")
        
        # ---- Interaction ----------------------------

	if(!poisson) 
          print(s$interaction, family=FALSE)

        # ---- Warnings issued ----------------------------

        probs <- s$problems
        if(!is.null(probs) && is.list(probs) && (length(probs) > 0)) 
          lapply(probs,
                 function(x) {
                   if(is.list(x) && !is.null(p <- x$print))
                     cat(paste("Problem:\n", p, "\n\n"))
                 })

        if(s$old)
          warning(paste("Model fitted by old spatstat version", s$version))
        
	return(invisible(NULL))
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


getglmfit <- function(object) {
  verifyclass(object, "ppm")
  glmfit <- object$internal$glmfit
  if(is.null(glmfit))
    return(NULL)
  if(object$method != "mpl")
    glmfit$coefficients <- object$coef
  return(glmfit)
}

    
# ??? method for 'effects' ???

valid.ppm <- function(object, na.value=TRUE) {
  verifyclass(object, "ppm")
  inte <- object$interaction
  if(is.null(inte))
    return(TRUE)
  checker <- inte$valid
  if(is.null(checker))
    return(na.value)
  Vnames <- object$internal$Vnames
  coeffs <- coef(object)
  Icoeffs <- coeffs[Vnames]
  return(checker(Icoeffs, inte))
}


logLik.ppm <- function(object, ...) {
  if(!is.poisson.ppm(object)) 
    warning(paste("log likelihood is not available for non-Poisson model;",
                  "log-pseudolikelihood returned"))
  ll <- object$maxlogpl
  attr(ll, "df") <- length(coef(object))
  class(ll) <- "logLik"
  return(ll)
}

#
# method for model.matrix

model.matrix.ppm <- function(object, ...) {
  gf <- getglmfit(object)
  if(is.null(gf)) {
    newobject <- update(object, forcefit=TRUE)
    gf <- getglmfit(newobject)
    if(is.null(gf))
      stop("internal error: unable to extract a glm fit")
  }
  mm <- model.matrix(gf, ...)
  return(mm)
}

model.images <- function(object, W=as.owin(object), ...) {
  X <- data.ppm(object)
  # make a quadscheme with a dummy point at every pixel
  Q <- pixelquad(X, W)
  # construct Berman-Turner frame
  needed <- c("trend", "interaction", "covariates", "correction", "rbord")
  bt <- do.call("bt.frame", append(list(Q), object[needed]))
  # compute model matrix
  mf <- model.frame(bt$fmla, bt$glmdata, ...)
  mm <- model.matrix(bt$fmla, mf, ...)
  # retain only the entries for dummy points (pixels)
  mm <- mm[!is.data(Q), ]
  # create template image
  Z <- as.im(attr(Q, "M"))
  # make images
  imagenames <- colnames(mm)
  result <- lapply(imagenames,
                   function(nama, Z, mm) {
                     values <- mm[, nama]
                     im(values, xcol=Z$xcol, yrow=Z$yrow,
                       lev=Z$lev, unitname=unitname(Z))
                   },
                   Z=Z, mm=mm)
  names(result) <- imagenames
  class(result) <- c("listof", class(result))
  return(result)
}
