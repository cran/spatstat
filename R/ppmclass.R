#
#	ppmclass.R
#
#	Class 'ppm' representing fitted point process models.
#
#
#	$Revision: 2.44 $	$Date: 2010/07/16 05:19:06 $
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



        tv <- s$trend$value

        if(length(tv) == 0)
          cat("\n[No trend coefficients]\n")
        else {
          cat(paste("\n", s$trend$label, ":", sep=""))
          if(is.list(tv)) {
            cat("\n")
            for(i in seq(tv))
              print(tv[[i]])
          } else if(is.numeric(tv) && length(tv) == 1 && is.null(names(tv))) {
            # single number: append to end of current line
            cat("\t", paste(tv), "\n")
          } else {
            # some other format 
            cat("\n")
            print(tv)
          }
        }

        cat("\n")
        
        # ---- Interaction ----------------------------

	if(!poisson) 
          print(s$interaction, family=FALSE)

        # ---- Warnings issued in mpl.prepare  ---------------------
        
        probs <- s$problems
        if(!is.null(probs) && is.list(probs) && (length(probs) > 0)) 
          lapply(probs,
                 function(x) {
                   if(is.list(x) && !is.null(p <- x$print))
                     cat(paste("Problem:\n", p, "\n\n"))
                 })

        if(s$old)
          warning(paste("Model fitted by old spatstat version", s$version))
        
        # ---- Algorithm status ----------------------------

        fitter <- s$fitter
        converged <- s$converged
        if(!is.null(fitter) && fitter %in% c("glm", "gam") && !converged)
          cat(paste("*** Fitting algorithm for", sQuote(fitter),
                    "did not converge ***\n"))

	return(invisible(NULL))
}

quad.ppm <- function(object, drop=FALSE) {
  verifyclass(object, "ppm")
  Q <- object$Q
  if(!drop || is.null(Q))
    return(Q)
  ok <- object$internal$glmdata$.mpl.SUBSET
  if(is.null(ok))
    return(Q)
  return(Q[ok])
}

data.ppm <- function(object) { 
  verifyclass(object, "ppm")
  object$Q$data
}

dummy.ppm <- function(object, drop=FALSE) { 
  return(quad.ppm(object, drop=drop)$dummy)
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

getglmdata <- function(object, drop=FALSE) {
  verifyclass(object, "ppm")
  gd <- object$internal$glmdata
  if(!drop) return(gd)
  return(gd[gd$.mpl.SUBSET,])
}

getglmsubset <- function(object) {
  gd <- object$internal$glmdata
  return(gd$.mpl.SUBSET)
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


logLik.ppm <- function(object, ..., warn=TRUE) {
  if(!is.poisson.ppm(object) && warn) 
    warning(paste("log likelihood is not available for non-Poisson model;",
                  "log-pseudolikelihood returned"))
  method <- object$method
  switch(method,
         mpl={
           ll <- object$maxlogpl
         },
         ho={
           # evaluate the log pseudolikelihood
           Q <- quad.ppm(object, drop=TRUE)
           Z <- is.data(Q)
           w <- w.quad(Q)
           cif <- fitted(object, type="cif", drop=TRUE)
           cifdata <- cif[Z]
           ll <- sum(log(cifdata[cifdata > 0])) - sum(w * cif)
         },
         stop(paste("Internal error: unrecognised ppm method:",
                    dQuote(method)))
         )
  attr(ll, "df") <- length(coef(object))
  class(ll) <- "logLik"
  return(ll)
}

# more methods

formula.ppm <- function(x, ...) {
  f <- x$trend
  if(is.null(f)) f <- ~1
  return(f)
}

terms.ppm <- function(x, ...) {
  terms(formula(x), ...)
}

extractAIC.ppm <- function (fit, scale = 0, k = 2, ...)
{
    edf <- length(coef(fit))
    aic <- AIC(fit)
    c(edf, aic + (k - 2) * edf)
}


#
# method for model.frame

model.frame.ppm <- function(formula, ...) {
  object <- formula
  gf <- getglmfit(object)
  if(is.null(gf)) {
    warning("Model re-fitted with forcefit=TRUE")
    object <- update(object, forcefit=TRUE)
    gf <- getglmfit(object)
  }
  gd <- getglmdata(object)
  model.frame(gf, data=gd, ...)
}

#
# method for model.matrix

model.matrix.ppm <- function(object, data=model.frame(object),
                             ..., keepNA=TRUE) {
  gf <- getglmfit(object)
  if(is.null(gf)) {
    warning("Model re-fitted with forcefit=TRUE")
    object <- update(object, forcefit=TRUE)
    gf <- getglmfit(object)
    if(is.null(gf))
      stop("internal error: unable to extract a glm fit")
  }
  if(!keepNA) {
    # extract model matrix of glm fit object
    # restricting to its 'subset' 
    mm <- model.matrix(gf, data, ...)
    return(mm)
  }
  # extract model matrix for all cases
  mm <- model.matrix(gf, data, ..., subset=NULL)
  cn <- colnames(mm)
  gd <- getglmdata(object, drop=FALSE)
  if(nrow(mm) != nrow(gd)) {
    # can occur if covariates include NA's or interaction is -Inf
    insubset <- getglmsubset(object)
    isna <- is.na(insubset) | !insubset
    if(sum(isna) + nrow(mm) == nrow(gd)) {
      # insert rows of NA's
      mmplus <- matrix( , nrow(gd), ncol(mm))
      mmplus[isna, ] <- NA
      mmplus[!isna, ] <- mm
      mm <- mmplus
    } else 
    stop("internal error: model matrix does not match glm data frame")
  }
  colnames(mm) <- cn
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
  mm <- mm[!is.data(Q), , drop=FALSE]
  # create template image
  Z <- as.im(attr(Q, "M"))
  # make images
  imagenames <- colnames(mm)
  result <- lapply(imagenames,
                   function(nama, Z, mm) {
                     values <- mm[, nama]
                     im(values, xcol=Z$xcol, yrow=Z$yrow,
                        unitname=unitname(Z))
                   },
                   Z=Z, mm=mm)
  result <- as.listof(result)
  names(result) <- imagenames
  return(result)
}
