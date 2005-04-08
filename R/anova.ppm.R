#
#   anova.ppm.R
#
#  $Revision: 1.5 $   $Date: 2005/04/12 22:20:13 $
#

anova.ppm <- function(object, ..., test=NULL, override=FALSE) {

  # list of models
  objex <- append(list(object), list(...))
  if(!all(unlist(lapply(objex, is.ppm))))
    stop("Arguments must all be \"ppm\" objects")

  # non-Poisson models?
  pois <- all(unlist(lapply(objex, is.poisson.ppm)))
  if(!pois) {
    whinge <- paste("Some of the fitted models are not Poisson processes:",
                    "p-values are not supported by any theory")
    if(override)
      warning(whinge)
    else
      stop(whinge)
  }

  # all models fitted by MPL?
  mplfit <- unlist(lapply(objex, function(x) { x$method=="mpl" }))
  if(!all(mplfit)) 
    stop(paste("Not all models fitted by maximum pseudolikelihood;",
               "comparison not possible"))

  # Extract glmfit objects 
  fitz <- lapply(objex, getglmfit)
  
  # Are any of them NULL (corresponding to uniform Poisson)
  not.glm <- unlist(lapply(fitz, is.null))
  if(any(not.glm)) {
    # force them to be fitted using glm
    objex[not.glm] <- lapply(objex[not.glm], update.ppm, forcefit=TRUE)
    fitz[not.glm] <- lapply(objex[not.glm], getglmfit)
  }
  
  # check whether all models fitted using GLM or all using GAM
  isgam <- unlist(lapply(fitz[!is.null(fitz)],
                         function(x) { inherits(x, "gam") }))
  if(any(isgam) && !all(isgam))
    warning("Some, but not all, models were fitted with use.gam=TRUE;",
            "anova may be incorrect.",
            "It is recommended to refit all models with use.gam=TRUE.")
  anovafun <- if(any(isgam)) "anova.gam" else "anova.glm"
  
  result <- do.call(anovafun, append(fitz, list(test=test, dispersion=1)))
  return(result)
}
