#
#
anova.ppm <- function(object, ..., test=NULL, override=FALSE) {
  objex <- append(list(object), list(...))
  if(!all(unlist(lapply(objex, is.ppm))))
    stop("Arguments must all be \"ppm\" objects")
  pois <- all(unlist(lapply(objex, is.poisson.ppm)))
  if(!pois)
    warning("Some of the fitted models are not Poisson processes:",
            "p-values are not supported by any theory")
  if(pois || override) {
    fitz <- lapply(objex, function(x) { x$internal$glmfit })
    do.call("anova.glm", append(fitz, list(test=test, dispersion=1)))
  } else
    NULL
}
