#
# blur.R
#
# apply Gaussian blur to an image
#
#    $Revision: 1.6 $   $Date: 2008/07/25 18:56:11 $
#
fillNA <- function(x, value=0) {
  stopifnot(is.im(x))
  v <- x$v
  v[is.na(v)] <- value
  x$v <- v
  return(x)
}
  
blur <- function(x, sigma=NULL, ..., normalise=FALSE, bleed=TRUE, varcov=NULL) {
  stopifnot(is.im(x))
  # determine smoothing kernel 
  sigma.given <- !is.null(sigma)
  varcov.given <- !is.null(varcov)
  if (sigma.given) {
    stopifnot(is.numeric(sigma))
    stopifnot(length(sigma) %in% c(1, 2))
    stopifnot(all(sigma > 0))
  }
  if (varcov.given)
    stopifnot(is.matrix(varcov) && nrow(varcov) == 2 && ncol(varcov) ==
              2)
  ngiven <- varcov.given + sigma.given
  switch(ngiven + 1,
         {
           sigma <- (1/8) * min(diff(x$xrange), diff(x$yrange))
         }, {
           if (sigma.given && length(sigma) == 2)
             varcov <- diag(sigma^2)
           if (!is.null(varcov))
             sigma <- NULL
         }, {
           stop(paste("Give only one of the arguments", sQuote("sigma"),
                      "and", sQuote("varcov")))
         })
  # replace NA's in image raster by zeroes 
  X <- fillNA(x, 0)
  # convolve with Gaussian
  Y <- second.moment.calc(X, sigma=sigma, varcov=varcov, what="smooth")
  # if no bleeding, we restrict data to the original boundary
  if(!bleed)
    Y$v[is.na(x$v)] <- NA
  # 
  if(!normalise)
    return(Y)
  # normalisation:
  # convert original image to window (0/1 image)
  Xone <- x
  isna <- is.na(x$v)
  Xone$v[isna] <- 0
  Xone$v[!isna] <- 1
  # convolve with Gaussian
  Ydenom <- second.moment.calc(Xone, sigma=sigma, ..., varcov=varcov, what="smooth")
  # normalise
  Z <- eval.im(Y/Ydenom)
  return(Z)
}
  
safelookup <- function(Z, X, factor=2, warn=FALSE) {
  # evaluates Z[X], replacing any NA's by blur(Z)[X]
  Zvals <- Z[X, drop=FALSE]
  if(any(isna <- is.na(Zvals))) {
    XX <- X[isna]
    pixdiam <- sqrt(Z$xstep^2 + Z$ystep^2)
    Zblur <- blur(Z, factor * pixdiam, bleed=TRUE, normalise=TRUE)
    Bvals <- Zblur[XX, drop=FALSE]
    if(any(is.na(Bvals)))
      stop("Internal error: pixel values were NA, even after blurring")
    Zvals[isna] <- Bvals
    warning(paste(sum(isna), "out of", X$n, "pixel values",
                  "were outside the pixel image domain",
                  "and were estimated by convolution"))
  }
  return(Zvals)
}
