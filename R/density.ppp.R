#
#  density.ppp.R
#
#  Method for 'density' for point patterns
#
#  $Revision: 1.1 $    $Date: 2007/08/17 14:29:20 $
#

ksmooth.ppp <- function(x, sigma, ..., edge=TRUE) {
  .Deprecated("density.ppp", package="spatstat")
  density.ppp(x, sigma, ..., edge=edge)
}

density.ppp <- function(x, sigma, ..., weights=NULL, edge=TRUE, varcov=NULL) {
  verifyclass(x, "ppp")
  sigma.given <- !missing(sigma) && !is.null(sigma)
  varcov.given <- !is.null(varcov)
  if(sigma.given) {
    stopifnot(is.numeric(sigma))
    stopifnot(length(sigma) %in% c(1,2))
    stopifnot(all(sigma > 0))
  }
  if(varcov.given)
    stopifnot(is.matrix(varcov) && nrow(varcov) == 2 && ncol(varcov)==2 )    
  ngiven <- varcov.given + sigma.given
  switch(ngiven+1,
         {
           # default
           w <- x$window
           sigma <- (1/8) * min(diff(w$xrange), diff(w$yrange))
         },
         {
           if(sigma.given && length(sigma) == 2) 
             varcov <- diag(sigma^2)
           if(!is.null(varcov))
             sigma <- NULL
         },
         {
           stop(paste("Give only one of the arguments",
                      sQuote("sigma"), "and", sQuote("varcov")))
         })
      
  smo <- second.moment.calc(x, sigma, what="smooth", ..., weights=weights, varcov=varcov)
  smo$v <- smo$v/(smo$xstep * smo$ystep)
  raw <- smo
  if(edge) {
    edg <- second.moment.calc(x, sigma, what="edge", ..., weights=weights, varcov=varcov)
    smo <- eval.im(smo/edg)
  }
  result <- smo[x$window, drop=FALSE]

  # internal use only
  spill <- list(...)$spill
  if(!is.null(spill)) {
    edg <- if(edge) im(edg, xcol=raw$xcol, yrow=raw$yrow) else NULL
    return(list(sigma=sigma, varcov=varcov, raw = raw, edg=edg))
  }

  # normal return
  return(result)
}

smooth.ppp <- function(X, ..., weights=rep(1,X$n)) {
  verifyclass(X, "ppp")
  if(is.marked(X)) {
    if(is.factor(marks(X)))
      warning("Factor values were converted to integers")
    values <- as.numeric(marks(X))
  }
  else
    values <- rep(1, X$n)
  numerator <-   density(X, ..., weights= values * weights)
  denominator <- density(X, ..., weights= weights)
  ratio <- eval.im(numerator/denominator)
  return(ratio)
}

  
