#
# randomOnLines.R
#
# $Revision: 1.3 $  $Date: 2010/03/06 02:25:04 $
#
# Generate random points on specified lines
#

runifpointOnLines <- function(n, L) {
  if(!is.numeric(n) || (length(n) != 1) || (n < 0) || (n %% 1 != 0))
    stop("n should be a single nonnegative integer")
  L <- as.psp(L)
  if(n == 0)
    return(ppp(numeric(0), numeric(0), window=as.owin(L)))
  len <- lengths.psp(L)
  cumlen <- cumsum(len)
  cum0len <- c(0, cumlen)
  # generate random positions
  uu <- runif(n, min=0, max=sum(len))
  # identify segment for each point
  kk <- findInterval(uu, cum0len, rightmost.closed=TRUE, all.inside=TRUE)
  # parametric position along segment
  tt <- (uu - cum0len[kk])/len[kk]
  tt[!is.finite(tt)] <- 0
  # map to (x,y)
  Ldf <- as.data.frame(L)
  dx <- with(Ldf, x1-x0)
  dy <- with(Ldf, y1-y0)
  x <- with(Ldf, x0[kk] + tt * dx[kk])
  y <- with(Ldf, y0[kk] + tt * dy[kk])
  # return
  out <- ppp(x, y, window=as.owin(L))
  return(out)
}

runifpoisppOnLines <- function(lambda, L) {
  if(!(is.numeric(lambda) && (length(lambda) == 1)
       && is.finite(lambda) && (lambda >= 0)))
    stop("lambda should be a single, finite, nonnegative number")
  L <- as.psp(L)
  mu <- lambda * sum(lengths.psp(L))
  n <- rpois(1, mu)
  out <- runifpointOnLines(n, L)
  return(out)
}

rpoisppOnLines <- function(lambda, L, lmax=NULL, ...) {
  if(!(is.numeric(lambda) || is.function(lambda) || is.im(lambda)))
    stop(paste(sQuote("lambda"),
               "must be a constant, a function or an image"))
  if(is.numeric(lambda)) 
    return(runifpoisppOnLines(lambda, L))

  if(is.im(lambda)) {
    if(!(lambda$type %in% c("real", "integer")))
        stop("lambda must be numeric-valued or integer-valued")
    slam <- summary(lambda)
    if(any(is.infinite(slam$range)))
      stop("Infinite pixel values not permitted")
    if(slam$min < 0)
      stop("Negative pixel values not permitted")
  }
  if(is.null(lmax)) {
    # compute lmax
    if(is.function(lambda)) {
      X <- pointsOnLines(L, np=10000)
      lambdaX <- lambda(X$x, X$y, ...)
      lmax <- max(lambdaX, na.rm=TRUE)
    } else if(is.im(lambda)) 
      lmax <- slam$max
    if(!is.finite(lmax))
      stop("Infinite values of lambda obtained")
  } 

  # rejection method
  Y <- runifpoisppOnLines(lmax, L)
  n <- Y$n
  if(n == 0)
    return(Y)
  # evaluate lambda at each simulated point
  if(is.function(lambda)) 
    lambdaY <- lambda(Y$x, Y$y, ...)
  else
    lambdaY <- safelookup(lambda, Y)
  lambdaY[is.na(lambdaY)] <- 0
  # accept/reject
  pY <- lambdaY/lmax
  retain <- (runif(n) < pY)
  Y <- Y[retain]
  return(Y)
}

      
  
  
  
