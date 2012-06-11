#
# randomOnLines.R
#
# $Revision: 1.5 $  $Date: 2012/06/06 09:58:28 $
#
# Generate random points on specified lines
#

runifpointOnLines <- function(n, L) {
  if(!is.numeric(n) || (length(n) != 1) || (n < 0) || (n %% 1 != 0))
    stop("n should be a single nonnegative integer")
  if(!is.psp(L))
    L <- as.psp(L)
  X <- datagen.runifpointOnLines(n, L)
  out <- ppp(X$x, X$y, window=as.owin(L))
  return(out)
}

datagen.runifpointOnLines <- function(n, L) {
  stopifnot(is.psp(L))
  if(n == 0) 
    return(data.frame(x=numeric(0),
                      y=numeric(0),
                      seg=integer(0),
                      tp=numeric(0)))
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
  # convert to (x,y)
  Ldf <- as.data.frame(L)
  dx <- with(Ldf, x1-x0)
  dy <- with(Ldf, y1-y0)
  x <- with(Ldf, x0[kk] + tt * dx[kk])
  y <- with(Ldf, y0[kk] + tt * dy[kk])
  #
  return(data.frame(x=x, y=y, seg=kk, tp=tt))
}

runifpoisppOnLines <- function(lambda, L) {
  if(!(is.numeric(lambda) && (length(lambda) == 1)
       && is.finite(lambda) && (lambda >= 0)))
    stop("lambda should be a single, finite, nonnegative number")
  if(!is.psp(L))
    L <- as.psp(L)
  X <- datagen.runifpoisppOnLines(lambda, L)
  out <- ppp(X$x, X$y, window=as.owin(L))
  return(out)
}

datagen.runifpoisppOnLines <- function(lambda, L) {
  stopifnot(is.psp(L))
  mu <- lambda * sum(lengths.psp(L))
  n <- rpois(1, mu)
  df <- datagen.runifpointOnLines(n, L)
  return(df)
}

rpoisppOnLines <- function(lambda, L, lmax=NULL, ...) {
  if(!(is.numeric(lambda) || is.function(lambda) || is.im(lambda)))
    stop(paste(sQuote("lambda"),
               "must be a constant, a function or an image"))
  if(!is.psp(L))
    L <- as.psp(L)
  X <- datagen.rpoisppOnLines(lambda, L, lmax=lmax, ...)
  out <- ppp(X$x, X$y, window=as.owin(L))
}

datagen.rpoisppOnLines <- function(lambda, L, lmax=NULL, ...)  {
  stopifnot(is.psp(L))
  if(is.numeric(lambda)) 
    return(datagen.runifpoisppOnLines(lambda, L))
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
  # Lewis-Shedler (rejection) method
  Y <- datagen.runifpoisppOnLines(lmax, L)
  n <- nrow(Y)
  if(n == 0)
    return(Y)
  # evaluate lambda at each simulated point
  if(is.function(lambda)) 
    lambdaY <- lambda(Y$x, Y$y, ...)
  else
    lambdaY <- safelookup(lambda, as.ppp(Y, W=as.owin(L)))
  lambdaY[is.na(lambdaY)] <- 0
  # accept/reject
  pY <- lambdaY/lmax
  retain <- (runif(n) < pY)
  Y <- Y[retain, , drop=FALSE]
  return(Y)
}

      
  
  
  
