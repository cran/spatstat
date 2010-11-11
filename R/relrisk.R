#
#    relrisk.R
#
#   Estimation of relative risk
#
#  $Revision: 1.4 $  $Date: 2010/11/08 08:46:42 $
#

relrisk <- function(X, sigma=NULL, ..., varcov=NULL, at="pixels") {
  stopifnot(is.ppp(X))
  stopifnot(is.multitype(X))
  Y <- split(X)
  ntypes <- length(Y)
  if(ntypes == 1)
    stop("Data contains only one type of points")
  marx <- marks(X)
  imarks <- as.integer(marx)
  lev <- levels(marx)
  # trap arguments
  dotargs <- list(...)
  isbwarg <- names(dotargs) %in% c("method", "nh")
  bwargs <- dotargs[isbwarg]
  dargs  <- dotargs[!isbwarg]
  # bandwidth
  if(is.null(sigma) && is.null(varcov)) {
    sigma <- do.call(bw.relrisk, append(list(X), bwargs))
  }
  # compute probabilities
  if(ntypes == 2) {
    # 1 = control, 2 = case
    # compute densities
    Deach <- do.call(density.splitppp,
                     append(list(Y, sigma=sigma, varcov=varcov, at=at),
                            dargs))
    Dall <- do.call(density.ppp,
                    append(list(X, sigma=sigma, varcov=varcov, at=at),
                           dargs))
    # compute probability of case
    switch(at,
           pixels= {
             Dcase <- Deach[[2]]
             result <- eval.im(Dcase/Dall)
           },
           points={
             result <- numeric(npoints(X))
             iscase <- (imarks == 2)
             result[iscase] <- Deach[[2]]/Dall[iscase]
             result[!iscase] <- 1 - Deach[[1]]/Dall[!iscase]
           })
  } else {
    # several types
    switch(at,
           pixels={
             Deach <- do.call(density.splitppp,
                              append(list(Y, sigma=sigma, varcov=varcov, at=at),
                                     dargs))
             Dall <- do.call(density.ppp,
                             append(list(X, sigma=sigma, varcov=varcov, at=at),
                                    dargs))
             result <- as.listof(lapply(Deach,
                                        function(d, dall) { eval.im(d/dall) },
                                        dall = Dall))
           },
           points = {
             npts <- npoints(X)
             # dummy variable matrix
             dumm <- matrix(0, npts, ntypes)
             dumm[cbind(seq(npts), imarks)] <- 1
             colnames(dumm) <- lev
             # compute probability of each type
             Z <- X %mark% dumm
             result <- do.call(smooth.ppp,
                               append(list(Z, sigma=sigma, varcov=varcov,
                                           at="points"),
                                      dargs))
           })
  }
  attr(result, "sigma") <- sigma
  attr(result, "varcov") <- varcov
  return(result)
}

bw.stoyan <- function(X, co=0.15) {
  # Stoyan's rule of thumb
  stopifnot(is.ppp(X))
  n <- npoints(X)
  W <- as.owin(X)
  a <- area.owin(W)
  stoyan <- co/sqrt(5 * n/a)
  return(stoyan)
}


bw.relrisk <- function(X, method="likelihood",nh=32) {
  stopifnot(is.ppp(X))
  stopifnot(is.multitype(X))
  Y <- split(X)
  ntypes <- length(Y)
  if(ntypes == 1)
    stop("Data contains only one type of points")
  marx <- marks(X)
  method <- pickoption("method", method,
                       c(likelihood="likelihood",
                         leastsquares="leastsquares",
                         ls="leastsquares",
                         LS="leastsquares",
                         weightedleastsquares="weightedleastsquares",
                         wls="weightedleastsquares",
                         WLS="weightedleastsquares"))
  # cross-validated bandwidth selection
  n <- npoints(X)
  W <- as.owin(X)
  a <- area.owin(W)
  d <- diameter(W)
  # Stoyan's rule of thumb applied to the least common type
  nmin <- max(1, min(table(marx)))
  stoyan <- 0.15/sqrt(nmin/a)
  # determine a range of bandwidth values 
  hmin <- max(min(nndist(unique(X))), stoyan/5)
  hmax <- max(d/4, stoyan * 5)
  h <- exp(seq(log(hmin), log(hmax), length=nh))
  # 
  # compute cross-validation criterion
  cv <- numeric(nh)
  if(method != "likelihood") {
    # dummy variables for each type
    imarks <- as.integer(marx)
    if(ntypes == 2) {
      # 1 = control, 2 = case
      indic <- (imarks == 2)
      y01   <- as.integer(indic)
    } else {
      indic <- matrix(FALSE, n, ntypes)
      indic[cbind(seq(n), imarks)] <- TRUE
      y01  <- indic * 1
    }
    X01 <- X %mark% y01
  }
  switch(method,
         likelihood={
           # for efficiency, only compute the estimate of p_j(x_i)
           # when j = m_i = mark of x_i.
           Dthis <- numeric(n)
           for(i in seq(h)) {
             Dall <- density.ppp(X, sigma=h[i], at="points", edge=FALSE)
             Deach <- density.splitppp(Y, sigma=h[i], at="points", edge=FALSE)
             split(Dthis, marx) <- Deach
             pthis <- Dthis/Dall
             cv[i] <- -mean(log(pthis))
           }
         },
         leastsquares={
           for(i in seq(h)) {
             phat <- smooth.ppp(X01, sigma=h[i], at="points", leaveoneout=TRUE)
             cv[i] <- mean((y01 - phat)^2)
           }
         },
         weightedleastsquares={
           # need initial value of h from least squares
           h0 <- bw.relrisk(X, "leastsquares", nh=ceiling(nh/2))
           phat0 <- smooth.ppp(X01, sigma=h0, at="points", leaveoneout=TRUE)
           var0 <- phat0 * (1-phat0)
           var0 <- pmax(0, 1e-6)
           for(i in seq(h)) {
             phat <- smooth.ppp(X01, sigma=h[i], at="points", leaveoneout=TRUE)
             cv[i] <- mean((y01 - phat)^2/var0)
           }
         })
  # optimize
  iopt <- which.min(cv)
  hopt <- h[iopt]
  #
  result <- hopt
  attr(result, "h") <- h
  attr(result, "cv") <- cv
  attr(result, "method") <- method
  class(result) <- "bw.relrisk"
  return(result)
}

print.bw.relrisk <- function(x, ...) {
  print(as.numeric(x), ...)
  return(invisible(NULL))
}

plot.bw.relrisk <- function(x, ...) {
  xname <- deparse(substitute(x))
  h <- attr(x, "h")
  cv <- attr(x, "cv")
  meth <- attr(x, "method")
  ylab <- paste(meth, "CV")
  do.call("plot.default",
          resolve.defaults(list(x=h, y=cv),
                           list(...),
                           list(main=xname, xlab="h", ylab=ylab),
                           list(type="l")))
  abline(v=as.numeric(x), lty=2)
  return(invisible(NULL))
}


which.max.im <- function(x) {
  stopifnot(is.list(x))
  n <- length(x)
  if(n == 0)
    return(list())
  if(!all(unlist(lapply(x, is.im))))
    stop("x should be a list of images")
  nama <- names(x)
  xmax <- x[[1]]
  wmax <- eval.im(as.integer(xmax == xmax))
  if(n > 1) {
    for(i in 2:n) {
      xi <- x[[i]]
      xmaxnew <- eval.im(pmax(xi, xmax))
      wmaxnew <- eval.im(ifelse(xi > xmax, i, wmax))
      xmax <- xmaxnew
      wmax <- wmaxnew
    }
  }
  wmax <- eval.im(factor(wmax, levels=1:n))
  if(!is.null(nama))
    levels(wmax) <- nama
  return(wmax)
}
