#
#    relrisk.R
#
#   Estimation of relative risk
#
#  $Revision$  $Date$
#

relrisk <- function(X, sigma=NULL, ..., varcov=NULL, at="pixels") {
  stopifnot(is.ppp(X))
  stopifnot(is.multitype(X))
  Y <- split(X)
  ntypes <- length(Y)
  if(ntypes == 1)
    stop("Data contains only one type of points")
  # trap arguments
  dotargs <- list(...)
  isbwarg <- names(dotargs) %in% c("method", "nh")
  bwargs <- dotargs[isbwarg]
  dargs  <- dotargs[!isbwarg]
  # bandwidth
  if(is.null(sigma) && is.null(varcov)) {
    sigma <- do.call(bw.relrisk, append(list(X), bwargs))
  }
  # compute densities
  Deach <- do.call(density.splitppp,
                   append(list(Y, sigma=sigma, varcov=varcov, at=at),
                          dargs))
  Dall <- do.call(density.ppp,
                  append(list(X, sigma=sigma, varcov=varcov, at=at),
                          dargs))
  # compute probabilities
  if(ntypes == 2) {
    # 1 = control, 2 = case
    # compute probability of case
    switch(at,
           pixels= {
             Dcase <- Deach[[2]]
             result <- eval.im(Dcase/Dall)
           },
           points={
             result <- numeric(npoints(X))
             iscase <- (as.integer(marks(X)) == 2)
             result[iscase] <- Deach[[2]]/Dall[iscase]
             result[!iscase] <- 1 - Deach[[1]]/Dall[!iscase]
           })
  } else {
    # several types;
    # compute probability of each type
    switch(at,
         pixels={
           result <- as.listof(lapply(Deach,
                                     function(d, dall) { eval.im(d/dall) },
                                     dall=Dall))
         },
         points={
           Dallsplit <- split(Dall, marks(X))
           Yprob <- list()
           for(i in seq(length(Y))) 
             Yprob[[i]] <- Y[[i]] %mark% (Deach[[i]]/Dallsplit[[i]])
           Xnew <- X
           split(Xnew, un=FALSE) <- Yprob
           result <- marks(Xnew)
         })
  }
  attr(result, "sigma") <- sigma
  attr(result, "varcov") <- varcov
  return(result)
}


bw.relrisk <- function(X, method="likelihood",nh=32) {
  stopifnot(is.ppp(X))
  stopifnot(is.multitype(X))
  Y <- split(X)
  ntypes <- length(Y)
  if(ntypes == 1)
    stop("Data contains only one type of points")
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
  # Stoyan's rule of thumb
  stoyan <- 0.15/sqrt(n/a)
  # determine a range of bandwidth values 
  hmin <- max(min(nndist(unique(X))), stoyan/5)
  hmax <- max(d/4, stoyan * 5)
  h <- exp(seq(log(hmin), log(hmax), length=nh))
  # dummy variables for each type
  if(ntypes == 2) {
    # 1 = control, 2 = case
    indic <- (as.integer(marks(X)) == 2)
    y01   <- as.integer(indic)
  } else {
    indic <- outer(as.integer(marks(X)), 1:ntypes, "==")
    y01  <- indic * 1
  }
  X01 <- X %mark% y01
  # compute cross-validation criterion
  cv <- numeric(nh)
  switch(method,
         likelihood={
           for(i in seq(h)) {
             phat <- smooth.ppp(X01, sigma=h[i], at="points", leaveoneout=TRUE)
             cv[i] <- -mean(log(phat[indic])) 
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

