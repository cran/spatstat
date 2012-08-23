#
#  lohboot.R
#
#  $Revision: 1.4 $   $Date: 2012/08/22 01:42:47 $
#
#  Loh's bootstrap CI's for local pcf, local K etc
#

lohboot <-
  function(X,
           fun=c("pcf", "Kest", "pcfinhom", "Kinhom"),
           ..., nsim=200, confidence=0.95, type=7) {
  stopifnot(is.ppp(X))
  fun <- match.arg(fun)
  # validate confidence level
  stopifnot(confidence > 0.5 && confidence < 1)
  alpha <- 1 - confidence
  probs <- c(alpha/2, 1-alpha/2)
  rank <- nsim * probs[2]
  if(abs(rank - round(rank)) > 0.001)
    warning(paste("confidence level", confidence,
                  "corresponds to a non-integer rank", paren(rank),
                  "so quantiles will be interpolated"))
  #
  n <- npoints(X)
  localfun <- switch(fun,
                     pcf=localpcf,
                     Kest=localK,
                     pcfinhom=localpcfinhom,
                     Kinhom=localKinhom)
  f <- localfun(X, ...)
  # parse edge correction info
  correction <- attr(f, "correction")
  switch(correction,
         none      = { ctag <- "un";    cadj <- "uncorrected" },
         border    = { ctag <- "bord";  cadj <- "border-corrected" },
         translate = { ctag <- "trans"; cadj <- "translation-corrected" },
         isotropic = { ctag <- "iso";   cadj <- "Ripley isotropic corrected" })
  # first n columns are the local pcfs for the n points of X
  y <- as.matrix(as.data.frame(f))[, 1:n]
  # average them
  ymean <- rowMeans(y, na.rm=TRUE)
  # resample
  ystar <- matrix(, nrow=nrow(y), ncol=nsim)
  for(i in 1:nsim) {
    # resample n points with replacement
    ind <- sample(n, replace=TRUE)
    # average their local pcfs
    ystar[,i] <- rowMeans(y[,ind], na.rm=TRUE)
  }
  # compute quantiles
  hilo <- apply(ystar, 1, quantile,
                probs=probs, na.rm=TRUE, type=type)
  # create fv object
  df <- data.frame(r=f$r,
                   theo=f$theo,
                   ymean,
                   lo=hilo[1,],
                   hi=hilo[2,])
  colnames(df)[3] <- ctag
  CIlevel <- paste(100 * confidence, "%% confidence", sep="")
  desc <- c("distance argument r",
            "theoretical Poisson %s",
            paste(cadj, "estimate of %s"),
            paste("lower", CIlevel, "limit for %s"),
            paste("upper", CIlevel, "limit for %s"))
  clabl <- paste("hat(%s)[", ctag, "](r)", sep="")
  labl <- c("r", "%s[pois](r)", clabl, "%s[loCI](r)", "%s[hiCI](r)")
  switch(fun,
         pcf={ fname <- "g" ; ylab <- quote(g(r)) },
         Kest={ fname <- "K" ; ylab <- quote(K(r)) },
         pcfinhom={ fname <- "g[inhom]" ; ylab <- quote(g[inhom](r)) },
         Kinhom={ fname <- "K[inhom]" ; ylab <- quote(K[inhom](r)) })
  g <- fv(df, "r", ylab, ctag, , c(0, max(f$r)), labl, desc, fname=fname)
  # default is to display them all
  formula(g) <- . ~ r
  fvnames(g, ".") <- c(ctag, "hi", "lo", "theo")
  unitname(g) <- unitname(X)
  g
}


    
  
  
  
  
