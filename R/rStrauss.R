rStrauss <- function(beta, gamma=1, R=0, W=owin()) {
  verifyclass(W, "owin")
  if(W$type != "rectangle")
    stop("W must be a rectangle")

  if(beta<=0)
    stop("beta should be positive")
  if(gamma<0 | gamma>1)
    stop("gamma should be between 0 and 1 (inclusive)")
  if(R<0)
    stop("interaction range R should be non-negative")

  nothing <- runif(1)

  xmin <- W$xrange[1]
  xmax <- W$xrange[2]
  ymin <- W$yrange[1]
  ymax <- W$yrange[2]
  
  nmax <- qpois(1-10^(-12),
                lambda=beta*(xmax-xmin)*(ymax-ymin))
  nout <- 0
  z <- .C("PerfectStrauss",
          beta = as.double(beta),
          gamma = as.double(gamma),
          R = as.double(R),
          xmin = as.double(xmin),
          xmax = as.double(xmax),
          ymin = as.double(ymin),
          ymax = as.double(ymax),
          nmax = as.integer(nmax),
          X = as.double(numeric(nmax)),
          Y = as.double(numeric(nmax)),
          nout=as.integer(nout),
          PACKAGE="spatstat")
  nout <- z$nout
  if(nout<0)
    stop(paste("PerfectStrauss failed! nout>",nmax))

  return(ppp(z$X[1:nout], z$Y[1:nout], window=W))
}

