#
#  Perfect Simulation 
#
#  $Revision: 1.10 $ $Date: 2011/12/17 09:29:25 $
#
#  rStrauss
#  rHardcore
#  rDiggleGratton

rStrauss <- function(beta, gamma=1, R=0, W=owin()) {
  if(!missing(W)) {
    verifyclass(W, "owin")
    if(W$type != "rectangle")
      stop("W must be a rectangle")
  }

  check.1.real(beta)
  check.1.real(gamma)
  check.1.real(R)

  check.finite(beta)
  check.finite(gamma)
  check.finite(R)
  
  stopifnot(beta > 0)
  stopifnot(gamma >= 0)
  stopifnot(gamma <= 1)
  stopifnot(R >= 0)

  nothing <- runif(1)

  xrange <- W$xrange
  yrange <- W$yrange
  storage.mode(beta) <- storage.mode(gamma) <- storage.mode(R) <- "double"
  storage.mode(xrange) <- storage.mode(yrange) <- "double"
  
  z <- .Call("PerfectStrauss",
             beta,
             gamma,
             R,
             xrange,
             yrange,
             PACKAGE="spatstat")

  X <- z[[1]]
  Y <- z[[2]]
  nout <- z[[3]]

  if(nout<0)
    stop("internal error: copying failed in PerfectStrauss")

  return(ppp(X[1:nout], Y[1:nout], window=W, check=FALSE))
}

#  Perfect Simulation of Hardcore process

rHardcore <- function(beta, R=0, W=owin()) {
  if(!missing(W)) {
    verifyclass(W, "owin")
    if(W$type != "rectangle")
      stop("W must be a rectangle")
  }

  check.1.real(beta)
  check.1.real(R)

  check.finite(beta)
  check.finite(R)

  stopifnot(beta > 0)
  stopifnot(R    >= 0)

  nothing <- runif(1)

  xrange <- W$xrange
  yrange <- W$yrange
  storage.mode(beta) <- storage.mode(R) <- "double"
  storage.mode(xrange) <- storage.mode(yrange) <- "double"
  
  z <- .Call("PerfectHardcore",
             beta,
             R,
             xrange,
             yrange,
             PACKAGE="spatstat")

  X <- z[[1]]
  Y <- z[[2]]
  nout <- z[[3]]

  if(nout<0)
    stop("internal error: copying failed in PerfectHardcore")

  return(ppp(X[1:nout], Y[1:nout], window=W, check=FALSE))
}

#
#  Perfect Simulation of Diggle-Gratton process
#

rDiggleGratton <- function(beta, delta, rho, kappa=1, W=owin()) {
  if(!missing(W)) {
    verifyclass(W, "owin")
    if(W$type != "rectangle")
      stop("W must be a rectangle")
  }

  check.1.real(beta)
  check.1.real(delta)
  check.1.real(rho)
  check.1.real(kappa)

  check.finite(beta)
  check.finite(delta)
  check.finite(rho)
  check.finite(kappa)

  stopifnot(beta > 0)
  stopifnot(delta >= 0)
  stopifnot(rho   >= 0)
  stopifnot(delta <= rho)
  stopifnot(kappa >= 0)

  nothing <- runif(1)

  xrange <- W$xrange
  yrange <- W$yrange
  storage.mode(beta) <- "double"
  storage.mode(delta) <- storage.mode(rho) <- storage.mode(kappa) <- "double"
  storage.mode(xrange) <- storage.mode(yrange) <- "double"
  
  z <- .Call("PerfectDiggleGratton",
             beta,
             delta,
             rho,
             kappa,
             xrange,
             yrange,
             PACKAGE="spatstat")

  X <- z[[1]]
  Y <- z[[2]]
  nout <- z[[3]]

  if(nout<0)
    stop("internal error: copying failed in PerfectDiggleGratton")

  return(ppp(X[1:nout], Y[1:nout], window=W, check=FALSE))
}


#
#  Perfect Simulation of Diggle-Gates-Stibbard process
#

rDGS <- function(beta, rho, W=owin()) {
  if(!missing(W)) {
    verifyclass(W, "owin")
    if(W$type != "rectangle")
      stop("W must be a rectangle")
  }

  check.1.real(beta)
  check.1.real(rho)

  check.finite(beta)
  check.finite(rho)

  stopifnot(beta > 0)
  stopifnot(rho  >= 0)

  nothing <- runif(1)

  xrange <- W$xrange
  yrange <- W$yrange
  storage.mode(beta) <- "double"
  storage.mode(rho) <- "double"
  storage.mode(xrange) <- storage.mode(yrange) <- "double"
  
  z <- .Call("PerfectDGS",
             beta,
             rho,
             xrange,
             yrange,
             PACKAGE="spatstat")

  X <- z[[1]]
  Y <- z[[2]]
  nout <- z[[3]]

  if(nout<0)
    stop("internal error: copying failed in PerfectDGS")

  return(ppp(X[1:nout], Y[1:nout], window=W, check=FALSE))
}


