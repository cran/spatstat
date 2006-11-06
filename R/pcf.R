#
#   pcf.R
#
#   $Revision: 1.23 $   $Date: 2006/10/31 07:45:35 $
#
#
#   calculate pair correlation function
#   from point pattern (pcf.ppp)
#   or from estimate of K or Kcross (pcf.fv)
#   or from fasp object
#
#
pcf <- function(X, ...) {
  UseMethod("pcf")
}

pcf.ppp <- function(X, ..., r=NULL,
                    kernel="epanechnikov", bw=NULL, stoyan=0.15,
                    correction=c("translate", "ripley"))
{
  verifyclass(X, "ppp")
  r.override <- !is.null(r)

  win <- X$window
  area <- area.owin(win)
  lambda <- X$n/area

  stopifnot(is.character(correction))
  if(!all(ok <- correction %in% c("translate", "ripley")))
    stop(paste("unrecognised correction",
               if(sum(!ok) > 1) "s", ": ",
               paste(correction[!ok], collapse=", "),
               sep=""))
  
  if(is.null(bw) && kernel=="epanechnikov") {
    # Stoyan & Stoyan 1995, eq (15.16), page 285
    h <- stoyan /sqrt(lambda)
    # conversion to standard deviation
    bw <- h/sqrt(5)
  }

  ########## r values ############################
  # handle arguments r and breaks 

  rmaxdefault <- rmax.rule("K", win, lambda)        
  breaks <- handle.r.b.args(r, NULL, win, rmaxdefault=rmaxdefault)
  if(!(breaks$even))
    stop("r values must be evenly spaced")
  # extract r values
  r <- breaks$r
  rmax <- breaks$max
  # recommended range of r values for plotting
  alim <- c(0, min(rmax, rmaxdefault))

  # arguments for 'density'
  from <- 0
  to <- max(r)
  nr <- length(r)
  #################################################
  
  
  # compute pairwise distances
  
  close <- closepairs(X, max(r))
  dIJ <- close$d
  XI <- ppp(close$xi, close$yi, window=win)

  # how to process the distances
  
  doit <- function(w, d, out, symb, desc, key, otherargs, lambda, area) {
    kden <- do.call("densityhack",
                    resolve.defaults(list(x=d, weights=w), otherargs))
                                     
    r <- kden$x
    y <- kden$y
    g <- y/(2 * pi * r * (lambda^2) * area)
    if(is.null(out)) {
      df <- data.frame(r=r, theo=rep(1,length(r)), g)
      colnames(df)[3] <- key
      out <- fv(df, "r", substitute(g(r), NULL), key, , alim, c("r","gPois(r)", symb),
                c("distance argument r", "theoretical Poisson g(r)", desc))
    } else {
      df <- data.frame(g)
      colnames(df) <- key
      out <- bind.fv(out, df, symb, desc, key)
    }
    return(out)
  }

  otherargs <- resolve.defaults(list(kernel=kernel, bw=bw),
                                list(...),
                                list(n=nr, from=from, to=to))

  ###### compute #######

  out <- NULL

  if(any(correction=="translate")) {
    # translation correction
    XJ <- ppp(close$xj, close$yj, window=win)
    edgewt <- edge.Trans(XI, XJ, paired=TRUE)
    out <- doit(edgewt, dIJ, out, "gTrans(r)",
                "translation-corrected estimate of g", "trans", otherargs,
                lambda, area)
  }
  if(any(correction=="ripley")) {
    # Ripley isotropic correction
    edgewt <- edge.Ripley(XI, matrix(dIJ, ncol=1))
    out <- doit(edgewt, dIJ, out, "gRipley(r)",
                "Ripley-corrected estimate of g", "ripl", otherargs,
                lambda, area)
  }
  
  # sanity check
  if(is.null(out)) {
    warning("Nothing computed - no edge corrections chosen")
    return(NULL)
  }
  
  # which corrections have been computed?
  nama2 <- names(out)
  corrxns <- rev(nama2[nama2 != "r"])

  # default is to display them all
  attr(out, "fmla") <- deparse(as.formula(paste(
                       "cbind(",
                        paste(corrxns, collapse=","),
                        ") ~ r")))
  units(out) <- units(X)
  return(out)
}


"pcf.fasp" <- function(X, ..., method="c") {
  verifyclass(X, "fasp")
  Y <- X
  Y$title <- paste("Array of pair correlation functions",
                   if(!is.null(X$dataname)) "for",
                   X$dataname)
  # go to work on each function
  for(i in seq(X$fns)) {
    Xi <- X$fns[[i]]
    PCFi <- pcf.fv(Xi, ..., method=method)
    Y$fns[[i]] <- as.fv(PCFi)
    if(is.fv(PCFi))
      Y$default.formula[[i]] <- attr(PCFi, "fmla")
  }
  return(Y)
}


"pcf.fv" <-
function(X, ..., method="c") {
  verifyclass(X, "fv")
  
  callmatched <- function(fun, argue) {
    formalnames <- names(formals(fun))
    formalnames <- formalnames[formalnames != "..."]
    do.call("fun", argue[names(argue) %in% formalnames])
  }

  # extract r and the recommended estimate of K
  r <- X[[attr(X, "argu")]]
  K <- X[[attr(X, "valu")]]
  alim <- attr(X, "alim")

  # remove NA's
  ok <- !is.na(K)
  K <- K[ok]
  r <- r[ok]
  switch(method,
         a = {
           ss <- callmatched(smooth.spline,
                             list(x=r, y=K, ...))
           dK <- predict(ss, r, deriv=1)$y
           g <- dK/(2 * pi * r)
         },
         b = {
           y <- K/(2 * pi * r)
           y[!is.finite(y)] <- 0
           ss <- callmatched(smooth.spline,
                             list(x=r, y=y, ...))
           dy <- predict(ss, r, deriv=1)$y
           g <- dy + y/r
         },
         c = {
           z <- K/(pi * r^2)
           z[!is.finite(z)] <- 1
           ss <- callmatched(smooth.spline,
                             list(x=r, y=z, ...))
           dz <- predict(ss, r, deriv=1)$y
           g <- (r/2) * dz + z
         },
         d = {
           z <- sqrt(K)
           z[!is.finite(z)] <- 0
           ss <- callmatched(smooth.spline,
                             list(x=r, y=z, ...))
           dz <- predict(ss, r, deriv=1)$y
           g <- z * dz/(pi * r)
         },
         stop(paste("unrecognised method", sQuote(method)))
         )

  # pack result into "fv" data frame
  Z <- fv(data.frame(r=r, pcf=g, theo=rep(1, length(r))),
          "r", substitute(pcf(r), NULL), "pcf", cbind(pcf, theo) ~ r, alim,
          c("r", "pcf(r)", "1"),
          c("distance argument r",
            "estimate of pair correlation function pcf(r)",
            "theoretical Poisson value, pcf(r) = 1"))
  units(Z) <- units(X)
  return(Z)
}
