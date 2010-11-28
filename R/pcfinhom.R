#
#   pcfinhom.R
#
#   $Revision: 1.2 $   $Date: 2010/11/21 04:22:52 $
#
#   inhomogeneous pair correlation function of point pattern 
#
#

pcfinhom <- function(X, lambda=NULL, ..., r=NULL,
                     kernel="epanechnikov", bw=NULL, stoyan=0.15,
                     correction=c("translate", "Ripley"),
                     sigma=NULL, varcov=NULL)
{
  verifyclass(X, "ppp")
  r.override <- !is.null(r)

  win <- X$window
  area <- area.owin(win)
  npts <- npoints(X)
  
  correction.given <- !missing(correction)
  correction <- pickoption("correction", correction,
                           c(isotropic="isotropic",
                             Ripley="isotropic",
                             translate="translate",
                             best="best"),
                           multi=TRUE)

  correction <- implemented.for.K(correction, win$type, correction.given)
  
  if(is.null(bw) && kernel=="epanechnikov") {
    # Stoyan & Stoyan 1995, eq (15.16), page 285
    h <- stoyan /sqrt(npts/area)
    # conversion to standard deviation
    bw <- h/sqrt(5)
  }

  ########## intensity values #########################

  if(is.null(lambda)) {
    # No intensity data provided
    # Estimate density by leave-one-out kernel smoothing
    lambda <- density(X, ..., sigma=sigma, varcov=varcov,
                      at="points", leaveoneout=TRUE)
  } else {
    # lambda values provided
    if(is.vector(lambda)) 
      check.nvector(lambda, npts)
    else if(is.im(lambda)) 
      lambda <- safelookup(lambda, X)
    else if(is.function(lambda)) 
      lambda <- lambda(X$x, X$y)
    else stop(paste(sQuote("lambda"),
                    "should be a vector, a pixel image, or a function"))
  }
  # compute reciprocal
  reciplambda <- 1/lambda
  
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

  ########## smoothing parameters for pcf ############################  
  # arguments for 'density'

  from <- 0
  to <- max(r)
  nr <- length(r)
  
  denargs <- resolve.defaults(list(kernel=kernel, bw=bw),
                              list(...),
                              list(n=nr, from=from, to=to))
  
  #################################################
  
  # compute pairwise distances
  
  close <- closepairs(X, max(r))
  dIJ <- close$d
  I <- close$i
  J <- close$j
  XI <- ppp(close$xi, close$yi, window=win, check=FALSE)
  wIJ <- reciplambda[I] * reciplambda[J]

  # initialise fv object
  
  df <- data.frame(r=r, theo=rep(1,length(r)))
  out <- fv(df, "r",
            substitute(ginhom(r), NULL), "theo", ,
            alim,
            c("r","%s[Pois](r)"),
            c("distance argument r", "theoretical Poisson %s"),
            fname="ginhom")

  ###### compute #######

  if(any(correction=="translate")) {
    # translation correction
    XJ <- ppp(close$xj, close$yj, window=win, check=FALSE)
    edgewt <- edge.Trans(XI, XJ, paired=TRUE)
    gT <- sewpcf(dIJ, edgewt * wIJ, denargs, area)$g
    out <- bind.fv(out,
                   data.frame(trans=gT),
                   "%s[Trans](r)",
                   "translation-corrected estimate of %s",
                   "trans")
  }
  if(any(correction=="isotropic")) {
    # Ripley isotropic correction
    edgewt <- edge.Ripley(XI, matrix(dIJ, ncol=1))
    gR <- sewpcf(dIJ, edgewt * wIJ, denargs, area)$g
    out <- bind.fv(out,
                   data.frame(iso=gR),
                   "%s[Ripley](r)",
                   "isotropic-corrected estimate of %s",
                   "iso")
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
  unitname(out) <- unitname(X)
  return(out)
}

