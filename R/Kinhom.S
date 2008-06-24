#
#	Kinhom.S	Estimation of K function for inhomogeneous patterns
#
#	$Revision: 1.32 $	$Date: 2008/06/18 21:07:01 $
#
#	Kinhom()	compute estimate of K_inhom
#
#
#       Reference:
#            Non- and semiparametric estimation of interaction
#	     in inhomogeneous point patterns
#            A.Baddeley, J.Moller, R.Waagepetersen
#            Statistica Neerlandica 54 (2000) 329--350.
#
# -------- functions ----------------------------------------
#	Kinhom()	compute estimate of K
#                       using various edge corrections
#
#       Kwtsum()         internal routine for border correction
#
# -------- standard arguments ------------------------------	
#	X		point pattern (of class 'ppp')
#
#	r		distance values at which to compute K	
#
#       lambda          vector of intensity values for points of X
#
# -------- standard output ------------------------------
#      A data frame (class "fv") with columns named
#
#	r:		same as input
#
#	trans:		K function estimated by translation correction
#
#	iso:		K function estimated by Ripley isotropic correction
#
#	theo:		K function for Poisson ( = pi * r ^2 )
#
#	border:		K function estimated by border method
#			(denominator = sum of weights of points)
#
#       bord.modif:	K function estimated by border method
#			(denominator = area of eroded window)
#
# ------------------------------------------------------------------------

"Linhom" <- function(...) {
  K <- Kinhom(...)
  L <- eval.fv(sqrt(K/pi))
  # relabel the fv object
  attr(L, "ylab") <- substitute(Linhom(r), NULL)
  nama <- names(L)
  desc <- attr(L, "desc")
  labl <- attr(L, "labl")
  TAGS <- c("theo", "un", "bord.modif", "border", "trans", "iso")
  DESC <- c("theoretical Poisson L(r)",
             "uncorrected estimate of Linhom(r)",
             "modified border-corrected estimate of Linhom(r)",
             "border-corrected estimate of Linhom(r)",
             "translation-corrected estimate of Linhom(r)",
             "Ripley isotropic correction estimate of Linhom(r)")
  LABL <- c("Lpois(r)",
             "Lun(r)",
             "Lbord*(r)",
             "Lbord(r)",
             "Ltrans(r)",
             "Liso(r)")

  for(i in seq(length(TAGS)))
    if(!is.na(m <- match(TAGS[i], nama))) {
      desc[m] <- DESC[i]
      labl[m] <- LABL[i]
    }
  attr(L, "desc") <- desc
  attr(L, "labl") <- labl
  return(L)  
}

"Kinhom"<-
  function (X, lambda, ..., r = NULL, breaks = NULL, 
         correction=c("border", "bord.modif", "isotropic", "translate"),
            nlarge = 1000, 
            lambda2=NULL,
            sigma=NULL, varcov=NULL)
{
    verifyclass(X, "ppp")
    nlarge.given <- !missing(nlarge)
    rfixed <- !missing(r) || !missing(breaks)
    
    W <- X$window
    npoints <- X$n
    area <- area.owin(W)

    rmaxdefault <- rmax.rule("K", W, npoints/area)
    breaks <- handle.r.b.args(r, breaks, W, rmaxdefault=rmaxdefault)
    r <- breaks$r
    rmax <- breaks$max

    if(missing(lambda)) {
      # Estimate density by leave-one-out kernel smoothing
      kerndens <- density(X, ..., sigma=sigma, varcov=varcov, spill=TRUE)
      sigma  <- kerndens$sigma
      varcov <- kerndens$varcov
      raw    <- kerndens$raw
      edg    <- kerndens$edg
      # Leave-one-out in numerator
      SqrtDetSigma <- if(!is.null(sigma)) sigma^2 else sqrt(det(varcov))
      peakvalue <- 1/(2 * pi * SqrtDetSigma)
      lambda <- safelookup(raw, X) - rep(peakvalue, npoints)
      # Edge correction if invoked
      if(!is.null(edg)) 
        lambda <- lambda/safelookup(edg, X)
      # validate
      if(length(lambda) != npoints) 
        stop(paste("Internal error: incorrect number of lambda values",
                   "in leave-one-out method:",
                   "length(lambda) = ", length(lambda),
                   "!=", npoints, "= npoints"))
      if(any(is.na(lambda))) {
        nwrong <- sum(is.na(lambda))
        stop(paste("Internal error:", nwrong, "NA or NaN",
                   ngettext(nwrong, "value", "values"),
                   "generated in leave-one-out method"))
      }
    } else if(is.vector(lambda)) {
      # vector of lambda values provided 
      if(length(lambda) != npoints) 
        stop(paste("The length of", sQuote("lambda"),
                   "should equal the number of data points"))
      if(any(is.na(lambda)))
        stop(paste("Some values of", sQuote("lambda"), "are NA or NaN"))
    } else if(is.im(lambda)) {
      # lambda image provided
      lambda <- safelookup(lambda, X)
    } else
    stop(paste(sQuote("lambda"),
               "should be a vector or a pixel image"))

    
    if(!is.null(lambda2)) {
      if(!is.matrix(lambda2))
        stop("lambda2 should be a matrix")
      if(nrow(lambda2) != ncol(lambda2))
        stop("lambda2 should be a square matrix")
      if(nrow(lambda2) != npoints)
        stop("Dimensions of lambda2 do not match number of data points")
    }
    
    # match corrections
    correction.given <- !missing(correction)
    correction.name <- c("border", "bord.modif", "isotropic", "Ripley", "translate")
    correction.list <- c("border", "bord.modif", "isotropic", "isotropic", "translate")
    correction.id <- pmatch(correction, correction.name)
    if(any(nbg <- is.na(correction.id)))
      stop(paste("unrecognised correction",
                 if(sum(nbg) > 1) "s",
                 ": ",
                 paste(correction[nbg], collapse=", "),
                 sep=""))
    correction <- correction.list[correction.id]
    
    # available selection of edge corrections depends on window
    if(W$type == "mask") {
      iso <- (correction == "isotropic") | (correction == "Ripley")
      if(any(iso)) {
        if(correction.given)
          warning("Isotropic correction not implemented for binary masks")
        correction <- correction[!iso]
      }
    }

    # recommended range of r values
    alim <- c(0, min(rmax, rmaxdefault))
        
  ###########################################
  # Efficient code for border method
  # Usable only if r values are evenly spaced from 0 to rmax
  # Invoked automatically if number of points is large

    can.do.fast <- breaks$even  && missing(lambda2)
    borderonly <- all(correction == "border" | correction == "bord.modif")
    large.n    <- (npoints >= nlarge)
    will.do.fast <- can.do.fast && (borderonly || large.n)
    asked      <- borderonly || (nlarge.given && large.n)

    if(will.do.fast && !asked)
      message(paste("number of data points exceeds",
                    nlarge, "- computing border estimate only"))

    if(asked && !can.do.fast) {
      whynot <-
        if(!(breaks$even)) "r values not evenly spaced" else
        if(!missing(lambda)) "matrix lambda2 was given" else NULL
      warning(paste(c("cannot use efficient code", whynot), sep="; "))
    }
    
    if(will.do.fast) {
    # restrict r values to recommended range, unless specifically requested
      if(!rfixed) 
        r <- seq(0, alim[2], length=length(r))
      return(Kborder.engine(X, max(r), length(r), correction, 1/lambda))
    }

  ###########################################
  # Slower code
  ###########################################
        
        
    # this will be the output data frame
    K <- data.frame(r=r, theo= pi * r^2)
    desc <- c("distance argument r", "theoretical Poisson Kinhom(r)")
    K <- fv(K, "r", substitute(Kinhom(r), NULL),
            "theo", , alim, c("r","Kpois(r)"), desc)
        
    # identify all close pairs
    rmax <- max(r)
    close <- closepairs(X, rmax)
    dIJ <- close$d
    # compute weights for these pairs
    I <- close$i
    J <- close$j
    wI <- 1/lambda[I]
    wIJ <- 
      if(is.null(lambda2))
        1/(lambda[I] * lambda[J])
      else 
        1/lambda2[cbind(I,J)]
    # 
    XI <- X[I]

    # compute edge corrected estimates
    if(any(correction == "border" | correction == "bord.modif")) {
      # border method
      # Compute distances to boundary
      b <- bdist.points(X)
      I <- close$i
      bI <- b[I]
      # apply reduced sample algorithm
      RS <- Kwtsum(dIJ, bI, wIJ, b, w=1/lambda, breaks)
      if(any(correction == "border")) {
        Kb <- RS$ratio
        K <- bind.fv(K, data.frame(border=Kb), "Kbord(r)",
                     "border-corrected estimate of Kinhom(r)",
                     "border")
      }
      if(any(correction == "bord.modif")) {
        Kbm <- RS$numerator/eroded.areas(W, r)
        K <- bind.fv(K, data.frame(bord.modif=Kbm), "Kbord*(r)",
                     "modified border-corrected estimate of Kinhom(r)",
                     "bord.modif")
      }
    }
    if(any(correction == "translate")) {
      # translation correction
      XJ <- X[J]
      edgewt <- edge.Trans(XI, XJ, paired=TRUE)
      allweight <- edgewt * wIJ
      wh <- whist(dIJ, breaks$val, allweight)
      Ktrans <- cumsum(wh)/area
      rmax <- diameter(W)/2
      Ktrans[r >= rmax] <- NA
      K <- bind.fv(K, data.frame(trans=Ktrans), "Ktrans(r)",
                   "translation-correction estimate of Kinhom(r)",
                   "trans")
    }
    if(any(correction == "isotropic" | correction == "Ripley")) {
      # Ripley isotropic correction
      edgewt <- edge.Ripley(XI, matrix(dIJ, ncol=1))
      allweight <- edgewt * wIJ
      wh <- whist(dIJ, breaks$val, allweight)
      Kiso <- cumsum(wh)/area
      rmax <- diameter(W)/2
      Kiso[r >= rmax] <- NA
      K <- bind.fv(K, data.frame(iso=Kiso), "Kiso(r)",
                   "Ripley isotropic correction estimate of Kinhom(r)",
                   "iso")
    }

    # default is to display them all
    attr(K, "fmla") <- . ~ r
    unitname(K) <- unitname(X)
    return(K)
}


Kwtsum <- function(dIJ, bI, wIJ, b, w, breaks) {
  #
  # "internal" routine to compute border-correction estimates of Kinhom
  #
  # dIJ:  vector containing pairwise distances for selected I,J pairs
  # bI:   corresponding vector of boundary distances for I
  # wIJ:  product weight for selected I, J pairs
  #
  # b:    vector of ALL distances to window boundary
  # w:   weights for ALL points
  #
  # breaks : breakpts object
  #

  stopifnot(length(dIJ) == length(bI))
  stopifnot(length(bI) == length(wIJ))
  stopifnot(length(w) == length(b))
  
  # determine which distances d_{ij} were observed without censoring
  uncen <- (dIJ <= bI)
  #
  # histogram of noncensored distances
  nco <- whist(dIJ[uncen], breaks$val, wIJ[uncen])
  # histogram of censoring times for noncensored distances
  ncc <- whist(bI[uncen], breaks$val, wIJ[uncen])
  # histogram of censoring times (yes, this is a different total size)
  cen <- whist(b, breaks$val, w)
  # total weight of censoring times beyond rightmost breakpoint
  uppercen <- sum(w[b > max(breaks$val)])
  # go
  RS <- reduced.sample(nco, cen, ncc, show=TRUE, uppercen=uppercen)
  # extract results
  numerator   <- RS$numerator
  denominator <- RS$denominator
  ratio        <- RS$numerator/RS$denominator
  # check
  if(length(numerator) != breaks$ncells)
    stop("internal error: length(numerator) != breaks$ncells")
  if(length(denominator) != breaks$ncells)
    stop("internal error: length(denom.count) != breaks$ncells")
  return(list(numerator=numerator, denominator=denominator, ratio=ratio))
}
	
