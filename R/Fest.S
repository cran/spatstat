#
#	Fest.S
#
#	S function empty.space()
#	Computes estimates of the empty space function
#
#	$Revision: 4.23 $	$Date: 2009/07/24 05:25:13 $
#
"Fest" <- 	
"empty.space" <-
function(X, ..., eps = NULL, r=NULL, breaks=NULL,
         correction=c("rs", "km", "cs")) {
  verifyclass(X, "ppp")
  
# Intensity estimate
  W <- X$window
  npoints <-X$n
  lambda <- npoints/area.owin(W)
  
# First discretise
  dwin <- as.mask(W, eps)
  dX <- ppp(X$x, X$y, window=dwin, check=FALSE)
#        
# histogram breakpoints 
#
  rmaxdefault <- rmax.rule("F", dwin, lambda)
  breaks <- handle.r.b.args(r, breaks, dwin, eps,
                                  rmaxdefault=rmaxdefault)
  rvals <- breaks$r
  rmax  <- breaks$max
  
# choose correction(s)
  correction.given <- !missing(correction) && !is.null(correction)
  if(is.null(correction))
    correction <- c("rs", "km", "cs")
  correction <- pickoption("correction", correction,
                           c(none="none",
                             border="rs",
                             rs="rs",
                             KM="km",
                             km="km",
                             Kaplan="km",
                             cs="cs",
                             ChiuStoyan="cs",
                             Hanisch="cs",
                             best="km"),
                           multi=TRUE)
  
# initialise fv object
  df <- data.frame(r=rvals, theo=1-exp(-lambda * pi * rvals^2))
  Z <- fv(df, "r", substitute(F(r), NULL), "theo", . ~ r,
          c(0,rmax),
          c("r", "%s[pois](r)"), 
          c("distance argument r", "theoretical Poisson %s"),
          fname="F")
  zeroes <- rep(0, length(rvals))
#
#  compute distances and censoring distances
  if(X$window$type == "rectangle") {
    # original data were in a rectangle
    # output of exactdt() is sufficient
    e <- exactdt(dX)
    dist <- e$d
    bdry <- e$b
  } else {
    # window is irregular..
    # Distance transform & boundary distance for all pixels
    e <- exactdt(dX)
    b <- bdist.pixels(dX$window, coords=FALSE)
    # select only those pixels inside mask
    mm <- dwin$m
    dist <- e$d[mm]
    bdry <- b[mm]
  }
  
# censoring indicators
  d <- (dist <= bdry)
#  observed distances
  o <- pmin(dist, bdry)

### start calculating estimates of F
  
  if("none" %in% correction) {
    #  UNCORRECTED e.d.f. of empty space distances
    if(npoints == 0)
      edf <- zeroes
    else {
      hh <- hist(dist[dist <= rmax],breaks=breaks$val,plot=FALSE)$counts
      edf <- cumsum(hh)/length(dist)
    }
    Z <- bind.fv(Z, data.frame(raw=edf), "%s[raw](r)",
                 "uncorrected estimate of %s", "raw")
  }
  
  if("cs" %in% correction) {
    # Chiu-Stoyan correction
    if(npoints == 0)
      cs <- zeroes
    else {
      #  uncensored distances
      x <- dist[d]
      #  weights
      a <- eroded.areas(W, rvals)
      # calculate Hanisch estimator
      h <- hist(x[x <= rmax], breaks=breaks$val, plot=FALSE)$counts
      H <- cumsum(h/a)
      cs <- H/max(H[is.finite(H)])
    }
    # add to fv object
    Z <- bind.fv(Z, data.frame(cs=cs),
                 "%s[cs](r)", 
                 "Chiu-Stoyan estimate of %s",
                 "cs")
    # modify recommended plot range
    attr(Z, "alim") <- range(rvals[cs <= 0.9])
  }

  if(any(correction %in% c("rs", "km"))) {
    # calculate Kaplan-Meier and border correction (Reduced Sample) estimators
    if(npoints == 0)
      result <- data.frame(rs=zeroes, km=zeroes, hazard=zeroes)
    else {
      result <- km.rs(o, bdry, d, breaks)
      result <- as.data.frame(result[c("rs", "km", "hazard")])
    }
    # add to fv object
    Z <- bind.fv(Z, result,
                 c("%s[bord](r)", "%s[km](r)", "hazard(r)"),
                 c("border corrected estimate of %s",
                   "Kaplan-Meier estimate of %s",
                   "Kaplan-Meier estimate of hazard function lambda(r)"),
                 "km")
    
    # modify recommended plot range
    attr(Z, "alim") <- range(rvals[result$km <= 0.9])
  }
  
  nama <- names(Z)
  attr(Z, "dotnames") <- rev(nama[!(nama %in% c("r", "hazard"))])
  unitname(Z) <- unitname(X)
  return(Z)
}

	
