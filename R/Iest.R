#	Iest.R
#
#	I function
#
#	$Revision: 1.10 $	$Date: 2009/06/21 01:43:07 $
#
#
#
Iest <- function(X, ..., eps=NULL, r = NULL, breaks = NULL, correction=NULL) {

  X <- as.ppp(X)
  if(!is.multitype(X))
    stop("Only applicable to multitype point patterns")
  marx <- marks(X, dfok=FALSE)
  ntypes <- length(levels(marx))

  Y <- unmark(split(X))
  
  # relative proportions 
  ni <- unlist(lapply(Y, function(Z) { Z$n }))
  fi <- ni/sum(ni)

  # J function of pattern regardless of type
  Jdotdot <- Jest(unmark(X), correction=correction)
  rvals <- Jdotdot$r
  
  # J function of subpattern of each type i
  Jii <- lapply(Y, Jest, r=rvals, correction=correction)
  nrvals <- unlist(lapply(Jii, function(x) { length(x$r) }))
  if(length(unique(nrvals)) != 1 || nrvals[1] != length(rvals))
    stop("Internal error: J function objects have different lengths")

  # initialise fv object
  alim <- attr(Jdotdot, "alim")
  Z <- fv(data.frame(r=rvals, theo=0),
          "r", substitute(I(r), NULL), "theo",
          . ~ r, alim,
          c("r", "%spois(r)"),
          c("distance argument r", "theoretical Poisson %s"),
          fname="I")
  
  # Estimates of each type
  extractit <- function(Z, what) { Z[[what]] }
  extract <- function(Zlist, what) { unlist(lapply(Zlist, extractit, what=what)) }
  namii <- unlist(lapply(Jii, names))
  namdd <- names(Jdotdot)
  bothnames <- namii[namii %in% namdd]
  
  if("un" %in% bothnames) {
    Jun <- matrix(extract(Jii, "un"), nrow=ntypes, byrow=TRUE)
    Iun <- apply(fi * Jun, 2, sum) - Jdotdot$un
    Z <- bind.fv(Z, data.frame(un=Iun), "%sun(r)",
                 "uncorrected estimate of %s", "un")
  }
  if("rs" %in% bothnames) {
    Jrs <- matrix(extract(Jii, "rs"), nrow=ntypes, byrow=TRUE)
    Irs <- apply(fi * Jrs, 2, sum) - Jdotdot$rs    
    Z <- bind.fv(Z, data.frame(rs=Irs), "%srs(r)",
                 "border corrected estimate of %s", "rs")
  }
  if("han" %in% bothnames) {
    Jhan <- matrix(extract(Jii, "han"), nrow=ntypes, byrow=TRUE)
    Ihan <- apply(fi * Jhan, 2, sum) - Jdotdot$han
    Z <- bind.fv(Z, data.frame(han=Ihan), "%shan(r)",
                 "Hanisch-style estimate of %s", "han")
  }
  if("km" %in% bothnames) {
    Jkm <- matrix(extract(Jii, "km"), nrow=ntypes, byrow=TRUE)
    Ikm <- apply(fi * Jkm, 2, sum) - Jdotdot$km
    Z <- bind.fv(Z, data.frame(km=Ikm), "%skm(r)",
                 "Kaplan-Meier estimate of %s", "km")
  }
  unitname(Z) <- unitname(X)
  return(Z)
}


