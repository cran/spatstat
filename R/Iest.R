#	Iest.R
#
#	I function
#
#	$Revision: 1.4 $	$Date: 2006/10/18 06:06:02 $
#
#
#
Iest <- function(X, eps=NULL, r = NULL, breaks = NULL) {

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
  Jdotdot <- Jest(unmark(X))
  rvals <- Jdotdot$r
  
  # J function of subpattern of each type i
  Jii <- lapply(Y, Jest, r=rvals)
  nrvals <- unlist(lapply(Jii, function(x) { length(x$r) }))
  if(length(unique(nrvals)) != 1 || nrvals[1] != length(rvals))
    stop("Internal error: J function objects have different lengths")

  # Estimates of each type
  extractit <- function(Z, what) { Z[[what]] }
  extract <- function(Zlist, what) { unlist(lapply(Zlist, extractit, what=what)) }
  Jrs <- matrix(extract(Jii, "rs"), nrow=ntypes, byrow=TRUE)
  Jkm <- matrix(extract(Jii, "km"), nrow=ntypes, byrow=TRUE)
  Jun <- matrix(extract(Jii, "un"), nrow=ntypes, byrow=TRUE)

  # Calculate
  Irs <- apply(fi * Jrs, 2, sum) - Jdotdot$rs
  Ikm <- apply(fi * Jkm, 2, sum) - Jdotdot$km
  Iun <- apply(fi * Jun, 2, sum) - Jdotdot$un

  rslt <- data.frame(r=rvals, theo=rep(0, length(rvals)),
                     rs=Irs, km=Ikm, un=Iun)
  alim <- attr(Jdotdot, "alim")
  labl <- c("r", "Ipois(r)", "Ibord(r)", "Ikm(r)", "Iun(r)")
  desc <- c("distance argument r",
            "theoretical Poisson I(r) = 0",
            "border corrected estimate of I(r)",
            "Kaplan-Meier estimate of I(r)",
            "uncorrected estimate of I(r)")
  Z <- fv(rslt, "r", "I(r)", "km", cbind(km, rs, un, theo) ~ r,
          alim, labl, desc)
  units(Z) <- units(X)
  return(Z)
}


