#	Jest.S
#
#	Usual invocation to compute J function
#	if F and G are not required 
#
#	$Revision: 4.12 $	$Date: 2009/04/07 07:01:28 $
#
#
#
"Jest" <-
function(X, eps=NULL, r=NULL, breaks=NULL) {
  X <- as.ppp(X)
  W<- X$window
  rmaxdefault <- rmax.rule("J", W)
  brks <- handle.r.b.args(r, breaks, W, rmaxdefault=rmaxdefault)$val
  FF <- Fest(X, eps, breaks=brks)
  G <- Gest(X, breaks=brks)
  ratio <- function(a, b, c) {
    result <- a/b
    result[ b == 0 ] <- c
    result
  }
  Jrs <- ratio(1-G$rs, 1-FF$rs, NA)
  Jkm <- ratio(1-G$km, 1-FF$km, NA)
  Jun <- ratio(1-G$raw, 1-FF$raw, NA)
  theo <- rep(1, length(FF$r))
  
  rslt <- data.frame(r=FF$r, theo=theo, un=Jun, rs=Jrs, km=Jkm)
# convert to class "fv"
  alim <- range(rslt$r[FF$km <= 0.9])
  labl <- c("r", "%spois(r)", "%sun(r)", "%sbord(r)", "%skm(r)")
  desc <- c("distance argument r",
            "theoretical Poisson %s",
            "uncorrected estimate of %s",
            "border corrected estimate of %s",
            "Kaplan-Meier estimate of %s")
  Z <- fv(rslt, "r", substitute(J(r), NULL),
          "km", . ~ r, alim, labl, desc, fname="J")

# add more info        
  attr(Z, "F") <- FF
  attr(Z, "G") <- G

  unitname(Z) <- unitname(X)
  return(Z)
}

