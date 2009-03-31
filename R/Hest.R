#
#    Hest.R
#
#  Contact distribution for a random set
#
#
Hest <- function(X, ...) {
  if(!(is.ppp(X) || is.psp(X) || is.owin(X)))
    stop("X should be an object of class ppp, psp or owin")
# compute distance map
  D <- distmap(X, ...)
  B <- attr(D, "bdry")
# histogram breakpoints 
  dmax <- summary(D)$max
  breaks <- handle.r.b.args(NULL, NULL, as.owin(D), NULL, rmaxdefault=dmax)
#  extract distances and censoring distances
  dist <- as.vector(as.matrix(D))
  bdry <- as.vector(as.matrix(B))
  ok <- !is.na(dist) && !is.na(bdry)
  dist <- dist[ok]
  bdry <- bdry[ok]
# censoring indicators
  d <- (dist <= bdry)
#  observed distances
  o <- pmin(dist, bdry)
# calculate Kaplan-Meier and border corrected estimates
  result <- km.rs(o, bdry, d, breaks)
# add uncorrected estimate - use with care!
  rightmost <- breaks$max
  hh <- hist(dist[dist <= rightmost],breaks=breaks$val,plot=FALSE)$counts
  result$raw <- cumsum(hh)/length(dist)
# convert to data frame  
  result$breaks <- NULL
  result <- as.data.frame(result)
# adjust for zero probability
  if(is.owin(X)) {
    zeroprob <- area.owin(X)/area.owin(as.rectangle(X))
    result[ , c("rs", "km", "raw")] <-
      (result[ , c("rs", "km", "raw")] - zeroprob)/(1-zeroprob)
  }
# convert to class "fv"
  Z <- result[, c("r", "rs", "km", "hazard", "raw")]
  alim <- range(result$r[result$km <= 0.9])
  labl <- c("r", "Fbord(r)", "Fkm(r)",
            "lambda(r)", "Fraw(r)")
  desc <- c("distance argument r",
            "border corrected estimate of F(r)",
            "Kaplan-Meier estimate of F(r)",
            "Kaplan-Meier estimate of hazard function lambda(r)",
            "uncorrected estimate of F(r)")
  Z <- fv(Z, "r", substitute(F(r), NULL), "km", . ~ r, alim, labl, desc)
  attr(Z, "dotnames") <- c("km", "rs")
  unitname(Z) <- unitname(X)
  return(Z)
}

	

