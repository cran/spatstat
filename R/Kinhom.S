#
#	Kinhom.S	Estimation of K function for inhomogeneous patterns
#
#	$Revision: 1.2 $	$Date: 2004/01/12 11:00:45 $
#
#	Kinhom()	compute estimate of K_inhom
#
#       Currently uses border method and slow code...
#                       
#       Reference:
#            Non- and semiparametric estimation of interaction
#	     in inhomogeneous point patterns
#            A.Baddeley, J.Moller, R.Waagepetersen
#            Statistica Neerlandica 54 (2000) 329--350.
#
"Kinhom"<-
  function (X, lambda, r = NULL, breaks = NULL) {
    verifyclass(X, "ppp")
    breaks <- handle.r.b.args(r, breaks, X$window)
    r <- breaks$r
    if (!is.vector(lambda))
        stop("The argument \'lambda\' should be a vector")
    if (length(lambda) != X$n) 
        stop("The length of the vector \'lambda\' should equal the number of data points")
    d <- pairdist(X$x, X$y)
    diag(d) <- Inf
    b <- bdist.points(X)
    w <- 1/lambda
    RS <- Kwtsum(d, b, w, w, breaks)
    K <- RS$numerator/RS$denominator
    w <- X$window
    alim <- c(0, min(diff(w$xrange), diff(w$yrange))/4)
    out <- fv(data.frame(r = r, theo = pi*r**2, border=K),
              "r", "Kinhom(r)", "border",
              cbind(border, theo) ~ r, alim,
              c("r", "Kpois(r)", "Kinhom(r)"),
              c("distance argument r",
                "theoretical Poisson K(r)",
                "border-corrected estimate of Kinhom(r)"))
    return(out)
}


	
Kwtsum <- function(d, b, wi, wj, breaks) {
  #
  # "internal" routine to compute border-correction estimate of Kinhom
  #
  # d : matrix of pairwise distances
  #                  (to exclude diagonal entries, set diag(d) = Inf)
  # b : column vector of distances to window boundary
  # wi : vector of weights on rows of d (estimate of 1/lambda_i)
  # wj : vector of weights on columns of d (estimate of 1/lambda_j)
  # 
  # breaks : breakpts object
  #

  if(length(wi) != nrow(d))
    stop("length of \'wi\' doesn\'t match number of rows of \'d\'")
  
  if(length(wj) != ncol(d))
    stop("length of \'wj\' doesn\'t match number of columns of \'d\'")
  
  r <- breaks$r
  nr <- length(r)
  numerator <- numeric(nr)
  denominator <- numeric(nr)

  wiwj <- outer(wi, wj, "*")
       
  for(k in 1:nr) {
    close <- (d <= r[k]) # assuming d[i,j] = Inf when X[i] == X[j]
    numi <- matrowsum(close * wiwj)
    bok <- (b > r[k]) 
    numerator[k] <- sum(numi[bok])
    denominator[k] <- sum(wi[bok])
  }
	
  return(list(RS=numerator/denominator,
              numerator=numerator,
              denominator=denominator))
}
