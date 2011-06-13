#
# linearK
#
# $Revision: 1.3 $ $Date: 2011/06/11 01:43:36 $
#
# K function for point pattern on linear network
#
#
linearK <- function(X, r=NULL, ...) {
  stopifnot(inherits(X, "lpp"))
  # extract info about pattern
  sX <- summary(X)
  np <- sX$npoints
  lengthL <- sX$totlength
  # compute K
  denom <- np * (np - 1)/lengthL
  # extract linear network
  L <- X$domain
  # extract points
  Y <- as.ppp(X)
  W <- Y$window
  # determine r values
  rmaxdefault <- 0.98 * circumradius(L)
  breaks <- handle.r.b.args(r, NULL, W, rmaxdefault=rmaxdefault)
  r <- breaks$r
  rmax <- breaks$max
  #
  ylab <- substitute(K[net](r), NULL)
  fname <- "K[net]"
  #
  if(np < 2) {
    # no pairs to count: return zero function
    zeroes <- rep(0, length(r))
    df <- data.frame(r = r, est = zeroes)
    K <- fv(df, "r", ylab, 
            "est", . ~ r, c(0, rmax),
            c("r", "%s(r)"),
            c("distance argument r", "estimated %s"),
            fname = fname)
  } else {
    # compute pairwise distances  
    D <- pairdist(X)
    # compile into Okabe-Yamada K function 
    K <- compileK(D, r, denom=denom)
    # set appropriate y axis label
    K <- rebadge.fv(K, new.ylab=ylab, new.fname=fname)
  }
  unitname(K) <- unitname(X)
  return(K)
}

