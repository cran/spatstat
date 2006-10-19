#
# nncorr.R
#
# $Revision: 1.3 $  $Date: 2006/10/10 04:22:48 $
#

nncorr <- function(X, f = function(m1,m2) { m1 * m2}, ...,
                   use = "all.obs",
                   method = c("pearson", "kendall", "spearman")) {
  verifyclass(X, "ppp")
  if(!is.marked(X))
    stop("X does not have marks")
  if(missing(method) || is.null(method))
    method <- "pearson"
  f.is.default <- missing(f)
  stopifnot(is.function(f))
  # denominator
  m  <- marks(X, dfok=FALSE)
  Efmm <- mean(outer(m, m, f, ...))
  # border method
  nn <- nnwhich(X)
  ok <- (nndist(X) <= bdist.points(X))
  if(!any(ok))
    stop("Insufficient data")
  Y <- X[nn[ok]]
  X <- X[ok]
  Efmk <- mean(f(marks(X), marks(Y), ...))
  #
  answer <- c(unnormalised=Efmk,
              normalised=Efmk/Efmm)
  if(!is.multitype(X) && f.is.default) {
    classic <- cor(marks(X), marks(Y), use=use, method=method)
    answer <- c(answer, correlation=classic)
  }
  return(answer)
}
  
