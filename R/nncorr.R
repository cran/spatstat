#
# nncorr.R
#
# $Revision: 1.2 $  $Date: 2006/04/18 09:19:41 $
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
  m  <- X$marks
  Efmm <- mean(outer(m, m, f, ...))
  # border method
  nn <- nnwhich(X)
  ok <- (nndist(X) <= bdist.points(X))
  if(!any(ok))
    stop("Insufficient data")
  Y <- X[nn[ok]]
  X <- X[ok]
  Efmk <- mean(f(X$marks, Y$marks, ...))
  #
  answer <- c(unnormalised=Efmk,
              normalised=Efmk/Efmm)
  if(!is.factor(X$marks) && f.is.default) {
    classic <- cor(X$marks, Y$marks, use=use, method=method)
    answer <- c(answer, correlation=classic)
  }
  return(answer)
}
  
