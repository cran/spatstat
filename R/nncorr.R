#
# nncorr.R
#
# $Revision: 1.6 $  $Date: 2010/03/09 04:25:51 $
#

nnmean <- function(X) {
  stopifnot(is.ppp(X) && is.marked(X))
  m <- as.data.frame(marks(X))
  nw <- nnwhich(X)
  meannum <- function(x) { if(!is.numeric(x)) NA else mean(x) }
  unlist(lapply(m[nw,], meannum))/unlist(lapply(m, meannum))
}

nnvario <- function(X) {
  m <- onecolumn(marks(X))
  f <- function(m1,m2) { ((m1-m2)^2)/2 }
  nncorr(X, f, denominator=var(m))
}

nncorr <- function(X, f = function(m1,m2) { m1 * m2}, ...,
                   use = "all.obs",
                   method = c("pearson", "kendall", "spearman")) {
  verifyclass(X, "ppp")
  if(missing(method) || is.null(method))
    method <- "pearson"
  f.is.default <- missing(f) || is.null(f)
  stopifnot(is.function(f))
  m  <- onecolumn(marks(X))
  # denominator
  Efmm <-
    if(f.is.default)
      mean(m)^2
    else if(!is.null(denominator <- list(...)$denominator))
      denominator
    else 
      mean(outer(m, m, f, ...))
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
  
