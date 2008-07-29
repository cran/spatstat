#
#
#  pcfcross.R
#
#  kernel estimation of cross-type pair correlation function
#
#  Currently computed by differencing pcf
#
pcfcross <- function(X, i, j, ...) {
  stopifnot(is.multitype(X))
  if(missing(i)) i <- levels(marks(X))[1]
  if(missing(j)) j <- levels(marks(X))[2]
  # extract points of types i and j
  Xsplit <- split(X)
  Xi <- Xsplit[[i]]
  Xj <- Xsplit[[j]]
  if(i == j) 
    return(pcf(Xi, ...))
  Xall <- superimpose(Xi, Xj, W=X$window)
  # estimate intensities
  lambda.i <- summary(Xi)$intensity
  lambda.j <- summary(Xj)$intensity
  lambda.all <- lambda.i + lambda.j
  # kernel estimates of unmarked pcf's
  p.all <- pcf(Xall, ...)
  rr <- p.all$r
  p.ii   <- do.call("pcf",
                   resolve.defaults(list(Xi),
                                    list(...),
                                    list(r=rr)))
  p.jj   <- do.call("pcf",
                   resolve.defaults(list(Xj),
                                    list(...),
                                    list(r=rr)))
  # differencing
  p.ij <- eval.fv((p.all * lambda.all^2
                   - p.ii * lambda.i^2
                   - p.jj * lambda.j^2)/(2 * lambda.i * lambda.j))
  #
  attr(p.ij, "ylab") <- substitute(g[i,j](r), list(i=paste(i),j=paste(j)))
  return(p.ij)
}


  
  

