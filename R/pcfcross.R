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

pcfdot <- function(X, i, ...) {
  stopifnot(is.multitype(X))
  marx <- marks(X)
  if(missing(i)) i <- levels(marx)[1]
  # map i and not-i to two types
  marks(X) <- factor(ifelse(marx == i, "i", "n"), levels=c("i", "n"))
  # extract points of type i and not-i
  splitX <- split(X)
  Xi <- splitX[["i"]]
  Xn <- splitX[["n"]]
  Xall <- unmark(X)
  # estimate intensities
  lambda.i <- summary(Xi)$intensity
  lambda.n <- summary(Xn)$intensity
  lambda.all <- lambda.i + lambda.n
  # compute cross type pcf from i to not-i
  p.in <- pcfcross(X, "i", "n", ...)
  rr <- p.in$r
  # compute pcf of type i points using same parameters
  p.ii   <- do.call("pcf",
                   resolve.defaults(list(Xi),
                                    list(...),
                                    list(r=rr)))
  # add
  p.idot <- eval.fv((p.in * lambda.n + p.ii * lambda.i)/lambda.all)
  #
  attr(p.idot, "ylab") <- substitute(gdot[i](r), list(i=paste(i)))
  return(p.idot)
}

  

