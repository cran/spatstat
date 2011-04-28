#
#  dist2dpath.R
#
#   $Revision: 1.4 $    $Date: 2011/04/18 08:44:32 $
#
#       dist2dpath    compute shortest path distances
#

dist2dpath <- function(dist, method="C") {
  # given a matrix of distances between adjacent vertices
  # (value = Inf if not adjacent)
  # compute the matrix of shortest path distances
  stopifnot(is.matrix(dist) && isSymmetric(dist))
  stopifnot(all(diag(dist) == 0))
  #
  n <- nrow(dist)
  cols <- col(dist)
  #
  switch(method,
         interpreted={
           dpathnew <- dpath <- dist
           changed <- TRUE
           while(changed) {
             for(j in 1:n) 
               dpathnew[,j] <- apply(dpath + dist[j,][cols], 1, min)
             changed <- any(dpathnew != dpath)
             dpath <- dpathnew
           }
         },
         C={
           adj <- is.finite(dist)
           diag(adj) <- TRUE
           d <- dist
           d[!adj] <- -1
           z <- .C("dist2dpath",
                   nv=as.integer(n),
                   d=as.double(d),
                   adj=as.integer(adj),
                   dpath=as.double(numeric(n*n)),
                   niter=as.integer(integer(1)),
                   PACKAGE="spatstat")
           if(z$niter >= n)
             warning("C algorithm did not converge")
           dpath <- matrix(z$dpath, n, n)
           # value=-1 implies unreachable
           dpath[dpath < 0] <- Inf
         },
         stop(paste("Unrecognised method", sQuote(method))))
  return(dpath)
}
