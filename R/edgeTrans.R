#
#        edgeTrans.R
#
#    $Revision: 1.2 $    $Date: 2002/05/06 09:38:12 $
#
#    Translation edge correction weights
#
#  edge.Trans(X)      compute translation correction weights
#                     for each pair of points from point pattern X 
#
#  edge.Trans(X, Y, W)   compute translation correction weights
#                        for all pairs of points X[i] and Y[j]
#                        (i.e. one point from X and one from Y)
#                        in window W
#
#  To estimate the K-function see the idiom in "Kest.S"
#
#######################################################################

edge.Trans <- function(X, Y=X, W=X$window, exact=FALSE, trim=1000) {

  X <- as.ppp(X, W)

  W <- X$window
  x <- X$x
  y <- X$y
  xx <- Y$x
  yy <- Y$y

  # For irregular polygons, exact evaluation is very slow;
  # so use pixel approximation, unless exact=TRUE
  if(W$type == "polygonal" && !exact)
    W <- as.mask(W)

  switch(W$type,
         rectangle={
           # Fast code for this case
           wide <- diff(W$xrange)
           high <- diff(W$yrange)
           DX <- abs(outer(x,xx,"-"))
           DY <- abs(outer(y,yy,"-"))
           weight <- wide * high / ((wide - DX) * (high - DY))
         },
         polygonal={
           # This code is SLOW
           a <- area.owin(W)
           weight <- matrix(, nrow=X$n, ncol=Y$n)
           for(i in seq(X$n)) {
             for(j in seq(Y$n)) {
               shiftvector <- c(x[i],y[i]) - c(xx[j],yy[j])
               Wshift <- shift(W, shiftvector)
               b <- overlap.owin(W, Wshift)
               weight[i,j] <- a/b
             }
           }
         },
         mask={
           # make difference vectors
           DX <- outer(x,xx,"-")
           DY <- outer(y,yy,"-")
           # compute set covariance of window
           g <- setcov(W)
           # evaluate set covariance at these vectors
           gvalues <- lookup.im(g, as.vector(DX), as.vector(DY))
           # reshape
           gvalues <- matrix(gvalues, nrow=X$n, ncol=Y$n)
           weight <- area.owin(W)/gvalues
         }
         )
  weight <- matrix(pmin(weight, trim), nrow=X$n, ncol=Y$n)
  return(weight)
}
