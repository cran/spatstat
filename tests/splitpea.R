#
#  tests/splitpea.R
#
#  Check behaviour of split.ppp etc
#
#  Thanks to Marcelino de la Cruz
#
#  $Revision: 1.6 $  $Date: 2011/08/08 06:34:09 $
#

require(spatstat)
if(require(gpclib)) spatstat.options(gpclib=TRUE)

W <- square(8)
X <- ppp(c(2.98, 4.58, 7.27, 1.61, 7.19),
         c(7.56, 5.29, 5.03, 0.49, 1.65),
         window=W)
Z <- quadrats(W, 4, 4)
Yall <- split(X, Z, drop=FALSE)
Ydrop <- split(X, Z, drop=TRUE)

P <- Yall[[1]]
if(!all(inside.owin(P$x, P$y, P$window)))
  stop("Black hole detected when drop=FALSE")
P <- Ydrop[[1]]
if(!all(inside.owin(P$x, P$y, P$window)))
  stop("Black hole detected when drop=TRUE")

Ydrop[[1]] <- P[1]
split(X, Z, drop=TRUE) <- Ydrop

# test NA handling
Zbad <- quadrats(square(4), 2, 2)
Ybdrop <- split(X, Zbad, drop=TRUE)
Yball  <- split(X, Zbad, drop=FALSE)

# From Marcelino
set.seed(1)
W<- square(10) # the big window
puntos<- rpoispp(0.5, win=W)
data(letterR)
r00 <- letterR
r05 <- shift(letterR,c(0,5))
r50 <- shift(letterR,c(5,0))
r55 <- shift(letterR,c(5,5))
tessr4 <- tess(tiles=list(r00, r05,r50,r55))
puntosr4 <- split(puntos, tessr4, drop=TRUE)
split(puntos, tessr4, drop=TRUE) <- puntosr4
